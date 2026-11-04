"""Lazy, process-local access to immutable release files.

No data files or native handles are opened at import time. Restart workers when
publishing a release update; caches deliberately assume immutable input files.
"""
from collections import OrderedDict
from contextlib import closing
from functools import cached_property
from io import BytesIO
from pathlib import Path
import json
import sqlite3
from threading import RLock
import zipfile
import zlib

import numpy as np
import pandas as pd

try:
    from .association import AssociationDataError, CovariateModel
    from .catalog import query_catalog, source_signature
except ImportError:  # gunicorn --chdir api app:api
    from association import AssociationDataError, CovariateModel
    from catalog import query_catalog, source_signature


SCALES = ("tpm", "normalized")
MODALITIES = ("expression", "isoforms", "splicing", "alt_TSS", "alt_polyA", "stability")


def phenotype_index_path(root, tissue, modality, scale="normalized"):
    return root / "api/phenotypes" / tissue / f"{modality}.{scale}.sqlite3"


def expression_path(root, tissue, scale):
    if scale not in SCALES:
        raise ValueError(f"Unknown expression scale: {scale}")
    if scale == "tpm":
        path = root / "expr" / f"expr.tpm.{tissue}.v4_rn8.bed.gz"
        if path.is_file():
            return path
    normalization = "unnorm" if scale == "tpm" else "norm"
    return root / "phenos" / f"phenos.{tissue}.expression.{normalization}.v4_rn8.bed.gz"


class locked_cached_property(cached_property):
    """cached_property without duplicate cold loads under threaded workers."""
    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        with instance._lock:
            return super().__get__(instance, owner)


class DataStore:
    def __init__(self, root, expression_cache_bytes=128 * 1024**2, archive_cache_size=4):
        if expression_cache_bytes < 0 or archive_cache_size < 1:
            raise ValueError("Cache bytes must be nonnegative and archive cache size must be positive.")
        self.root = Path(root).resolve()
        self.expression_cache_bytes = expression_cache_bytes
        self.archive_cache_size = archive_cache_size
        self._expression = OrderedDict()
        self._expression_bytes = 0
        self._archives = OrderedDict()
        self._covariates = OrderedDict()
        self._lock = RLock()

    @locked_cached_property
    def tissue_info(self):
        return pd.read_csv(self.root / "tissue_info.v4.tsv", sep="\t").to_dict("records")

    @locked_cached_property
    def tissues(self):
        return [row["tissueSiteDetailId"] for row in self.tissue_info]

    @locked_cached_property
    def datasets(self):
        return {row["tissueSiteDetailId"]: row["dataset"] for row in self.tissue_info}

    @locked_cached_property
    def genes(self):
        return pd.read_csv(self.root / "gene.v4.tsv", sep="\t", index_col="geneId",
                           dtype={"chromosome": str}).fillna("")

    @locked_cached_property
    def gene_lookup(self):
        lookup = {}
        for gene in self.genes.index:
            key = gene.casefold()
            lookup[key] = gene if key not in lookup else None
        return lookup

    def resolve_gene(self, gene):
        return gene if gene in self.genes.index else self.gene_lookup.get(gene.casefold())

    @locked_cached_property
    def median_expression(self):
        return pd.read_csv(self.root / "expr/medianGeneExpression.v4.tsv.gz",
                           sep="\t", index_col="geneId")

    @locked_cached_property
    def top_expression(self):
        return pd.read_csv(self.root / "expr/topExpressedGene.v4.tsv", sep="\t")

    def catalog_records(self, collection, gene, tissues=None, modality=None, limit=None, offset=0):
        return query_catalog(self.root, collection, gene, tissues, modality, limit, offset)

    def expression(self, tissue, gene, scale="tpm"):
        """Return an independent Series; never expose a mutable cached frame."""
        path = expression_path(self.root, tissue, scale)
        index = phenotype_index_path(self.root, tissue, "expression", scale)
        if index.is_file():
            with closing(sqlite3.connect(index.as_uri() + "?mode=ro", uri=True)) as conn:
                meta = dict(conn.execute("SELECT key, value FROM metadata"))
                if (meta.get("phenotype_format") == "phenotypes-v1"
                        and meta.get("tissue") == tissue and meta.get("modality") == "expression"
                        and meta.get("scale") == scale
                        and json.loads(meta["source"]) == source_signature(path)):
                    row = conn.execute("SELECT data FROM phenotype WHERE phenotype_id = ?", (gene,)).fetchone()
                    if row is None:
                        return None
                    samples = json.loads(meta["samples"])
                    values = np.frombuffer(zlib.decompress(row[0]), dtype="<f8").copy()
                    return pd.Series(values, index=samples, name=gene)
        # Serialize cold fallback loads to avoid duplicate large matrices.
        with self._lock:
            key = (tissue, scale)
            if key in self._expression:
                frame, size = self._expression.pop(key)
                self._expression[key] = (frame, size)
            else:
                frame = pd.read_csv(path, sep="\t", index_col="phenotype_id",
                                    usecols=lambda col: col not in {"#chr", "start", "end"})
                size = int(frame.memory_usage(deep=True).sum())
                while self._expression and self._expression_bytes + size > self.expression_cache_bytes:
                    _, (_, removed) = self._expression.popitem(last=False)
                    self._expression_bytes -= removed
                if size <= self.expression_cache_bytes:
                    self._expression[key] = (frame, size)
                    self._expression_bytes += size
            return frame.loc[gene].astype(float).copy() if gene in frame.index else None

    def covariate_model(self, tissue, modality):
        with self._lock:
            key = (tissue, modality)
            if key in self._covariates:
                model = self._covariates.pop(key)
            else:
                path = self.root / "covar" / f"covar.{tissue}.{modality}.v4.tsv"
                # Check the raw header before pandas can rename duplicate IDs.
                with path.open() as handle:
                    samples = handle.readline().rstrip("\n\r").split("\t")[1:]
                if len(samples) != len(set(samples)) or any(not sample for sample in samples):
                    raise AssociationDataError("Duplicate or empty covariate sample IDs.")
                model = CovariateModel.from_frame(pd.read_csv(path, sep="\t", index_col=0).T)
            self._covariates[key] = model
            while len(self._covariates) > 8:
                self._covariates.popitem(last=False)
            return model

    def nominal_associations(self, tissue, modality, phenotypes, genotypes):
        """Infer with released modality-specific covariates, never plot-scale TPM."""
        model = self.covariate_model(tissue, modality)
        if (not phenotypes.columns.is_unique or set(phenotypes.columns) != set(model.samples)
                or not set(model.samples).issubset(genotypes.index)):
            raise AssociationDataError("Phenotype, covariate and VCF samples do not match.")
        results = model.test(phenotypes.loc[:, list(model.samples)].to_numpy(),
                             genotypes.reindex(model.samples).to_numpy(dtype=float))
        # v4's PLINK mapping inputs use --maf 0.01 --mac 2, before tensorQTL's
        # mean imputation. Check the tissue cohort, not a genotype plot group.
        cohort_path = self.root / "rat_ids" / f"rat_ids.{tissue}.v4.txt"
        cohort = cohort_path.read_text().split() if cohort_path.is_file() else model.samples
        if not set(cohort).issubset(genotypes.index):
            raise AssociationDataError("Mapping cohort samples are absent from the VCF.")
        called = genotypes.reindex(cohort).to_numpy(dtype=float)
        called = called[np.isfinite(called)]
        mac = min(called.sum(), 2 * len(called) - called.sum()) if len(called) else 0
        if len(called) and (mac < 2 or mac / (2 * len(called)) < .01):
            results = [{**result, "pValue": None, "slope": None, "slopeSE": None,
                        "associationStatus": "variant_below_mapping_filter"} for result in results]
        return results

    def _phenotype_rows(self, tissue, modality, gene, with_data=True):
        """Fetch only this gene's normalized vectors, loci and group threshold.

        Multimodal queries require prepared indexes: never fall back to loading
        a full modality matrix into a web worker.
        """
        path = phenotype_index_path(self.root, tissue, modality)
        if not path.is_file():
            raise AssociationDataError("Build normalized phenotype indexes with api.prepare_phenotypes.")
        source = self.root / "phenos" / f"phenos.{tissue}.{modality}.norm.v4_rn8.bed.gz"
        thresholds = self.root / "xqtl" / f"top_assoc.{modality}.v4_rn8.tsv"
        with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as conn:
            meta = dict(conn.execute("SELECT key, value FROM metadata"))
            if (meta.get("phenotype_format") != "phenotypes-v1"
                    or meta.get("modality") != modality
                    or meta.get("tissue") != tissue
                    or meta.get("scale") != "normalized"
                    or json.loads(meta["source"]) != source_signature(source)
                    or json.loads(meta["threshold_source"]) != source_signature(thresholds)):
                raise AssociationDataError("Missing or stale phenotype index; run api.prepare_phenotypes.")
            samples = json.loads(meta["samples"])
            columns = "phenotype_id, chromosome, start, end" + (", data" if with_data else "")
            rows = conn.execute(f"SELECT {columns} FROM phenotype WHERE gene_id = ? ORDER BY phenotype_id", (gene,)).fetchall()
            threshold = conn.execute("SELECT value FROM threshold WHERE gene_id = ?", (gene,)).fetchone()
        return rows, samples, threshold[0] if threshold else None

    def phenotype_metadata(self, tissue, modality, gene):
        """List measured phenotypes without reading or decompressing sample blobs."""
        rows, samples, threshold = self._phenotype_rows(tissue, modality, gene, with_data=False)
        return [{"phenotypeId": phenotype, "geneId": gene, "modality": modality,
                 "tissueSiteDetailId": tissue, "chromosome": chrom,
                 "start": start + 1, "end": end, "sampleCount": len(samples),
                 "pValueThreshold": threshold} for phenotype, chrom, start, end in rows]

    def phenotypes(self, tissue, modality, gene):
        rows, samples, threshold = self._phenotype_rows(tissue, modality, gene)
        vectors, loci = {}, {}
        for phenotype, chrom, start, end, blob in rows:
            vectors[phenotype] = np.frombuffer(zlib.decompress(blob), dtype="<f8").copy()
            # tensorQTL converts BED start to 1-based coordinates.
            loci[phenotype] = (chrom, start + 1, end)
        return pd.DataFrame.from_dict(vectors, orient="index", columns=samples), loci, threshold

    def _archive(self, path):
        # Called only while holding _lock. Bound open handles and ZIP directories.
        if path in self._archives:
            archive = self._archives.pop(path)
        else:
            archive = zipfile.ZipFile(path)
        self._archives[path] = archive
        while len(self._archives) > self.archive_cache_size:
            _, old = self._archives.popitem(last=False)
            old.close()
        return archive

    def significant_associations(self, gene):
        with self._lock:
            archive = self._archive(self.root / "eqtl/all_gene_esnps.v4.zip")
            try:
                contents = archive.read(f"all_gene_esnps.v4/{gene}.txt")
            except KeyError:
                return pd.DataFrame()
        frame = pd.read_csv(BytesIO(contents), sep="\t", dtype={"chromosome": str})
        frame["geneId"] = gene
        # GTEx's legacy 'nes' is the pipeline's regression slope, not log2 aFC.
        frame["slope"] = frame["nes"]
        return frame

    def close(self):
        with self._lock:
            for archive in self._archives.values():
                archive.close()
            self._archives.clear()
            self._expression.clear()
            self._expression_bytes = 0
            self._covariates.clear()
