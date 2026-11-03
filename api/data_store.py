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


SCALES = ("tpm", "normalized")
ASSOCIATION_COLUMNS = {
    "tissue": "tissueSiteDetailId", "gene_id": "geneId",
    "phenotype_id": "phenotypeId", "variant_id": "variantId",
    "pval_beta": "pValueBeta", "log2_aFC": "log2aFC", "chrom": "chromosome",
    "pos": "pos", "af": "altAlleleFrequency", "tss_distance": "tssDistance",
    "pval_nominal": "pValue", "slope": "slope", "slope_se": "slopeSE",
    "qval": "qValue", "rank": "rank", "modality": "modality",
    "ref": "ref", "alt": "alt",
}


def expression_path(root, tissue, scale):
    if scale not in SCALES:
        raise ValueError(f"Unknown expression scale: {scale}")
    if scale == "tpm":
        path = root / "expr" / f"expr.tpm.{tissue}.v4_rn8.bed.gz"
        if path.is_file():
            return path
    normalization = "unnorm" if scale == "tpm" else "norm"
    return root / "phenos" / f"phenos.{tissue}.expression.{normalization}.v4_rn8.bed.gz"


def source_signature(path):
    stat = path.stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


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

    @locked_cached_property
    def exons(self):
        return pd.read_csv(self.root / "exon.v4.tsv", sep="\t", dtype={"chromosome": str}).set_index("geneId", drop=False)

    def _associations(self, filename):
        frame = pd.read_csv(self.root / filename, sep="\t",
                            usecols=lambda col: col in ASSOCIATION_COLUMNS,
                            dtype={"chrom": str})
        return frame.rename(columns=ASSOCIATION_COLUMNS).set_index("geneId", drop=False)

    @locked_cached_property
    def eqtls(self):
        return self._associations("eqtl/eqtls_indep.v4_rn8.tsv")

    @locked_cached_property
    def xqtls(self):
        return self._associations("xqtl/xqtls_indep.cross_modality.v4_rn8.tsv")

    @locked_cached_property
    def thresholds(self):
        return pd.read_csv(self.root / "eqtl/top_assoc.v4_rn8.tsv", sep="\t",
                           usecols=["tissue", "gene_id", "pval_nominal_threshold"],
                           index_col=["tissue", "gene_id"])["pval_nominal_threshold"]

    @staticmethod
    def by_gene(frame, gene):
        return frame.loc[[gene]] if gene in frame.index else frame.iloc[:0]

    def expression(self, tissue, gene, scale="tpm"):
        """Return an independent Series; never expose a mutable cached frame."""
        path = expression_path(self.root, tissue, scale)
        index = self.root / "api/expression" / f"{tissue}.{scale}.sqlite3"
        if index.is_file():
            with closing(sqlite3.connect(index.as_uri() + "?mode=ro", uri=True)) as conn:
                meta = dict(conn.execute("SELECT key, value FROM metadata"))
                if json.loads(meta["source"]) == source_signature(path):
                    row = conn.execute("SELECT data FROM expression WHERE gene_id = ?", (gene,)).fetchone()
                    if row is None:
                        return None
                    samples = json.loads(meta["samples"])
                    values = np.frombuffer(zlib.decompress(row[0]), dtype="<f8").copy()
                    return pd.Series(values, index=samples, name=gene)
        # Older releases work without preprocessing. Serialize cold loads to avoid
        # concurrent requests loading the same large matrix more than once.
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

    def cis_pvalue(self, tissue, gene, variant):
        with self._lock:
            archive = self._archive(self.root / "cis_pvals" / f"{tissue}.v4.zip")
            try:
                with archive.open(f"{tissue}.v4/{gene}.txt") as member:
                    next(member, None)
                    for line in member:
                        variant_id, value = line.decode().rstrip().split("\t")
                        if variant_id == variant:
                            return float(value)
            except KeyError:
                pass
        return None

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
