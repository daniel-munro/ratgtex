"""RatGTEx v4 read API.

Configuration and data preparation: see api/README.md.
"""
from collections import OrderedDict
import gzip
import hashlib
import math
import os
from pathlib import Path
from threading import RLock
from time import monotonic

from flask import Flask, jsonify, request
from flask.json.provider import DefaultJSONProvider
from flask_cors import CORS
import numpy as np
import pandas as pd
import pysam
from scipy.cluster.hierarchy import linkage, to_tree
from werkzeug.exceptions import HTTPException

if __package__:
    from .data_store import DataStore, SCALES, MODALITIES
    from .association import AssociationDataError
else:  # Existing gunicorn --chdir api app:api and python api/app.py commands.
    from data_store import DataStore, SCALES, MODALITIES
    from association import AssociationDataError


class APIError(Exception):
    def __init__(self, message, status=400):
        self.message, self.status = message, status


def clean_json(value):
    """JSON null for missing/nonfinite values, never nonstandard NaN/Infinity."""
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


class ReleaseJSONProvider(DefaultJSONProvider):
    sort_keys = False

    def dumps(self, obj, **kwargs):
        kwargs["allow_nan"] = False
        return super().dumps(clean_json(obj), **kwargs)


def genotype(gt):
    """Alternate-allele dosage for a fully called diploid biallelic genotype."""
    if gt is None or len(gt) != 2 or any(allele not in (0, 1) for allele in gt):
        return None
    return sum(gt)


def variant_record(variant_id, vcf):
    """Match the exact 1-based VCF position, not an overlapping indel."""
    try:
        chrom, position = variant_id.split(":")
        position = int(position)
        if position < 1:
            return None
        records = [record for record in vcf.fetch(chrom, position - 1, position)
                   if record.pos == position]
        return records[0] if len(records) == 1 else None
    except (ValueError, KeyError, IndexError, TypeError):
        return None


def get_newick(node, newick, parentdist, leaf_names):
    """Serialize a scipy hierarchy using the portal's existing Newick format."""
    if node.is_leaf():
        return "%s:%g%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    newick = "):%g%s" % (parentdist - node.dist, newick) if newick else ");"
    newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
    newick = get_newick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
    return "(%s" % newick


def row_tree(frame):
    if len(frame) < 2 or frame.shape[1] == 0 or not np.isfinite(frame.to_numpy()).all():
        return "Not enough data"
    tree = to_tree(linkage(frame, method="average", optimal_ordering=True))
    return get_newick(tree, "", tree.dist, frame.index)


def format_per_tissue_gene_info(info, tissues):
    for gene in info:
        gene["hasEqtl"] = bool(gene["hasEqtl"])
        gene["hasTwas"] = bool(gene.get("hasTwas", False))
        gene["statusInTissue"] = [
            {"tissueSiteDetailId": tissue,
             "expressed": bool(gene.pop("expr_" + tissue, False)),
             "testedEqtl": bool(gene.pop("testedEqtl_" + tissue, False)),
             "eqtl": bool(gene.pop("eqtl_" + tissue, False))}
            for tissue in tissues
        ]


def parameter(name, default=None):
    values = request.args.getlist(name)
    if len(values) > 1:
        raise APIError(f"Use a single {name} parameter (comma-separated for lists).")
    value = values[0].strip() if values else default
    if value == "":
        raise APIError(f"{name} must not be empty.")
    return value


def identifiers(name, maximum, required=True):
    value = parameter(name)
    if value is None:
        if required:
            raise APIError(f"Missing required parameter: {name}")
        return None
    items = [item.strip() for item in value.split(",")]
    if len(items) > maximum:
        raise APIError(f"{name} accepts at most {maximum} identifiers.")
    if any(not item for item in items):
        raise APIError(f"{name} contains an empty identifier.")
    return list(dict.fromkeys(items))


def boolean(name, default=False):
    value = parameter(name)
    if value is None:
        return default
    if value.lower() not in ("true", "false"):
        raise APIError(f"{name} must be true or false.")
    return value.lower() == "true"


def integer(name, default, minimum=0, maximum=1000):
    try:
        value = int(parameter(name, str(default)))
    except ValueError:
        raise APIError(f"{name} must be an integer.") from None
    if not minimum <= value <= maximum:
        raise APIError(f"{name} must be between {minimum} and {maximum}.")
    return value


def create_app(config=None, store=None):
    app = Flask(__name__)
    app.json = ReleaseJSONProvider(app)
    app.config.from_mapping(
        DATA_DIR=os.environ.get("RATGTEX_DATA_DIR", str(Path(__file__).resolve().parents[1] / "data/v4")),
        EXPRESSION_CACHE_BYTES=int(os.environ.get("RATGTEX_EXPRESSION_CACHE_MB", "128")) * 1024**2,
        RESPONSE_CACHE_BYTES=int(os.environ.get("RATGTEX_RESPONSE_CACHE_MB", "32")) * 1024**2,
        CACHE_SECONDS=300,
        REFERENCE_DATASET=os.environ.get("RATGTEX_REFERENCE_DATASET"),
    )
    if config:
        app.config.update(config)
    store = store if store is not None else DataStore(app.config["DATA_DIR"], app.config["EXPRESSION_CACHE_BYTES"])
    app.extensions["ratgtex_store"] = store
    CORS(app, resources={r"/api/*": {"origins": "*"}}, expose_headers=["ETag"])

    # Cache serialized successful responses within a byte budget, before gzip.
    # Every worker owns its caches; release replacement requires a worker restart.
    responses = OrderedDict()
    response_bytes = 0
    response_lock = RLock()

    def cache_key():
        return request.path, tuple(sorted(request.args.items(multi=True)))

    @app.before_request
    def cached_response():
        nonlocal response_bytes
        if len(request.query_string) > 16384:
            raise APIError("Query string exceeds 16384 bytes.")
        if request.method not in ("GET", "HEAD"):
            return None
        with response_lock:
            cached = responses.pop(cache_key(), None)
            if cached:
                body, expires = cached
                if expires > monotonic():
                    responses[cache_key()] = cached
                    return app.response_class(body, mimetype="application/json")
                response_bytes -= len(body)
        return None

    @app.after_request
    def cache_and_compress(response):
        nonlocal response_bytes
        if response.status_code != 200 or not response.is_json or request.method not in ("GET", "HEAD"):
            response.cache_control.no_store = True
            return response
        body = response.get_data()
        budget = app.config["RESPONSE_CACHE_BYTES"]
        if budget > 0 and len(body) <= budget:
            with response_lock:
                previous = responses.pop(cache_key(), None)
                if previous:
                    response_bytes -= len(previous[0])
                while responses and (response_bytes + len(body) > budget or len(responses) >= 256):
                    _, (removed, _) = responses.popitem(last=False)
                    response_bytes -= len(removed)
                responses[cache_key()] = (body, monotonic() + app.config["CACHE_SECONDS"])
                response_bytes += len(body)
        response.set_etag(hashlib.sha256(body).hexdigest(), weak=True)
        response.cache_control.public = True
        response.cache_control.max_age = app.config["CACHE_SECONDS"]
        response.vary.add("Accept-Encoding")
        response.make_conditional(request)
        if response.status_code == 200 and len(body) >= 1024 and request.accept_encodings["gzip"] > 0:
            compressed = gzip.compress(body, compresslevel=5, mtime=0)
            if len(compressed) < len(body):
                response.set_data(compressed)
                response.headers["Content-Encoding"] = "gzip"
        return response

    @app.errorhandler(APIError)
    def api_error(error):
        return jsonify(error=error.message, status=error.status), error.status

    @app.errorhandler(FileNotFoundError)
    def data_unavailable(error):
        app.logger.error("Release file unavailable: %s", error)
        return jsonify(error="Required release data are unavailable on this server.", status=503), 503

    @app.errorhandler(AssociationDataError)
    def association_unavailable(error):
        app.logger.error("Association inputs unavailable: %s", error)
        return jsonify(error=str(error), status=503), 503

    @app.errorhandler(HTTPException)
    def http_error(error):
        return jsonify(error=error.description, status=error.code), error.code

    def selected_tissues(required=False):
        # 'tissueId' is the preferred spelling; keep the GTEx name as an alias.
        if "tissueId" in request.args and "tissueSiteDetailId" in request.args:
            raise APIError("Use either tissueId or tissueSiteDetailId, not both.")
        key = "tissueId" if "tissueId" in request.args else "tissueSiteDetailId"
        selected = identifiers(key, 100, required=required)
        if selected is None:
            return store.tissues
        unknown = [tissue for tissue in selected if tissue not in store.tissues]
        if unknown:
            raise APIError("Unknown tissue ID(s): " + ", ".join(unknown))
        return selected

    def selected_genes(maximum=500):
        requested = identifiers("geneId", maximum)
        resolved, missing = [], []
        for name in requested:
            gene = store.resolve_gene(name)
            if gene is None:
                missing.append(name)
            elif gene not in resolved:
                resolved.append(gene)
        return resolved, missing

    def one_gene():
        genes, _ = selected_genes(1)
        if not genes:
            raise APIError("Gene not found.", 404)
        return genes[0]

    def one_tissue():
        tissues = selected_tissues(required=True)
        if len(tissues) != 1:
            raise APIError("This endpoint requires exactly one tissue.")
        return tissues[0]

    def expression_scale(default="tpm"):
        scale = parameter("expressionScale", default)
        if scale not in SCALES:
            raise APIError("expressionScale must be tpm or normalized.")
        return scale

    def records(frame):
        return frame.to_dict(orient="records")

    def page(frame, key):
        if "limit" not in request.args and "offset" not in request.args:
            return jsonify({key: records(frame)})
        limit = integer("limit", 100, 1, 1000)
        offset = integer("offset", 0, 0, 10**9)
        return jsonify({key: records(frame.iloc[offset:offset + limit]),
                        "pagination": {"total": len(frame), "limit": limit, "offset": offset}})

    @app.get("/api/v4/metadata")
    def metadata():
        return jsonify(
            dataRelease="v4", genomeAssembly="GRCr8", species="Rattus norvegicus",
            tissueIds=store.tissues, genotypeDatasetIds=list(dict.fromkeys(store.datasets.values())),
            coordinateSystems={"variantId": "chromosome:position; position is 1-based",
                               "genesAndExons": "1-based inclusive"},
            modalities=list(MODALITIES),
            nominalAssociation={"method": "Additive OLS with intercept and modality-specific covariates",
                                "phenotypeScale": "normalized", "missingGenotypes": "mean imputed for inference only",
                                "pValueSource": "computed", "covariateScope": "modality",
                                "multipleTesting": "Precomputed modality-specific gene/group thresholds; no request-time permutations",
                                "cisWindowBp": 1000000, "minimumMAF": .01, "minimumMAC": 2},
            expressionScales={"tpm": "Unnormalized expression in transcripts per million",
                              "normalized": "Quantile-normalized, inverse-normal transformed expression phenotype; no additional request-time covariate adjustment"},
            terminology={"geneId": "Rat gene symbol/annotation identifier, not an Ensembl identifier",
                         "tissueId": "Tissue identifier; tissueSiteDetailId is an accepted query alias",
                         "genotypePhenotypes": "Normalized RNA phenotype vectors and nominal associations for a variant",
                         "maf": "Minor allele frequency among fully called diploid samples",
                         "altAlleleFrequency": "Alternate allele frequency, which may exceed 0.5",
                         "nes": "Nominal regression slope on the normalized phenotype scale",
                         "pValueBeta": "Beta-approximated permutation p-value from the compiled pipeline output"},
            limits={"geneIds": 500, "expressionGeneIds": 10, "variantIds": 100, "clusterGenes": 200, "pageSize": 1000},
            endpoints=sorted(rule.rule for rule in app.url_map.iter_rules() if rule.rule.startswith("/api/")),
            documentation="/about/api/", openapi="/openapi.yaml",
        )

    @app.get("/api/v4/tissueInfo")
    def tissue_info():
        return jsonify(tissueInfo=store.tissue_info)

    @app.get("/api/v4/gene")
    def gene():
        genes, missing = selected_genes()
        info = records(store.genes.loc[genes].reset_index())
        format_per_tissue_gene_info(info, store.tissues)
        return jsonify(gene=info, notFound=missing)

    @app.get("/api/v4/genes")
    def genes_directory():
        query = (parameter("query") or "").casefold()
        limit = integer("limit", 100, 1, 1000)
        offset = integer("offset", 0, 0, 10**9)
        frame = store.genes
        if query:
            mask = frame.index.str.casefold().str.contains(query, regex=False)
            mask |= frame["description"].str.casefold().str.contains(query, regex=False)
            frame = frame.loc[mask]
        info = records(frame.iloc[offset:offset + limit].reset_index())
        format_per_tissue_gene_info(info, store.tissues)
        return jsonify(gene=info, pagination={"total": len(frame), "limit": limit, "offset": offset})

    @app.get("/api/v4/exon")
    def exon():
        gene = one_gene()
        return page(store.by_gene(store.exons, gene), "exon")

    @app.get("/api/v4/eqtl")
    def eqtl():
        gene = one_gene()
        frame = store.by_gene(store.eqtls, gene)
        if "tissueId" in request.args or "tissueSiteDetailId" in request.args:
            frame = frame.loc[frame.tissueSiteDetailId.isin(selected_tissues())]
        return page(frame, "eqtl")

    @app.get("/api/v4/xqtl")
    def xqtl():
        gene = one_gene()
        frame = store.by_gene(store.xqtls, gene)
        if "tissueId" in request.args or "tissueSiteDetailId" in request.args:
            frame = frame.loc[frame.tissueSiteDetailId.isin(selected_tissues())]
        modality = parameter("modality")
        if modality is not None:
            if modality not in set(store.xqtls.modality):
                raise APIError("Unknown modality: " + modality)
            frame = frame.loc[frame.modality == modality]
        return page(frame, "xqtl")

    @app.get("/api/v4/significantCisEqtls")
    def single_tissue_eqtl():
        frame = store.significant_associations(one_gene())
        if "tissueId" in request.args or "tissueSiteDetailId" in request.args:
            tissues = selected_tissues()
            if not frame.empty:
                frame = frame.loc[frame.tissueSiteDetailId.isin(tissues)]
        return page(frame, "significantCisEqtls")

    @app.get("/api/v4/geneExpression")
    def gene_expression():
        genes, missing = selected_genes(10)
        tissues = selected_tissues()
        scale = expression_scale()
        include_samples = boolean("includeSampleIds")
        infos = []
        for tissue in tissues:
            for gene in genes:
                expr = store.expression(tissue, gene, scale)
                if expr is None:
                    continue
                info = {"data": expr.to_numpy(), "datasetId": "ratgtex_v4", "genotypeDatasetId": store.datasets[tissue],
                        "geneId": gene, "geneSymbol": gene, "tissueSiteDetailId": tissue,
                        "unit": "TPM" if scale == "tpm" else "normalized expression",
                        "expressionScale": scale, "sampleCount": len(expr)}
                if include_samples:
                    info["sampleIds"] = expr.index.tolist()
                infos.append(info)
        return jsonify(geneExpression=infos, notFound=missing)

    @app.get("/api/v4/expressionSummary")
    def expression_summary():
        genes, missing = selected_genes(10)
        tissues = selected_tissues()
        scale = expression_scale()
        infos = []
        for tissue in tissues:
            for gene in genes:
                expr = store.expression(tissue, gene, scale)
                if expr is None:
                    continue
                values = expr.to_numpy(dtype=float)
                valid = values[np.isfinite(values)]
                q1, median, q3 = np.quantile(valid, [.25, .5, .75]) if len(valid) else [None] * 3
                infos.append({"geneId": gene, "tissueSiteDetailId": tissue, "expressionScale": scale,
                              "unit": "TPM" if scale == "tpm" else "normalized expression",
                              "sampleCount": len(values), "nonMissingCount": len(valid),
                              "mean": float(valid.mean()) if len(valid) else None,
                              "median": median, "q1": q1, "q3": q3,
                              "min": float(valid.min()) if len(valid) else None,
                              "max": float(valid.max()) if len(valid) else None,
                              "zeroCount": int(np.count_nonzero(valid == 0))})
        return jsonify(expressionSummary=infos, notFound=missing)

    @app.get("/api/v4/medianGeneExpression")
    def median_gene_expression():
        genes, missing = selected_genes()
        tissues = selected_tissues()
        cluster = boolean("cluster", True)
        if cluster and len(genes) > 200:
            raise APIError("Clustering accepts at most 200 genes; use cluster=false for larger requests.")
        available = [gene for gene in genes if gene in store.median_expression.index]
        frame = store.median_expression.loc[available, tissues]
        enough = not frame.empty
        clusters = {"gene": row_tree(frame) if cluster and enough else "Not enough data",
                    "tissue": row_tree(frame.T) if cluster and enough else "Not enough data"}
        result = frame.reset_index().melt(id_vars="geneId", var_name="tissueSiteDetailId", value_name="median")
        return jsonify(medianGeneExpression=records(result), clusters=clusters, notFound=missing)

    @app.get("/api/v4/topExpressedGene")
    def top_expressed_gene():
        tissue = one_tissue()
        frame = store.top_expression.loc[store.top_expression.tissueSiteDetailId == tissue]
        if boolean("filterMtGene"):
            frame = frame.loc[~frame.mtGene]
        limit = integer("limit", 50, 1, 50)
        return jsonify(topExpressedGene=records(frame.iloc[:limit].drop(columns="mtGene")))

    @app.get("/api/v4/variant")
    def variant():
        variants = identifiers("variantId", 100)
        datasets = list(dict.fromkeys(store.datasets.values()))
        default = app.config["REFERENCE_DATASET"] or ("ratgtex_v4_round11_2" if "ratgtex_v4_round11_2" in datasets else datasets[0])
        dataset = parameter("genotypeDatasetId", default)
        if dataset not in datasets:
            raise APIError("Unknown genotypeDatasetId: " + dataset)
        infos, missing = [], []
        with pysam.VariantFile(store.root / "geno" / f"{dataset}.rn8.vcf.gz", drop_samples=True) as vcf:
            for name in variants:
                rec = variant_record(name, vcf)
                if rec is None or not rec.alts or len(rec.alts) != 1:
                    missing.append(name)
                    continue
                infos.append({"alt": rec.alts[0], "chromosome": rec.contig, "pos": rec.pos,
                              "ref": rec.ref, "variantId": name, "genotypeDatasetId": dataset})
        return jsonify(variant=infos, notFound=missing)

    def variant_genotypes(dataset, variant):
        with pysam.VariantFile(store.root / "geno" / f"{dataset}.rn8.vcf.gz") as vcf:
            rec = variant_record(variant, vcf)
            if rec is None:
                raise APIError("Variant not found: " + variant, 404)
            if not rec.alts or len(rec.alts) != 1:
                raise APIError("Genotype phenotypes require a biallelic variant.", 422)
            return pd.Series({sample: genotype(rec.samples[sample]["GT"]) for sample in rec.samples}, dtype=object), rec.ref, rec.alts[0]

    def cis_statistics(stats, variant, locus):
        chrom, pos = variant.split(":")
        is_cis = chrom == locus[0] and locus[1] - 1000000 <= int(pos) <= locus[2] + 1000000
        if not is_cis:
            stats = {**stats, "pValue": None, "slope": None, "slopeSE": None, "associationStatus": "outside_cis_window"}
        return {**stats, "isCis": is_cis}

    @app.get("/api/v4/genotypePhenotypes")
    def genotype_phenotypes():
        gene, tissue = one_gene(), one_tissue()
        variant = identifiers("variantId", 1)[0]
        modality = parameter("modality", "expression")
        if modality not in MODALITIES:
            raise APIError("Unknown modality: " + modality)
        include_samples = boolean("includeSampleIds")
        selected = identifiers("phenotypeId", 1, required=False)
        values, loci, threshold = store.phenotypes(tissue, modality, gene)
        if selected:
            values = values.loc[values.index.isin(selected)]
        if values.empty:
            return jsonify(genotypePhenotypes=[], geneId=gene, tissueSiteDetailId=tissue,
                           variantId=variant, modality=modality)
        dataset = store.datasets[tissue]
        genotypes, ref, alt = variant_genotypes(dataset, variant)
        stats = store.nominal_associations(tissue, modality, values, genotypes)
        geno = [genotypes.get(sample) for sample in values.columns]
        counts = [geno.count(dosage) for dosage in (0, 1, 2)]
        called = sum(counts)
        frequency = (counts[1] + 2 * counts[2]) / (2 * called) if called else None
        results = [{"phenotypeId": phenotype, "data": row.to_list(),
                    **cis_statistics(stat, variant, loci[phenotype])}
                   for (phenotype, row), stat in zip(values.iterrows(), stats)]
        info = {"genotypePhenotypes": results, "geneId": gene, "variantId": variant,
                "modality": modality, "tissueSiteDetailId": tissue, "genotypeDatasetId": dataset,
                "genotypes": geno, "ref": ref, "alt": alt,
                "homoRefCount": counts[0], "hetCount": counts[1], "homoAltCount": counts[2],
                "calledSampleCount": called, "missingGenotypeCount": len(geno) - called,
                "altAlleleFrequency": frequency,
                "maf": min(frequency, 1 - frequency) if frequency is not None else None,
                "pValueThreshold": threshold, "thresholdScope": "modality_gene_group",
                "pValueSource": "computed", "covariateScope": "modality",
                "associationPhenotypeScale": "normalized", "unit": "normalized phenotype"}
        if include_samples:
            info["sampleIds"] = values.columns.tolist()
        return jsonify(info)

    return app


api = create_app()

if __name__ == "__main__":
    api.run()
