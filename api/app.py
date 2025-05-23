from flask import Flask, jsonify, json, request
from flask_cors import CORS
import numpy as np
import pandas as pd
import pysam
from scipy.cluster.hierarchy import linkage, to_tree
import zipfile


def genotype(gt: tuple) -> int:
    """Convert genotype tuple to dosage (0/1/2)"""
    return None if gt == (None, None) else gt[0] + gt[1]


def variant_record(variant_id, vcf):
    """Get record for one variant from VCF"""
    chrom, pos = variant_id.split(":")
    pos = int(pos)
    recs = list(vcf.fetch(chrom, pos - 1, pos, reopen=True))
    assert len(recs) == 1, f"Genotype retrieval error: {variant_id}"
    return recs[0]


def geno_matrix(ids, vcf):
    """Get genotype matrix for a list of SNPs
    Assumes SNPs are in close proximity on a chromosome, e.g. in a cis-window.
    """
    chrom = ids[0].split(":")[0]
    pos = [int(x.split(":")[1]) for x in ids]
    genos = {}
    for rec in vcf.fetch(chrom, min(pos) - 1, max(pos) + 1):
        if rec.id in ids:
            genos[rec.id] = [genotype(rec.samples[s]["GT"]) for s in vcf.header.samples]
    mat = np.array([genos[id] if id in genos else [None] * len(vcf.header.samples) for id in ids])
    return mat


def get_newick(node, newick, parentdist, leaf_names):
    """Save dendrogram in Newick format
    from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format/31878514#31878514
    """
    if node.is_leaf():
        return "%s:%g%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    if len(newick) > 0:
        newick = "):%g%s" % (parentdist - node.dist, newick)
    else:
        newick = ");"
    newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
    newick = get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
    newick = "(%s" % (newick)
    return newick


def row_tree(d):
    """Get Newick representation of matrix for clustering"""
    clust = linkage(d, method="average", optimal_ordering=True)
    tree = to_tree(clust)
    return get_newick(tree, "", tree.dist, d.index)


def validate_genes(ids, genes):
    """Return valid gene IDs for a list of gene IDs/names"""
    valid = []
    for id in ids:
        if id in genes.index:
            valid.append(id)
        else:
            id2 = id[0].upper() + id[1:].lower()
            if id2 in genes.index:
                valid.append(id2)
    return valid


def format_per_tissue_gene_info(info: list, tissues: list):
    """Collect per-tissue expression, eQTL, and sQTL indicators into a list"""
    for gene in info:
        gene["statusInTissue"] = []
        for tissue in tissues:
            item = {
                "tissueSiteDetailId": tissue,
                "expressed": gene["expr_" + tissue],
                "testedEqtl": gene["testedEqtl_" + tissue],
                "eqtl": gene["eqtl_" + tissue],
                "altSplice": gene["altSplice_" + tissue],
                "testedSqtl": gene["testedSqtl_" + tissue],
                "sqtl": gene["sqtl_" + tissue],
            }
            gene["statusInTissue"].append(item)
            del gene["expr_" + tissue]
            del gene["testedEqtl_" + tissue]
            del gene["eqtl_" + tissue]
            del gene["altSplice_" + tissue]
            del gene["testedSqtl_" + tissue]
            del gene["sqtl_" + tissue]


def cis_pval(tissue, gene, variant):
    """Return nominal p-value for a given cis-window variant"""
    with zipfile.ZipFile(f"../data/cis_pvals/{tissue}.v3.zip", "r") as archive:
        fname = f"{tissue}.v3/{gene}.txt"
        if fname in archive.namelist():
            df = pd.read_csv(archive.open(fname), sep="\t", index_col="variant_id")
            if variant in df.index:
                return df.loc[variant, "pval_nominal"]
        return None


def single_tissue(gene):
    """Return table of significant cis-eSNPs for a gene"""
    with zipfile.ZipFile(f"../data/singleTissueEqtl.v3.zip", "r") as archive:
        fname = f"singleTissueEqtl.v3/{gene}.txt"
        if fname in archive.namelist():
            d = pd.read_csv(archive.open(fname), sep="\t", dtype={"chromosome": str})
            d["geneId"] = gene
            return d
        return None


def load_expr(fname):
    df = pd.read_csv(fname, sep="\t", dtype={"#chr": str}, index_col="gene_id")
    df.drop(columns=["#chr", "start", "end"], inplace=True)
    return df


def load_eqtls(fname):
    eqtls = pd.read_csv(fname, sep="\t")
    eqtls = eqtls[
        [
            "tissue",
            "gene_id",
            "variant_id",
            "ref",
            "alt",
            "pval_beta",
            "log2_aFC",
        ]
    ]
    eqtls = eqtls.rename(
        columns={
            "tissue": "tissueSiteDetailId",
            "gene_id": "geneId",
            "variant_id": "variantId",
        }
    )
    return eqtls


def load_sqtls(fname):
    sqtls = pd.read_csv(fname, sep="\t")
    sqtls = sqtls[
        [
            "tissue",
            "phenotype_id",
            "gene_id",
            "variant_id",
            "ref",
            "alt",
            "pval_beta",
        ]
    ]
    sqtls = sqtls.rename(
        columns={
            "tissue": "tissueSiteDetailId",
            "phenotype_id": "phenotypeId",
            "gene_id": "geneId",
            "variant_id": "variantId",
        }
    )
    return sqtls


df = pd.read_csv("../data/tissueInfo.v3.txt", sep="\t")
tissueInfo = df.to_dict(orient="records")

topExpr = pd.read_csv("../data/topExpressedGene.v3.txt", sep="\t")

df = pd.read_csv("../data/gene.v3.txt", sep="\t", index_col="geneId", dtype={"chromosome": str})
genes = df.fillna("")

tissues = [tissue["tissueSiteDetailId"] for tissue in tissueInfo]
dataset = {tissue["tissueSiteDetailId"]: tissue["dataset"] for tissue in tissueInfo}

med_expr = pd.read_csv(
    "../data/medianGeneExpression.v3.txt.gz", sep="\t", index_col="geneId"
)

tpm = {}
for tissue in tissues:
    fname = f"../data/expr/expr.tpm.{tissue}.v3_rn7.bed.gz"
    tpm[tissue] = load_expr(fname)

iqn = {}
for tissue in tissues:
    fname = f"../data/expr/expr.iqn.{tissue}.v3_rn7.bed.gz"
    iqn[tissue] = load_expr(fname)

vcf = {}
for dset in set(dataset.values()):
    vcf[dset] = pysam.VariantFile(f"../data/geno/{dset}.rn7.vcf.gz")
ref_vcf = vcf["ratgtex_v3_round10_5"]

exons = pd.read_csv("../data/exon.v3.txt", sep="\t", dtype={"chromosome": str})

top_assoc = pd.read_csv(
    "../data/eqtl/top_assoc.v3_rn7.txt", sep="\t", index_col=["tissue", "gene_id"]
)
top_assoc = top_assoc[["pval_nominal_threshold"]]
eqtls = load_eqtls("../data/eqtl/eqtls_indep.v3_rn7.txt")
sqtls = load_sqtls("../data/splice/sqtls_indep.v3_rn7.txt")

api = Flask(__name__)
CORS(api)


@api.route("/api/v3/dyneqtl", methods=["GET"])
def dyneqtl():
    variant = request.args.get("variantId")
    gene = request.args.get("geneId")
    tissue = request.args.get("tissueSiteDetailId")
    expr = iqn[tissue].loc[gene, :]
    rec = variant_record(variant, vcf[dataset[tissue]])
    assert len(rec.alts) == 1, f"Multiple alt alleles: {variant}"
    gt = rec.samples
    geno = [genotype(gt[ind]["GT"]) for ind in expr.index]
    # ignoring error, nes, tStatistic, timing
    counts = [int(np.sum(np.array(geno) == x)) for x in [0, 1, 2]]
    pval = cis_pval(tissue, gene, variant)
    thresh = top_assoc.loc[(tissue, gene), "pval_nominal_threshold"]
    info = {
        "data": list(expr),
        "geneId": gene,
        "genotypes": geno,
        "hetCount": counts[1],
        "homoAltCount": counts[2],
        "homoRefCount": counts[0],
        "maf": (counts[1] + 2 * counts[2]) / int(2 * np.sum(counts)),
        "ref": rec.ref,  # I added ref and alt to API since they aren't in our variant IDs
        "alt": rec.alts[0],
        "pValue": pval,
        "pValueThreshold": thresh,
        "tissueSiteDetailId": tissue,
        "variantId": variant,
    }
    return jsonify(info)


@api.route("/api/v3/eqtl", methods=["GET"])
def eqtl():
    gene = request.args.get("geneId")
    d = eqtls.loc[eqtls["geneId"] == gene, :]
    info = d.to_dict(orient="records")
    return jsonify({"eqtl": info})


@api.route("/api/v3/exon", methods=["GET"])
def exon():
    gene = request.args.get("geneId")
    d = exons.loc[exons["geneId"] == gene, :]
    d = d.to_dict(orient="records")
    return jsonify({"exon": d})


@api.route("/api/v3/gene", methods=["GET"])
def gene():
    ids = request.args.get("geneId").split(",")
    ids = validate_genes(ids, genes)
    d = genes.loc[ids, :].reset_index()  # Include geneId in dict
    info = d.to_dict(orient="records")
    format_per_tissue_gene_info(info, tissues)
    return jsonify({"gene": info})


@api.route("/api/v3/geneExpression", methods=["GET"])
def gene_exp():
    gene = request.args.get("geneId")
    infos = []
    for tissue in tissues:
        info = {
            "data": list(tpm[tissue].loc[gene, :]),
            "datasetId": "ratgtex_v1",
            "geneId": gene,
            "tissueSiteDetailId": tissue,
            "unit": "TPM",
        }
        infos.append(info)
    return jsonify({"geneExpression": infos})


@api.route("/api/v3/ld", methods=["GET"])
def ld():
    gene = request.args.get("geneId")
    d = single_tissue(gene)
    if d is None:
        return jsonify({"ld": []})
    ids = d["variantId"].unique()
    geno = geno_matrix(ids, ref_vcf)
    # ldmat = np.corrcoef(geno) ** 2
    geno = pd.DataFrame(geno.T, dtype=float)  # Pandas corr allows missing values
    ldmat = geno.corr().to_numpy() ** 2
    ldmat = ldmat.round(3)
    lds = []
    for i in range(len(ids) - 1):
        for j in range(i + 1, len(ids)):
            ld = ldmat[i, j] if not np.isnan(ldmat[i, j]) else None
            lds.append([ids[i], ids[j], ld])
    return jsonify({"ld": lds})


@api.route("/api/v3/medianGeneExpression", methods=["GET"])
def med_gene_exp():
    ids = request.args.get("geneId").split(",")
    ids = [x for x in ids if x in med_expr.index]
    if request.args.get("tissueSiteDetailId") is None:
        d = med_expr.loc[ids, :]
    else:
        tissues = request.args.get("tissueSiteDetailId").split(",")
        d = med_expr.loc[ids, tissues]
    gene_tree = row_tree(d)
    tissue_tree = row_tree(d.T)
    d = d.reset_index().melt(
        id_vars="geneId", var_name="tissueSiteDetailId", value_name="median"
    )
    d = d.to_dict(orient="records")
    info = {
        "clusters": {"gene": gene_tree, "tissue": tissue_tree},
        "medianGeneExpression": d,
    }
    return jsonify(info)


@api.route("/api/v3/singleTissueEqtl", methods=["GET"])
def single_tissue_eqtl():
    gene = request.args.get("geneId")
    d = single_tissue(gene)
    if d is None:
        return jsonify({"singleTissueEqtl": []})
    # Format p-values to 3 decimal places
    if 'pValue' in d.columns:
        d['pValue'] = d['pValue'].apply(lambda x: float(f"{x:.3e}") if pd.notnull(x) else x)
    info = d.to_dict(orient="records")
    return jsonify({"singleTissueEqtl": info})


@api.route("/api/v3/sqtl", methods=["GET"])
def sqtl():
    gene = request.args.get("geneId")
    d = sqtls.loc[sqtls["geneId"] == gene, :]
    info = d.to_dict(orient="records")
    return jsonify({"sqtl": info})


@api.route("/api/v3/tissueInfo", methods=["GET"])
def tissue_info():
    return jsonify({"tissueInfo": tissueInfo})


@api.route("/api/v3/topExpressedGene", methods=["GET"])
def top_expressed_gene():
    tissue = request.args.get("tissueSiteDetailId")
    filterMt = request.args.get("filterMtGene", type=json.loads, default=False)
    x = topExpr.loc[topExpr["tissueSiteDetailId"] == tissue, :]
    if filterMt:
        x = x.loc[~x["mtGene"], :]
    x = x.iloc[:50, :]
    x = x.drop(columns="mtGene")
    x = x.to_dict(orient="records")
    return jsonify({"topExpressedGene": x})


@api.route("/api/v3/variant", methods=["GET"])
def variant():
    ids = request.args.get("variantId").split(",")
    infos = []
    for variant in ids:
        # Ignoring b37VariantId, datasetId, maf01, shorthand, snpId
        rec = variant_record(variant, ref_vcf)
        assert len(rec.alts) == 1, f"Multiple alt alleles: {variant}"
        info = {
            "alt": rec.alts[0],
            "chromosome": rec.contig,
            "pos": rec.pos,
            "ref": rec.ref,
            "variantId": variant,
        }
        infos.append(info)
    return jsonify({"variant": infos})


if __name__ == "__main__":
    api.run()
