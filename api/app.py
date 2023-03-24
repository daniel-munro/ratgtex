from flask import Flask, g, jsonify, json, request
from flask_cors import CORS
import numpy as np
import os
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
    chrom = chrom.replace("chr", "")
    pos = int(pos)
    recs = list(vcf.fetch(chrom, pos - 1, pos, reopen=True))
    assert len(recs) == 1, f"Genotype retrieval error: {variant_id}"
    return recs[0]


def geno_matrix(ids, vcf):
    """Get genotype matrix for a list of SNPs
    Assumes SNPs are in close proximity on a chromosome, e.g. in a cis-window.
    """
    chrom = ids[0].split(":")[0].replace("chr", "")
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
            x = list(genes.loc[genes["geneSymbol"] == id, :].index)
            if len(x) > 0:
                valid.append(x[0])
            else:
                id2 = id[0].upper() + id[1:].lower()
                x = list(genes.loc[genes["geneSymbol"] == id2, :].index)
                if len(x) > 0:
                    valid.append(x[0])
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


def cis_pval(tissue, genome, gene, variant):
    """Return nominal p-value for a given cis-window variant"""
    with zipfile.ZipFile(f"../data/cis_pvals/{tissue}.{genome}.zip", "r") as archive:
        fname = f"{tissue}.{genome}/{gene}.txt"
        if fname in archive.namelist():
            df = pd.read_csv(archive.open(fname), sep="\t", index_col="variant_id")
            if variant in df.index:
                return df.loc[variant, "pval_nominal"]
        return None


def single_tissue(genome, gene):
    """Return table of significant cis-eSNPs for a gene"""
    with zipfile.ZipFile(f"../data/{genome}.singleTissueEqtl.zip", "r") as archive:
        fname = f"{genome}.singleTissueEqtl/{gene}.txt"
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
            "gene_name",
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
            "gene_name": "geneSymbol",
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
            "gene_name",
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
            "gene_name": "geneSymbol",
            "variant_id": "variantId",
        }
    )
    return sqtls


tissueInfo = {}
topExpr = {}
genes = {}
tissues = {}
dataset = {}
med_expr = {}
tpm = {}
iqn = {}
vcf = {}
ref_vcf = {}
exons = {}
top_assoc = {}
eqtls = {}
sqtls = {}

for genome in ["rn6", "rn7"]:
    df = pd.read_csv(f"../data/{genome}.tissueInfo.txt", sep="\t")
    tissueInfo[genome] = df.to_dict(orient="records")

    topExpr[genome] = pd.read_csv(f"../data/{genome}.topExpressedGene.txt", sep="\t")

    df = pd.read_csv(f"../data/{genome}.gene.txt", sep="\t", index_col="geneId", dtype={"chromosome": str})
    genes[genome] = df.fillna("")

    tissues[genome] = [tissue["tissueSiteDetailId"] for tissue in tissueInfo[genome]]
    dataset[genome] = {tissue["tissueSiteDetailId"]: tissue["dataset"] for tissue in tissueInfo[genome]}

    med_expr[genome] = pd.read_csv(
        f"../data/{genome}.medianGeneExpression.txt.gz", sep="\t", index_col="geneId"
    )

    tpm[genome] = {}
    for tissue in tissues[genome]:
        fname = f"../data/expr/{tissue}.{genome}.expr.tpm.bed.gz"
        tpm[genome][tissue] = load_expr(fname)

    iqn[genome] = {}
    for tissue in tissues[genome]:
        fname = f"../data/expr/{tissue}.{genome}.expr.iqn.bed.gz"
        iqn[genome][tissue] = load_expr(fname)

    vcf[genome] = {}
    for dset in set(dataset[genome].values()):
        vcf[genome][dset] = pysam.VariantFile(f"../data/geno/{dset}.{genome}.vcf.gz")
    ref_vcf[genome] = vcf[genome]["BLA_NAcc2_PL2"]

    exons[genome] = pd.read_csv(f"../data/{genome}.exon.txt", sep="\t", dtype={"chromosome": str})

    top_assoc[genome] = pd.read_csv(
        f"../data/eqtl/{genome}.top_assoc.txt", sep="\t", index_col=["tissue", "gene_id"]
    )  # Just for pval_nominal_threshold
    eqtls[genome] = load_eqtls(f"../data/eqtl/{genome}.eqtls_indep.txt")
    sqtls[genome] = load_sqtls(f"../data/splice/{genome}.sqtls_indep.txt")

api = Flask(__name__)
CORS(api)
# api.config["APPLICATION_ROOT"] = "/api/v1" # doesn't work??


@api.route("/api/v2/dyneqtl", methods=["GET"])
def dyneqtl():
    variant = request.args.get("variantId")
    gene = request.args.get("geneId")
    tissue = request.args.get("tissueSiteDetailId")
    genome = request.args.get("genome")
    expr = iqn[genome][tissue].loc[gene, :]
    rec = variant_record(variant, vcf[genome][dataset[genome][tissue]])
    assert len(rec.alts) == 1, f"Multiple alt alleles: {variant}"
    gt = rec.samples
    # indivs = [x.split("_")[0] for x in expr.index]
    geno = [genotype(gt[ind]["GT"]) for ind in expr.index]
    # ignoring error, nes, tStatistic, timing
    counts = [int(np.sum(np.array(geno) == x)) for x in [0, 1, 2]]
    pval = cis_pval(tissue, genome, gene, variant)
    thresh = top_assoc[genome].loc[(tissue, gene), "pval_nominal_threshold"]
    info = {
        "data": list(expr),
        "geneId": gene,
        "geneSymbol": genes[genome].loc[gene, "geneSymbol"],
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


@api.route("/api/v2/eqtl", methods=["GET"])
def eqtl():
    gene = request.args.get("geneId")
    genome = request.args.get("genome")
    d = eqtls[genome].loc[eqtls[genome]["geneId"] == gene, :]
    info = d.to_dict(orient="records")
    return jsonify({"eqtl": info})


@api.route("/api/v2/exon", methods=["GET"])
def exon():
    gene = request.args.get("geneId")
    genome = request.args.get("genome")
    d = exons[genome].loc[exons[genome]["geneId"] == gene, :]
    d = d.to_dict(orient="records")
    return jsonify({"exon": d})


@api.route("/api/v2/gene", methods=["GET"])
def gene():
    ids = request.args.get("geneId").split(",")
    genome = request.args.get("genome")
    ids = validate_genes(ids, genes[genome])
    d = genes[genome].loc[ids, :].reset_index()  # Include geneId in dict
    info = d.to_dict(orient="records")
    format_per_tissue_gene_info(info, tissues[genome])
    return jsonify({"gene": info})


@api.route("/api/v2/geneExpression", methods=["GET"])
def gene_exp():
    gene = request.args.get("geneId")
    genome = request.args.get("genome")
    symbol = genes[genome].loc[gene, "geneSymbol"]
    infos = []
    for tissue in tissues[genome]:
        info = {
            "data": list(tpm[genome][tissue].loc[gene, :]),
            "datasetId": "ratgtex_v1",
            "geneId": gene,
            "geneSymbol": symbol,
            "tissueSiteDetailId": tissue,
            "unit": "TPM",
        }
        infos.append(info)
    return jsonify({"geneExpression": infos})


@api.route("/api/v2/ld", methods=["GET"])
def ld():
    gene = request.args.get("geneId")
    genome = request.args.get("genome")
    d = single_tissue(genome, gene)
    if d is None:
        return jsonify({"ld": []})
    d["pos"] = [int(x.split(":")[1]) for x in d["variantId"]]
    d = d.sort_values(by="pos")
    ids = d["variantId"].unique()
    geno = geno_matrix(ids, ref_vcf[genome])
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


@api.route("/api/v2/medianGeneExpression", methods=["GET"])
def med_gene_exp():
    ids = request.args.get("geneId").split(",")
    genome = request.args.get("genome")
    ids = [x for x in ids if x in med_expr[genome].index]
    if request.args.get("tissueSiteDetailId") is None:
        d = med_expr[genome].loc[ids, :]
    else:
        tissues = request.args.get("tissueSiteDetailId").split(",")
        d = med_expr[genome].loc[ids, tissues]
    gene_tree = row_tree(d)
    tissue_tree = row_tree(d.T)
    d = d.reset_index().melt(
        id_vars="geneId", var_name="tissueSiteDetailId", value_name="median"
    )
    d = d.merge(genes[genome]["geneSymbol"].reset_index(), how="left", on="geneId")
    d = d.to_dict(orient="records")
    info = {
        "clusters": {"gene": gene_tree, "tissue": tissue_tree},
        "medianGeneExpression": d,
    }
    return jsonify(info)


@api.route("/api/v2/singleTissueEqtl", methods=["GET"])
def single_tissue_eqtl():
    gene = request.args.get("geneId")
    genome = request.args.get("genome")
    d = single_tissue(genome, gene)
    if d is None:
        return jsonify({"singleTissueEqtl": []})
    d["geneSymbol"] = genes[genome].loc[gene, "geneSymbol"]
    d["geneSymbolUpper"] = d["geneSymbol"]
    info = d.to_dict(orient="records")
    return jsonify({"singleTissueEqtl": info})


@api.route("/api/v2/sqtl", methods=["GET"])
def sqtl():
    gene = request.args.get("geneId")
    genome = request.args.get("genome")
    d = sqtls[genome].loc[sqtls[genome]["geneId"] == gene, :]
    info = d.to_dict(orient="records")
    return jsonify({"sqtl": info})


@api.route("/api/v2/tissueInfo", methods=["GET"])
def tissue_info():
    genome = request.args.get("genome")
    return jsonify({"tissueInfo": tissueInfo[genome]})


@api.route("/api/v2/topExpressedGene", methods=["GET"])
def top_expressed_gene():
    tissue = request.args.get("tissueSiteDetailId")
    genome = request.args.get("genome")
    filterMt = request.args.get("filterMtGene", type=json.loads, default=False)
    x = topExpr[genome].loc[topExpr[genome]["tissueSiteDetailId"] == tissue, :]
    if filterMt:
        x = x.loc[~x["mtGene"], :]
    x = x.iloc[:50, :]
    x = x.drop(columns="mtGene")
    x = x.to_dict(orient="records")
    return jsonify({"topExpressedGene": x})


@api.route("/api/v2/variant", methods=["GET"])
def variant():
    ids = request.args.get("variantId").split(",")
    genome = request.args.get("genome")
    infos = []
    for variant in ids:
        # Ignoring b37VariantId, datasetId, maf01, shorthand, snpId
        rec = variant_record(variant, ref_vcf[genome])
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
