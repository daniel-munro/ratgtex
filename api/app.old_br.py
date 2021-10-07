from flask import Flask, g, jsonify, json, request
from flask_cors import CORS
import numpy as np
import os
import pandas as pd
import pysam
from scipy.cluster.hierarchy import linkage, to_tree
import zipfile


def genotype(gt: tuple) -> int:
    return None if gt == (None, None) else gt[0] + gt[1]


def variant_record(variant_id):
    chrom, pos = variant_id.split(":")
    chrom = chrom.replace("chr", "")
    pos = int(pos)
    recs = list(vcf.fetch(chrom, pos - 1, pos, reopen=True))
    assert len(recs) == 1, f"Genotype retrieval error: {variant_id}"
    return recs[0]


def geno_matrix(ids):
    """Assumes SNPs are in close proximity on a chromosome, e.g. in a cis-window."""
    chrom = ids[0].split(":")[0].replace("chr", "")
    pos = [int(x.split(":")[1]) for x in ids]
    genos = {}
    for rec in vcf.fetch(chrom, min(pos) - 1, max(pos) + 1):
        if rec.id in ids:
            genos[rec.id] = [genotype(rec.samples[s]["GT"]) for s in vcf.header.samples]
    mat = np.array([genos[id] for id in ids])
    return mat


def get_newick(node, newick, parentdist, leaf_names):
    """from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format/31878514#31878514"""
    if node.is_leaf():
        return "%s:%g%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%g%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


def row_tree(d):
    clust = linkage(d, method="average", optimal_ordering=True)
    tree = to_tree(clust)
    return get_newick(tree, "", tree.dist, d.index)


def validate_genes(ids, genes):
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


def load_tpm(path):
    tpm = {}
    expr = pd.read_csv(path, sep="\t")
    samples = pd.read_csv("../data/ref/metadata.csv")
    samples = samples.loc[samples["QC_pass"] == "pass", :]
    expr = expr.loc[:, expr.columns.isin(samples["library"])]
    tis_conv = {"Acbc": "NAcc", "IL": "IL", "LHB": "LHb", "PL": "PL", "VoLo": "OFC"}
    tis = pd.Series([tis_conv[x.split("_")[1]] for x in expr.columns])
    for tissue in tis.unique():
        tpm[tissue] = expr.loc[:, list(tis == tissue)]
    return tpm


def cis_pval(tissue, gene, variant):
    # g.db = sqlite3.connect("../data/cis_pvals.db")
    # res = g.db.execute(f'SELECT * FROM {tissue} WHERE gene_id = "{gene}" AND variant_id = "{variant}"')
    with zipfile.ZipFile(f"../data/cis_pvals/{tissue}.zip", "r") as archive:
        fname = f"{tissue}/{gene}.txt"
        if fname in archive.namelist():
            df = pd.read_csv(archive.open(fname), sep="\t", index_col="variant_id")
            if variant in df.index:
                return df.loc[variant, "pval_nominal"]
        return None


def single_tis(gene):
    with zipfile.ZipFile(f"../data/singleTissueEqtl.zip", "r") as archive:
        fname = f"singleTissueEqtl/{gene}.txt"
        if fname in archive.namelist():
            d = pd.read_csv(archive.open(fname), sep="\t", dtype={"chromosome": str})
            d["geneId"] = gene
            return d
        return None


tissueInfo = pd.read_csv("../data/tissueInfo.txt", sep="\t")
tissueInfo = tissueInfo.to_dict(orient="records")

topExpr = pd.read_csv("../data/topExpressedGene.txt", sep="\t")

genes = pd.read_csv("../data/gene.txt", sep="\t", index_col="geneId").fillna("")

tissues = ["IL", "LHb", "NAcc", "OFC", "PL"]
# med_expr = []
# for tissue in tissues:
#     tissue_sam = samples.library[samples.brain_region == tissue]
#     med_expr.append(tpm[tissue_sam].median(axis=1))
# med_expr = pd.concat(med_expr, axis=1)
# med_expr.columns = tissues
# med_expr["geneId"] = med_expr.index
med_expr = pd.read_csv("../data/medianGeneExpression.txt.gz", sep="\t", index_col="geneId")
tpm = load_tpm("../data/expr/ensembl-gene_raw-tpm.txt")
iqn = {}
for tissue in tissues:
    # iqn_file = f"../data/expr/ensembl-gene_inv-quant_PCresiduals_{tissue}.txt.gz"
    iqn_file = f"../data/expr/ensembl-gene_inv-quant_{tissue}.bed.gz"
    iqn[tissue] = pd.read_csv(iqn_file, sep="\t", index_col="gene_id")
    iqn[tissue].drop(columns=["#chr", "start", "end"], inplace=True)

vcf = pysam.VariantFile("../data/P50.rnaseq.88.unpruned.vcf.gz")

# Load all significant pairs for singleTissueEqtl
# all_sig = pd.read_csv("../data/singleTissueEqtl.txt.gz", sep="\t")

exons = pd.read_csv("../data/exon.txt", sep="\t", low_memory=False)

top_assoc = pd.read_csv("../data/eqtl/top_assoc.txt", sep="\t", index_col=["tissue", "gene_id"])  # Just for pval_nominal_threshold
eqtls = pd.read_csv("../data/eqtl/eqtls_indep.txt", sep="\t")
eqtls = eqtls[["tissue", "gene_id", "gene_name", "tss", "variant_id", "ref", "alt", "pval_beta", "log2_aFC"]]
eqtls = eqtls.rename(columns={"tissue": "tissueSiteDetailId", "gene_id": "geneId", "gene_name": "geneSymbol", "variant_id": "variantId"})

# cis_pvals = {}
# for tissue in tissues:
#     df = pd.read_csv(f"~/br/data/tensorqtl/all_cis_pvals/{tissue[0]}QCT.all_cis_pvals.txt.gz", sep="\t", index_col=["phenotype_id", "variant_id"])
#     df = df.iloc[:5, :]
#     x = {gene: {var: df.loc[(gene, var), 'pval_nominal'] for var in df.xs(gene).index.levels[0]} for gene in df.index.levels[0]}
#     cis_pvals[tissue] = x
# print(cis_pvals)

api = Flask(__name__)
CORS(api)
# api.config["APPLICATION_ROOT"] = "/api/v1" # doesn't work??


@api.route("/api/v1/dyneqtl", methods=["GET"])
def dyneqtl():
    variant = request.args.get("variantId")
    gene = request.args.get("geneId")
    tissue = request.args.get("tissueSiteDetailId")
    expr = iqn[tissue].loc[gene, :]
    rec = variant_record(variant)
    assert len(rec.alts) == 1, f"Multiple alt alleles: {variant}"
    gt = rec.samples
    indivs = [x.split("_")[0] for x in expr.index]
    geno = [genotype(gt[ind]["GT"]) for ind in indivs]
    # ignoring error, nes, tStatistic, timing
    counts = [int(np.sum(np.array(geno) == x)) for x in [0, 1, 2]]
    # pval = None
    # if genes.loc[gene, "hasEqtl"]:
    #     # print('a', variant, gene, tissue)
    #     d = pd.read_csv(f"../data/singleTissueEqtl/{gene}.txt.gz", sep="\t")
    #     d = d.loc[(d["tissueSiteDetailId"] == tissue) & (d["variantId"] == variant), :]
    #     if d.shape[0] > 0:
    #         pval = d["pValue"].iloc[0]
    #     # print('b', variant, gene, tissue)
    # if cis_pvals.index.contains((gene, variant)):
    #     pval = cis_pvals[gene][varia]
    # if gene in cis_pvals[tissue] and variant in cis_pvals[tissue][gene]:
    #     pval = cis_pvals[tissue][gene][variant]
    pval = cis_pval(tissue, gene, variant)
    # thresh = list(top_assoc["pval_nominal_threshold"].loc[(top_assoc["tissue"] == tissue) & (top_assoc["gene_id"] == gene)])[0]
    # print('c', variant, gene, tissue)
    thresh = top_assoc.loc[(tissue, gene), "pval_nominal_threshold"]
    # print('d', variant, gene, tissue)
    info = {
        "data": list(expr),
        "geneId": gene,
        "geneSymbol": genes.loc[gene, "geneSymbol"],
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


@api.route("/api/v1/exon", methods=["GET"])
def exon():
    gene = request.args.get("geneId")
    # ignoring: gencodeVersion, genomeBuild
    # chromosome, end, exonId, exonNumber, featureType, geneId, geneSymbol, start, strand, transcriptId
    d = exons.loc[exons["geneId"] == gene, :]
    d = d.to_dict(orient="records")
    return jsonify({"exon": d})


@api.route("/api/v1/gene", methods=["GET"])
def gene():
    # if request.args.get("geneId") is not None:
    #     ids = request.args.get("geneId").split(",")
    #     # TEMP:
    #     ids = [x for x in ids if x in genes.index]
    #     d = genes.loc[ids, :].reset_index() # Include geneId in dict
    # else:
    #     ids = request.args.get("geneSymbol").split(",")
    #     genes2 = genes.reset_index().set_index("geneSymbol")
    #     ids = [x for x in ids if x in genes2.index]
    #     d = genes2.loc[ids, :].reset_index()
    ids = request.args.get("geneId").split(",")
    ids = validate_genes(ids, genes)
    d = genes.loc[ids, :].reset_index() # Include geneId in dict
    info = d.to_dict(orient="records")
    return jsonify({"gene": info})


@api.route("/api/v1/geneExpression", methods=["GET"])
def gene_exp():
    gene = request.args.get("geneId")
    symbol = genes.loc[gene, "geneSymbol"]
    infos = []
    for tissue in tissues:
        info = {
            "data": list(tpm[tissue].loc[gene, :]),
            "datasetId": "ratgtex_v1",
            "geneId": gene,
            "geneSymbol": symbol,
            "tissueSiteDetailId": tissue,
            "unit": "TPM"
        }
        infos.append(info)
    return jsonify({"geneExpression": infos})


@api.route("/api/v1/ld", methods=["GET"])
def ld():
    gene = request.args.get("geneId")
    # d = pd.read_csv(f"../data/singleTissueEqtl/{gene}.txt.gz", sep="\t")
    d = single_tis(gene)
    if d is None:
        return jsonify({"ld": []})
    d["pos"] = [int(x.split(":")[1]) for x in d["variantId"]]
    d = d.sort_values(by="pos")
    ids = d["variantId"].unique()
    geno = geno_matrix(ids)
    ldmat = np.corrcoef(geno) ** 2
    ldmat = ldmat.round(3)
    lds = []
    for i in range(len(ids) - 1):
        for j in range(i + 1, len(ids)):
            # lds.append([",".join([ids[i], ids[j]]), ldmat[i, j]])
            lds.append([ids[i], ids[j], ldmat[i, j]])
    return jsonify({"ld": lds})


@api.route("/api/v1/medianGeneExpression", methods=["GET"])
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
    d = d.reset_index().melt(id_vars="geneId", var_name="tissueSiteDetailId", value_name="median")
    d = d.merge(genes["geneSymbol"].reset_index(), how="left", on="geneId")
    d = d.to_dict(orient="records")
    info = {"clusters": {"gene": gene_tree, "tissue": tissue_tree}, "medianGeneExpression": d}
    return jsonify(info)


@api.route("/api/v1/singleTissueEqtl", methods=["GET"])
def single_tissue_eqtl():
    gene = request.args.get("geneId")
    # ignoring datasetId, snpId
    # chromosome, geneId, geneSymbol, geneSymbolUpper, nes, pValue, pos, tissueSiteDetailId, variantId
    # if not os.path.exists(f"../data/singleTissueEqtl/{gene}.txt.gz"):
    #     return jsonify({"singleTissueEqtl": []})
    # d = pd.read_csv(f"../data/singleTissueEqtl/{gene}.txt.gz", sep="\t", dtype={"chromosome": str})
    d = single_tis(gene)
    if d is None:
        return jsonify({"singleTissueEqtl": []})
    d["geneSymbol"] = genes.loc[gene, "geneSymbol"]
    d["geneSymbolUpper"] = d["geneSymbol"]
    info = d.to_dict(orient="records")
    return jsonify({"singleTissueEqtl": info})


@api.route("/api/v1/eqtl", methods=["GET"])
def eqtl():
    gene = request.args.get("geneId")
    d = eqtls.loc[eqtls["geneId"] == gene, :]
    info = d.to_dict(orient="records")
    return jsonify({"eqtl": info})


@api.route("/api/v1/tissueInfo", methods=["GET"])
def tissue_info():
    # return open("../data/tissueInfo.json", "r").read()
    return jsonify({"tissueInfo": tissueInfo})


@api.route("/api/v1/topExpressedGene", methods=["GET"])
def top_expressed_gene():
    tissue = request.args.get("tissueSiteDetailId")
    filterMt = request.args.get("filterMtGene", type=json.loads, default=False)
    # return open(f"../data/topExpressedGene/{tissue}.json", "r").read()
    x = topExpr.loc[topExpr["tissueSiteDetailId"] == tissue, :]
    if filterMt:
        x = x.loc[~x["mtGene"], :]
    x = x.iloc[:50, :]
    x = x.drop(columns="mtGene")
    x = x.to_dict(orient="records")
    return jsonify({"topExpressedGene": x})


@api.route("/api/v1/variant", methods=["GET"])
def variant():
    ids = request.args.get("variantId").split(",")
    infos = []
    for variant in ids:
        # Ignoring b37VariantId, datasetId, maf01, shorthand, snpId
        rec = variant_record(variant)
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
