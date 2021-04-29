from flask import Flask, json, request
import pandas as pd
import numpy as np
# from scipy.cluster.hierarchy import linkage
import pysam


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


genes = pd.read_csv("../data/gene.txt", sep="\t", index_col="gencodeId").fillna("")

tissues = ["Acbc", "IL", "LHB", "PL", "VoLo"]
# med_expr = []
# for tissue in tissues:
#     tissue_sam = samples.library[samples.brain_region == tissue]
#     med_expr.append(tpm[tissue_sam].median(axis=1))
# med_expr = pd.concat(med_expr, axis=1)
# med_expr.columns = tissues
# med_expr["gencodeId"] = med_expr.index
iqn = {}
for tissue in tissues:
    iqn_file = f"../data/expr/ensembl-gene_inv-quant_PCresiduals_{tissue}.txt"
    iqn[tissue] = pd.read_csv(iqn_file, sep="\t")

vcf = pysam.VariantFile("../data/P50.rnaseq.88.unpruned.vcf.gz")

# Load all significant pairs for singleTissueEqtl
# all_sig = pd.read_csv("../data/singleTissueEqtl.txt.gz", sep="\t")

exons = pd.read_csv("../data/exon.txt", sep="\t")

api = Flask(__name__)


@api.route("/tissueInfo", methods=["GET"])
def tissue_info():
    return open("../data/tissueInfo.json", "r").read()


@api.route("/topExpressedGene", methods=["GET"])
def top_expressed_gene():
    tissue = request.args.get("tissueSiteDetailId")
    return open(f"../data/topExpressedGene/{tissue}.json", "r").read()


@api.route("/gene", methods=["GET"])
def gene():
    ids = request.args.get("geneId").split(",")
    # info = [{"gencodeId": x, "geneSymbol": x} for x in ids]
    # d = genes.loc[genes["gencodeId"] == gene, :]
    # TEMP:
    ids = [x for x in ids if x in genes.index]
    d = genes.loc[ids, :].reset_index() # Include gencodeId in dict
    info = d.to_dict(orient="records")
    print(info)
    return json.dumps({"gene": info})


@api.route("/medianGeneExpression", methods=["GET"])
def med_gene_exp():
    # ids = request.args.get("gencodeId").split(",")
    # d = med_expr.loc[ids, :]
    # # a = d.to_numpy()
    # gene_clust = linkage(d, method="average", optimal_ordering=True)
    # d = d.melt(id_vars="gencodeId", var_name="tissueSiteDetailId", value_name="median")
    # d = d.to_dict(orient="records")
    # info = {"clusters": {"gene": gene_clust, "tissue": tissue_clust}, "medianGeneExpression": d}
    # return json.dumps()
    ## Instead, precompute for each tissue. edit BatchGeneExpression.js to pass tissue name instead of gene list.
    tissue = request.args.get("tissue")
    return open(f"../data/medianGeneExpression/{tissue}.json", "r").read()


@api.route("/variant", methods=["GET"])
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
    return json.dumps({"variant": infos})


@api.route("/dyneqtl", methods=["GET"])
def dyneqtl():
    variant = request.args.get("variantId")
    gene = request.args.get("gencodeId")
    tissue = request.args.get("tissueSiteDetailId")
    expr = iqn[tissue].loc[gene, :]
    rec = variant_record(variant)
    assert len(rec.alts) == 1, f"Multiple alt alleles: {variant}"
    gt = rec.samples
    indivs = [x.split("_")[0] for x in expr.index]
    geno = [genotype(gt[ind]["GT"]) for ind in indivs]
    # ignoring error, nes, tStatistic, timing
    counts = [int(np.sum(np.array(geno) == x)) for x in [0, 1, 2]]
    info = {
        "data": list(expr),
        "gencodeId": gene,
        "geneSymbol": "geneXYZ",
        "genotypes": geno,
        "hetCount": counts[1],
        "homoAltCount": counts[2],
        "homoRefCount": counts[0],
        "maf": (counts[1] + 2 * counts[2]) / int(2 * np.sum(counts)),
        "ref": rec.ref,  # I added ref and alt to API since they aren't in our variant IDs
        "alt": rec.alts[0],
        "pValue": "???",
        "pValueThreshold": "???",
        "tissueSiteDetailId": tissue,
        "variantId": variant,
    }
    return json.dumps(info)


@api.route("/singleTissueEqtl", methods=["GET"])
def single_tissue_eqtl():
    gene = request.args.get("gencodeId")
    # ignoring datasetId, snpId
    # chromosome, gencodeId, geneSymbol, geneSymbolUpper, nes, pValue, pos, tissueSiteDetailId, variantId
    # d = all_sig.loc[all_sig["gencodeId"] == gene, :].copy()
    d = pd.read_csv(f"../data/singleTissueEqtl/{gene}.txt.gz", sep="\t")
    d["geneSymbol"] = d["gencodeId"]
    d["geneSymbolUpper"] = d["gencodeId"]
    info = d.to_dict(orient="records")
    return json.dumps({"singleTissueEqtl": info})


@api.route("/ld", methods=["GET"])
def ld():
    gene = request.args.get("gencodeId")
    # d = all_sig.loc[all_sig["gencodeId"] == gene, :].copy()
    d = pd.read_csv(f"../data/singleTissueEqtl/{gene}.txt.gz", sep="\t")
    d["pos"] = [int(x.split(":")[1]) for x in d["variantId"]]
    d = d.sort_values(by="pos")
    ids = d["variantId"].unique()
    geno = geno_matrix(ids)
    ldmat = np.corrcoef(geno) ** 2
    lds = []
    for i in range(len(ids) - 1):
        for j in range(i + 1, len(ids)):
            lds.append([",".join([ids[i], ids[j]]), round(ldmat[i, j], 6)])
    return json.dumps({"ld": lds})


@api.route("/exon", methods=["GET"])
def exon():
    gene = request.args.get("gencodeId")
    # ignoring: gencodeVersion, genomeBuild
    # chromosome, end, exonId, exonNumber, featureType, gencodeId, geneSymbol, start, strand, transcriptId
    d = exons.loc[exons["gencodeId"] == gene, :]
    d = d.to_dict(orient="records")
    return json.dumps({"exon": d})


if __name__ == "__main__":
    api.run()
