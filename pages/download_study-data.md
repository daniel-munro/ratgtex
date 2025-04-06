---
permalink: /download/study-data/
title: Download study data
layout: base
---

# Original study data

The processed data on the main Download page may differ from the results reported in each original study for a couple reasons:

1. Different processing decisions could be made for each study, while the data in this portal are processed uniformly to facilitate comparison across tissues.
2. The published study results are immutable, while new versions of data in this portal may be released, for example to improve quality or to use updated reference genome and annotations.

Therefore, we also provide data files here corresponding directly to the results in each paper.

## IL, LHb, NAcc1, OFC, PL1 {#IL_LHb_NAcc_OFC_PL}

[The regulatory landscape of multiple brain regions in outbred heterogeneous stock rats](https://academic.oup.com/nar/article/50/19/10882/6764417)

The `Rnor_6.0` genome was used for all data here.

- [`P50.rnaseq.88.unpruned.vcf.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/P50.rnaseq.88.unpruned.vcf.gz)
- [`metadata_p50_hao_chen_2014.csv`](/data/studies/IL_LHb_NAcc_OFC_PL/metadata_p50_hao_chen_2014.csv)
- [`ensembl-gene_raw-counts.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_raw-counts.txt)
- [`ensembl-gene_raw-tpm.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_raw-tpm.txt)
- [`ensembl-gene_log2_IL.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_log2_IL.bed.gz)
- [`ensembl-gene_log2_LHb.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_log2_LHb.bed.gz)
- [`ensembl-gene_log2_NAcc.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_log2_NAcc.bed.gz)
- [`ensembl-gene_log2_OFC.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_log2_OFC.bed.gz)
- [`ensembl-gene_log2_PL.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_log2_PL.bed.gz)
- [`ensembl-gene_inv-quant_IL.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_inv-quant_IL.bed.gz)
- [`ensembl-gene_inv-quant_LHb.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_inv-quant_LHb.bed.gz)
- [`ensembl-gene_inv-quant_NAcc.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_inv-quant_NAcc.bed.gz)
- [`ensembl-gene_inv-quant_OFC.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_inv-quant_OFC.bed.gz)
- [`ensembl-gene_inv-quant_PL.bed.gz`](/data/studies/IL_LHb_NAcc_OFC_PL/ensembl-gene_inv-quant_PL.bed.gz)
- [`top_assoc.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/top_assoc.txt)
- [`eqtls_indep.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/eqtls_indep.txt)
- [`ASE_aFC.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/ASE_aFC.txt)
- [`top_assoc_splice.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/top_assoc_splice.txt)
- [`sqtls_indep.txt`](/data/studies/IL_LHb_NAcc_OFC_PL/sqtls_indep.txt)

## Eye {#Eye}

[Genome-wide association study finds multiple loci associated with intraocular pressure in HS rats](https://www.frontiersin.org/articles/10.3389/fgene.2022.1029058/full)

The `Rnor_6.0` genome was used for all data here.

- [`eyes.vcf.gz`](/data/studies/Eye/eyes.vcf.gz)
- [`Eye.expr.log2.bed.gz`](/data/studies/Eye/Eye.expr.log2.bed.gz)
- [`Eye.expr.tpm.bed.gz`](/data/studies/Eye/Eye.expr.tpm.bed.gz)
- [`Eye.expr.iqn.bed.gz`](/data/studies/Eye/Eye.expr.iqn.bed.gz)
- [`Eye.expr.iqn.filtered.bed.gz`](/data/studies/Eye/Eye.expr.iqn.filtered.bed.gz)
- [`covar.txt`](/data/studies/Eye/covar.txt)
- [`Eye.cis_qtl_signif.txt.gz`](/data/studies/Eye/Eye.cis_qtl_signif.txt.gz)
- [`Eye.cis_qtl.txt.gz`](/data/studies/Eye/Eye.cis_qtl.txt.gz)
- [`Eye.cis_independent_qtl.txt.gz`](/data/studies/Eye/Eye.cis_independent_qtl.txt.gz)
- [`Eye.top_assoc.txt`](/data/studies/Eye/Eye.top_assoc.txt)
- [`Eye.top_eQTLs.txt`](/data/studies/Eye/Eye.top_eQTLs.txt)
- [`Eye.aFC.txt`](/data/studies/Eye/Eye.aFC.txt)

## Adipose, Liver {#Adipose_Liver}

[Genetic mapping of multiple traits identifies novel genes for adiposity, lipids, and insulin secretory capacity in outbred rats](https://diabetesjournals.org/diabetes/article-abstract/72/1/135/147723/Genetic-Mapping-of-Multiple-Traits-Identifies)

Phenotype and genotype data have been deposited in [RGD](https://rgd.mcw.edu/) (`RGD:153344611`).

# Additional studies

These studies used RatGTEx data and/or produced data that is relevant but hasn't been integrated into the RatGTEx portal.

### [Y and mitochondrial chromosomes in the heterogeneous stock rat population](https://doi.org/10.1093/g3journal/jkae213) {#Y_MT}

This study used the following files from RatGTEx (all are from v2 using `mRatBN7.2`):

- [`BLA.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/BLA.rn7.expr.log2.bed.gz)
- [`Brain.rn7.expr.iqn.bed.gz`](/data/studies/Y_MT/Brain.rn7.expr.iqn.bed.gz)
- [`Brain.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/Brain.rn7.expr.log2.bed.gz)
- [`Brain.rn7.expr.tpm.bed.gz`](/data/studies/Y_MT/Brain.rn7.expr.tpm.bed.gz)
- [`Eye.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/Eye.rn7.expr.log2.bed.gz)
- [`IL.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/IL.rn7.expr.log2.bed.gz)
- [`LHb.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/LHb.rn7.expr.log2.bed.gz)
- [`NAcc.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/NAcc.rn7.expr.log2.bed.gz)
- [`NAcc2.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/NAcc2.rn7.expr.log2.bed.gz)
- [`OFC.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/OFC.rn7.expr.log2.bed.gz)
- [`PL.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/PL.rn7.expr.log2.bed.gz)
- [`PL2.rn7.expr.log2.bed.gz`](/data/studies/Y_MT/PL2.rn7.expr.log2.bed.gz)

### [Bioenergetic-related gene expression in the hippocampus predicts internalizing vs. externalizing behavior in an animal model of temperament](https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2025.1469467/full) {#HPC_F2}

This study produced cis-eQTLs for hippocampus tissue from an F2 cross of selectively-bred rats using the `Rnor_6.0` genome.

- [`f2.eqtls_indep.txt`](/data/studies/HPC_F2/f2.eqtls_indep.txt) - conditionally independent cis-eQTLs
- [`f2.top_assoc.txt`](/data/studies/HPC_F2/f2.top_assoc.txt) - top association per gene, even if not significant  
- [`f2.cis_qtl_signif.txt.gz`](/data/studies/HPC_F2/f2.cis_qtl_signif.txt.gz) - all significantly associated cis-window SNPs
- [`colocs.tsv`](/data/studies/HPC_F2/colocs.tsv) - stats for all colocalization tests
- [`colocs_sig.tsv`](/data/studies/HPC_F2/colocs_sig.tsv) - stats for significant colocalization tests
