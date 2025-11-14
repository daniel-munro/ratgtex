---
permalink: /about/
title: About
layout: base
---

# Data Specifications {#data-specs}

The datasets provided in this portal conform to the following specifications to facilitate comparisons across tissues.
The raw data for the tissues in this portal come from multiple studies. The processed data, however, may differ from the results reported in each study for a couple reasons:

1. Different processing decisions could be made for each study, while the data in this portal are processed uniformly to facilitate comparison across tissues.
2. The published study results are immutable, while new versions of data in this portal may be released, for example to improve quality or to use updated reference genome and annotations.

Nevertheless, we also host the original results from individual studies on the [Download](/download/#studies) page when available.

### Processing

The [ratgtex-pipeline](https://github.com/daniel-munro/ratgtex-pipeline) repository contains code for preprocessing, QC, read mapping, and some QTL mapping code. [Pantry](https://github.com/PejLab/Pantry) was used to quantify RNA phenotypes and map xQTLs. The original RatGTEx pipeline is on [protocols.io](http://dx.doi.org/10.17504/protocols.io.rm7vzyk92lx1/v1) (this does not yet reflect the latest version that uses Pantry). The [ratgtex-server-data](https://github.com/daniel-munro/ratgtex-server-data) repository contains code that processes results into additional data files for download and for the API.

##### Genome assembly and gene annotations used for each data release

- v1:
  - Assembly: [Rnor_6.0 Ensembl](http://ftp.ensembl.org/pub/release-99/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz) (GCA_000001895.4)
  - Annotations: [Ensembl Rat Release 99](http://ftp.ensembl.org/pub/release-99/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.99.gtf.gz)
- v2:
  - Assembly: [mRatBN7.2 Ensembl](http://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz) (GCA_015227675.2)
  - Annotations: [Ensembl Rat Release 108](http://ftp.ensembl.org/pub/release-108/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.108.gtf.gz)
- v3:
  - Assembly: [mRatBN7.2 RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz) (GCA_015227675.2)
  - Annotations: [mRatBN7.2 RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz)
- v4:
  - Assembly [GRCr8 RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz) (GCA_036323735.1)
  - Annotations: [GRCr8 RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.gtf.gz)

##### Methods used for v4

- RNA-seq reads mapped using [STAR](https://pubmed.ncbi.nlm.nih.gov/23104886/) 2.7.11b
- Phenotyping:
  - TSS and polyA site usage quantified using [txrevise](https://pubmed.ncbi.nlm.nih.gov/30618377/) and [kallisto](https://pubmed.ncbi.nlm.nih.gov/27043002/)
  - Expression levels and isoform ratios quantified using [kallistpo](https://pubmed.ncbi.nlm.nih.gov/27043002/)
  - Intron excision ratios (alternative splicing) quantified using [RegTools](https://pubmed.ncbi.nlm.nih.gov/36949070/) and [LeafCutter](https://pubmed.ncbi.nlm.nih.gov/29229983/)
  - RNA stability quantified using [featureCounts](https://pubmed.ncbi.nlm.nih.gov/24227677/) to count exonic and intronic reads and using their ratio as a [model of mRNA stability](https://pubmed.ncbi.nlm.nih.gov/26098447/)
- xQTL mapping:
  - All RNA phenotypes were normalized by quantile normalization followed by inverse-normal transform. Those with more than 50% zeros across samples before normalization were removed.
  - Data-driven covariates were computed for use in xQTL mapping: The first five principal components from the genotype matrix and the first 20 principal components from the RNA phenotype matrix
  - xQTL mapping using [tensorQTL](https://pubmed.ncbi.nlm.nih.gov/31675989/)
    - cis-QTL mapping was performed on each modality of RNA phenotypes, and was also performed on the combined set of phenotypes (cross-modality mapping), which produced cis-xQTLs that are conditionally independent across modalities per gene per tissue
    - eQTL results in this portal are from the expression-only mapping, and both separate and cross-modality mapping results are available for download
  - [aFC](https://pubmed.ncbi.nlm.nih.gov/29021289/) (allelic fold change) was used to quantify the effect size of cis-eQTLs

### Formats

##### RNA phenotype tables (including gene expression)

`expr.{log2,tpm}.{tissue}.v4.bed.gz`
<br />
`phenos.{tissue}.{modality}.{unnorm,norm}.v4.bed.gz`

BED format, with four phenotype information columns plus one column per sample. The coordinates define the cis-window around the gene TSS for xQTL mapping.

- **#chr** Chromosome name, e.g. `chr1`
- **start** Transcription start site of the phenotype's gene. Coordinates in BED format are 0-based, so this is one less than the position as it is typically represented.
- **end** One higher than the start value. Thus, this would be the typical coordinate to represent the TSS.
- **phenotype_id** For modalities with one phenotype per gene (expression and stability), this is the gene symbol. Otherwise, it is an ID with format `{gene_id}__{other_phenotype_info}`
- **{sample}** Each additional column is labeled with the sample or rat ID and contains the raw or normalized phenotype values.

##### All significant cis-xQTL phenotype-SNP pairs

`cis_qtl_signif.{tissue}.modality.v4_rn8.txt.gz`

A compressed tab-separated table.

- **phenotype_id** Phenotype ID, see above
- **variant_id** SNP ID (e.g. `chr1:793394`)
- **tss_distance** SNP position - TSS position, oriented on the gene's strand
- **af** Alternative allele frequency
- **ma_samples** Number of samples with minor allele
- **ma_count** Total number of minor alleles (i.e. # het. + 2 \* # hom. alt.)
- **pval_nominal** Nominal p-value
- **slope** Coefficient for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **pval_nominal_threshold** Gene-tissue-specific nominal p-value significance threshold as determined by permutations and transcriptome-wide FDR 0.05

##### Strong associations from trans-eQTL mapping

`trans_qtl_pairs.{tissue}.v4.txt.gz`

All measured SNPs genome-wide were tested against expression of each gene, and pairs with p-value &lt; 1e-5 and &gt; 5 Mb TSS distance are included here. Rows are sorted by variant location.

- **variant_id** SNP ID (e.g. `chr1:793394`)
- **phenotype_id** Phenotype ID, see above
- **pval** Nominal p-value
- **slope** Coefficient (slope) for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **af** Alternative allele frequency

##### Conditionally independent cis-xQTLs

`eqtls_indep.v4_rn8.tsv`
<br />
`xqtls_indep.{modality}.v4_rn8.tsv`

A table of cis-QTLs from all tissues. Stepwise regression was used to test for cis-QTLs beyond and uncorrelated with the top cis-QTL per gene. `eqtls_indep.v4_rn8.tsv` and `xqtls_indep.cross_modality.v4_rn8.tsv` are the same results but with formatting differences as described below. For `xqtls_indep.cross_modality.v4_rn8.tsv`, all six modalities were combined and mapped as one set so that xQTLs are conditionally independent across modalities per gene per tissue.

- **tissue** Tissue abbreviation
- **modality** Modality of RNA phenotype. Present in `xqtls_indep.cross_modality.v4_rn8.tsv` only.
- **phenotype_id** Phenotype ID, see above. Not present in `eqtls_indep.v4_rn8.tsv`.
- **gene_id** Gene symbol
- **num_var** Number of cis-window variants tested for the specified gene in this tissue
- **variant_id** SNP ID (e.g. `chr1:793394`) for the top eQTL SNP (eSNP). In cases of tied top SNPs that are in 100% LD, one was randomly chosen among them.
- **chrom** Chromosome of the top SNP (e.g. `chr1`)
- **pos** Position of the top SNP (bp)
- **ref** Reference allele for the top SNP
- **alt** Alternative allele for the top SNP
- **af** Alternative allele frequency
- **tss_distance** SNP position - TSS position, oriented on the gene's strand
- **pval_nominal** Nominal p-value
- **slope** Coefficient for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **pval_beta** Empirical beta-approximated p-value
- **rank** The strongest cis-QTL found per gene per tissue has rank 1, next-strongest has rank 2, etc.
- **log2_aFC** eQTL effect size measured as allelic fold change. Present in `eqtls_indep.v4_rn8.tsv` only.

##### Top association per gene

`top_assoc.v4_rn8.tsv`
<br />
`top_assoc.{modality}.v4_rn8.tsv`

A table of the strongest variant-gene association per gene per tissue, even if not significant. `top_assoc.v4_rn8.tsv` and `top_assoc.expression.v4_rn8.tsv` are the same results but with formatting differences as described below. For `top_assoc.cross_modality.v4_rn8.txt`, all six modalities were combined and mapped as one set so that one top association of any modality per gene per tissue is shown.

- **tissue** Tissue abbreviation
- **modality** Modality of RNA phenotype. Present in `top_assoc.cross_modality.v4_rn8.tsv` only.
- **phenotype_id** Phenotype ID, see above. Not present in `top_assoc.v4_rn8.tsv`.
- **gene_id** Gene symbol
- **num_var** Number of cis-window variants tested for the specified gene in this tissue
- **variant_id** SNP ID (e.g. `chr1:793394`) for the top eQTL SNP (eSNP). In cases of tied top SNPs that are in 100% LD, one was randomly chosen among them.
- **chrom** Chromosome of the top SNP (e.g. `chr1`)
- **pos** Position of the top SNP (bp)
- **ref** Reference allele for the top SNP
- **alt** Alternative allele for the top SNP
- **af** Alternative allele frequency
- **tss_distance** SNP position - TSS position, oriented on the gene's strand
- **pval_nominal** Nominal p-value
- **slope** Coefficient for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **pval_beta** Empirical beta-approximated p-value
- **qval** Q-value used to control transcriptome-wide FDR per tissue (q < 0.05 was used to determine significant cis-eQTLs)
- **pval_nominal_threshold** Gene-tissue-specific nominal p-value significance threshold as determined by permutations and transcriptome-wide FDR
- **log2_aFC** eQTL effect size measured as allelic fold change. Present in `top_assoc.v4_rn8.tsv` only.
