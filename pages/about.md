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

#### Processing

The [ratgtex-pipeline](https://github.com/daniel-munro/ratgtex-pipeline) repository contains all RNA-Seq and eQTL/sQTL mapping pipeline code. The steps are also on [protocols.io](http://dx.doi.org/10.17504/protocols.io.rm7vzyk92lx1/v1). The [ratgtex-server-data](https://github.com/daniel-munro/ratgtex-server-data) repository contains code that processes those results into additional data files for download and for the API.

- Genome assembly and gene annotations used for each data release:
  - v1:
    - Assembly: [Rnor_6.0 Ensembl](http://ftp.ensembl.org/pub/release-99/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz) (GCA_000001895.4)
    - Annotations: [Ensembl Rat Release 99](http://ftp.ensembl.org/pub/release-99/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.99.gtf.gz)
  - v2:
    - Assembly: [mRatBN7.2 Ensembl](http://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz) (GCA_015227675.2)
    - Annotations: [Ensembl Rat Release 108](http://ftp.ensembl.org/pub/release-108/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.108.gtf.gz)
  - v3:
    - Assembly: [mRatBN7.2 RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz) (GCA_015227675.2)
    - Annotations: [mRatBN7.2 RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz)
- RNA-Seq read counts, transformed as log<sub>2</sub>(count+1), were used to compute cis-eQTL effect sizes (allelic fold change).
- Inverse-quantile normalized expression values were used for eQTL mapping.
- eQTL mapping: [tensorQTL](https://github.com/broadinstitute/tensorqtl). cis-eQTL mapping tested SNPs in cis-windows +/- 1 Mb from each gene's transcription start site.
- sQTL mapping: splice junction usage was quantified using leafcutter, and sQTLs were mapped with tensorQTL, grouping phenotypes by gene. cis-sQTL mapping tested SNPs in cis-windows +/- 1 Mb from each gene's transcription start site using, with multiple splice phenotypes per gene tested as a group.

#### Formats

##### Gene expression (`expr.{units}.{tissue}.{version}.bed.gz`)

Expression tables are provided in BED format, with four columns describing the gene plus one column per sample. Only the transcription start sites are specified, not the whole gene length.

- **#chr** Chromosome name, e.g. `chr1`
- **start** Transcription start site. Coordinates in BED format are 0-based, so this is one less than the position as it is typically represented.
- **end** One higher than the start value. Thus, this would be the typical coordinate to represent the TSS.
- **gene_id** Gene symbol
- **{sample}** Each additional column is labeled with the sample or rat ID and contains log2(count+1), TPM, or inverse-quantile normalized values.

##### All significant cis-eQTL gene-SNP pairs (`cis_qtl_signif.{tissue}.{version}.txt.gz`)

A compressed tab-separated table.

- **gene_id** Gene symbol
- **variant_id** SNP ID (e.g. `chr1:669562`)
- **tss_distance** SNP position - TSS position, oriented on the gene's strand
- **af** Alternative allele frequency
- **ma_samples** Number of samples with minor allele
- **ma_count** Total number of minor alleles (i.e. no. het. + 2 \* no. hom. alt.)
- **pval_nominal** Nominal p-value
- **slope** Coefficient for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **pval_nominal_threshold** Gene-tissue-specific nominal p-value significance threshold as determined by permutations and transcriptome-wide FDR

##### Strong associations from trans-eQTL mapping (`trans_qtl_pairs.{tissue}.{version}.txt.gz`)

All measured SNPs genome-wide were tested against expression of each gene, and pairs with p-value &lt; 1e-5 and &gt; 5 Mb TSS distance are included here. Rows are sorted by variant location.

- **variant_id** SNP ID (e.g. `chr1:669562`)
- **gene_id** Gene symbol
- **pval** Nominal p-value
- **slope** Coefficient (slope) for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **af** Alternative allele frequency

##### Conditionally independent cis-eQTLs (`eqtls_indep.{version}.txt`)

A table of cis-eQTLs from all tissues. Stepwise regression was used to test for cis-eQTLs beyond and uncorrelated with the top cis-eQTL per gene.

- **tissue** Tissue abbreviation
- **gene_id** Gene symbol
- **num_var** Number of cis-window variants tested for the specified gene in this tissue
- **variant_id** SNP ID (e.g. `chr1:669562`) for the top eQTL SNP (eSNP). In cases of tied top SNPs that are in 100% LD, one was randomly chosen among them.
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
- **rank** The strongest cis-eQTL found per gene-tissue has rank 1, next-strongest has rank 2, etc.
- **log2_aFC** eQTL effect size measured as allelic fold change

##### Top association per gene (`top_assoc.{version}.txt`)

A table of the strongest variant-gene association per gene per tissue, even if not significant.

- **tissue** Tissue abbreviation
- **gene_id** Gene symbol
- **num_var** Number of cis-window variants tested for the specified gene in this tissue
- **variant_id** SNP ID (e.g. `chr1:669562`) for the top eQTL SNP (eSNP). In cases of tied top SNPs that are in 100% LD, one was randomly chosen among them.
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
- **log2_aFC** Effect size measured as allelic fold change

<!-- Splicing -->

##### Splice phenotypes (`leafcutter.{tissue}.{version}.bed.gz`)

Splice phenotypes are provided in BED format, with four columns describing the splice junctions and genes plus one column per sample. Only the gene transcription start sites are specified, not the whole gene length or the splice junction itself, as these are used to determine the gene's cis-window for QTL mapping. See [leafcutter docs](http://davidaknowles.github.io/leafcutter/) for descriptions of these junctions and junction clusters.

- **#Chr** Chromosome name, e.g. `chr1`
- **start** Transcription start site. Coordinates in BED format are 0-based, so this is one less than the position as it is typically represented.
- **end** One higher than the start value. Thus, this would be the typical coordinate to represent the TSS.
- **ID** Splice junction ID (`chromosome:start:end:cluster_id:gene_id`)
- **{sample}** Each additional column is labeled with the sample or rat ID and contains inverse-quantile normalized values.

##### All significant cis-sQTL phenotype-SNP pairs (`cis_qtl_signif.{tissue}.{version}.txt.gz`)

A compressed tab-separated table.

- **phenotype_id** Splice junction ID (`chromosome:start:end:cluster_id:gene_id`)
- **variant_id** SNP ID (e.g. `chr1:669562`)
- **tss_distance** SNP position - TSS position, oriented on the gene's strand
- **af** Alternative allele frequency
- **ma_samples** Number of samples with minor allele
- **ma_count** Total number of minor alleles (i.e. no. het. + 2 \* no. hom. alt.)
- **pval_nominal** Nominal p-value
- **slope** Coefficient for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **pval_nominal_threshold** Gene-tissue-specific nominal p-value significance threshold as determined by permutations and transcriptome-wide FDR

##### Strong associations from trans-sQTL mapping (`splice.trans_qtl_pairs.{tissue}.{version}.txt.gz`)

All measured SNPs genome-wide were tested against each splice phenotype, and pairs with p-value &lt; 1e-5 and &gt; 5 Mb TSS distance are included here. Rows are sorted by variant location.

- **variant_id** SNP ID (e.g. `chr1:669562`)
- **phenotype_id** Splice junction ID (`chromosome:start:end:cluster_id:gene_id`)
- **pval** Nominal p-value
- **slope** Coefficient (slope) for the SNP genotype in the linear model
- **slope_se** Standard error of the slope
- **af** Alternative allele frequency

##### Conditionally independent cis-sQTLs (`sqtls_indep.{version}.txt`)

A table of cis-sQTLs from all tissues. Stepwise regression was used to test for cis-sQTLs beyond and uncorrelated with the top cis-sQTL per gene.

- **tissue** Tissue abbreviation
- **phenotype_id** Splice junction ID (`chromosome:start:end:cluster_id:gene_id`)
- **num_var** Number of cis-window variants tested for the specified phenotype in this tissue
- **variant_id** SNP ID (e.g. `chr1:669562`) for the top sQTL SNP (sSNP). In cases of tied top SNPs that are in 100% LD, one was randomly chosen among them.
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
- **gene_id** Gene symbol, used as phenotype groups for sQTL mapping
- **group_size** Number of splice phenotypes tested for the given gene
- **rank** The strongest cis-sQTL found per gene-tissue has rank 1, next-strongest has rank 2, etc.

##### Top splice association per gene (`top_assoc_splice.{version}.txt`)

A table of the strongest variant-gene association per gene per tissue, even if not significant.

- **tissue** Tissue abbreviation
- **phenotype_id** Splice junction ID (`chromosome:start:end:cluster_id:gene_id`)
- **num_var** Number of cis-window variants tested for the specified phenotype in this tissue
- **variant_id** SNP ID (e.g. `chr1:669562`) for the top sQTL SNP (sSNP). In cases of tied top SNPs that are in 100% LD, one was randomly chosen among them.
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
- **gene_id** Gene symbol, used as phenotype groups for sQTL mapping
- **group_size** Number of splice phenotypes tested for the given gene
- **qval** Q-value used to control transcriptome-wide FDR per tissue (q < 0.05 was used to determine significant cis-sQTLs)
- **pval_nominal_threshold** Gene-tissue-specific nominal p-value significance threshold as determined by permutations and transcriptome-wide FDR
