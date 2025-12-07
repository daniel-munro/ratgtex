---
permalink: /
layout: base
---

<h1 class="text-center mt-4">
    <img src="/assets/images/RatGTExPortal.png" class="img-fluid" width="500px" alt="RatGTEx Portal">
</h1>

This portal provides a multimodal collection of RNA phenotypes (gene expression, alternative splicing, etc.), xQTLs for those phenotypes, and TWAS results for multiple rat tissues. It is managed by the [NIDA Center of Excellence for Genetics, Genomics, and Epigenetics of Substance Use Disorders in Outbred Rats (P30DA060810)](https://ratgenes.org).

### <span class="badge text-bg-success">Data Release v4</span>

The data and visualizations on this site are for the RatGTEx v4 release. The main changes from v3 are:

- Three tissue datasets from one new project have been added: a new brain region, agranular insular cortex (IC), and additional samples for existing brain regions, nucleus accumbens core (NAcc) and posterior ventral tegmental area (pVTA).
- [Pantry](https://pantry.pejlab.org) is now used to extract RNA phenotypes and perform xQTL mapping, providing data for six RNA modalities, compared to only gene expression levels and alternative splicing in the previous versions.
- The genome assembly has been updated from `mRatBN7.2` to `GRCr8`.

Data from previous versions are still available from the [Download](/download/) page.

### Tissues

<table class="table table-sm w-auto">
    <caption>* = Datasets from multiple studies were combined and analyzed as one larger dataset, with sample sizes of each constituent dataset indicated to the right. Sample sizes reflect the final post-QC datasets.</caption>
    <thead>
        <tr><th>Tissue</th><th>Abbreviation</th><th>Samples</th><th>Originating studies</th></tr>
    </thead>
    <tbody>
        <tr><td>Adipose</td><td>Adipose</td><td>411</td><td><a href="#R01DK106386">R01DK106386</a></td></tr>
        <tr><td>Agranular insular cortex</td><td>IC</td><td>148</td><td><a href="#P50DA037844-y6-10-proj1">P50DA037844 Y6-10 Project 1</a></td></tr>
        <tr><td>Basolateral amygdala</td><td>BLA</td><td>189</td><td><a href="#U01DA046077">U01DA046077</a></td></tr>
        <tr><td>Brain hemisphere</td><td>Brain</td><td>342</td><td><a href="#pilot-twas">Pilot: Creating the dataset for TWAS in HS rats</a></td></tr>
        <tr><td>Eye</td><td>Eye</td><td>53</td><td><a href="#R01EY021200">R01EY021200</a></td></tr>
        <tr><td>Infralimbic cortex</td><td>IL</td><td>83</td><td><a href="#P50DA037844-y1-5">P50DA037844 Y1-5</a></td></tr>
        <tr><td>Lateral habenula</td><td>LHb</td><td>82</td><td><a href="#P50DA037844-y1-5">P50DA037844 Y1-5</a></td></tr>
        <tr><td>Liver</td><td>Liver</td><td>411</td><td><a href="#R01DK106386">R01DK106386</a></td></tr>
        <tr><td>Nucleus accumbens core</td><td>NAcc</td><td>570*</td><td><a href="#P50DA037844-y1-5">P50DA037844 Y1-5</a> (77), <a href="#P50DA037844-y6-10-proj1">P50DA037844 Y6-10 Project 1</a> (143), <a href="#P50DA037844-y6-10-proj2">P50DA037844 Y6-10 Project 2</a> (159), <a href="#U01DA046077">U01DA046077</a> (191)</td></tr>
        <tr><td>Orbitofrontal cortex</td><td>OFC</td><td>82</td><td><a href="#P50DA037844-y1-5">P50DA037844 Y1-5</a></td></tr>
        <tr><td>Prelimbic cortex</td><td>PL</td><td>404*</td><td><a href="#P50DA037844-y1-5">P50DA037844 Y1-5</a> (82), <a href="#P50DA037844-y6-10-proj2">P50DA037844 Y6-10 Project 2</a> (130), <a href="#U01DA046077">U01DA046077</a> (192)</td></tr>
        <tr><td>Posterior ventral tegmental area</td><td>pVTA</td><td>293*</td><td><a href="#P50DA037844-y6-10-proj1">P50DA037844 Y6-10 Project 1</a> (141), <a href="#P50DA037844-y6-10-proj2">P50DA037844 Y6-10 Project 2</a> (152)</td></tr>
        <tr><td>Rostromedial tegmental nucleus</td><td>RMTg</td><td>92</td><td><a href="#U01DA044468">U01DA044468</a></td></tr>
    </tbody>
</table>

### Data Use

For the tissues above whose originating study has been published (see below), you may use the data in your studies by citing the originating study. For the other tissues, contact us about whether, when, and how you may use them.

### Associated Studies

The raw data for the tissues in this portal come from multiple studies. The processed data, however, may differ from the results reported in each study for a couple reasons:

1. Different processing decisions could be made for each study, while the data in this portal are processed uniformly to facilitate comparison across tissues.
2. The published study results are immutable, while new versions of data in this portal may be released, for example to improve quality or to use updated reference genome and annotations.

Nevertheless, we also host the original results from individual studies on the [Download](/download/#studies) page when available.

More info on the [ratgenes.org eQTL page](https://ratgenes.org/research-projects/eqtl/). Additional publications will be posted when available.

##### <a id="P50DA037844-y1-5" href="https://ratgenes.org/research-projects/rp2/">P50DA037844 Y1-5, Project 2</a>: Socially-acquired nicotine self-administration

- Publication: [The regulatory landscape of multiple brain regions in outbred heterogeneous stock rats](https://academic.oup.com/nar/article/50/19/10882/6764417)
- Investigator: Hao Chen, UTHSC
- Prior experience/treatment: Naive
- Sex ratio: 43 Female, 45 Male
- Age: Mean 86.3 days, SD 3.1
- Tissues: IL, LHb, NAcc, OFC, PL

##### <a id="P50DA037844-y6-10-proj1" href="https://ratgenes.org/research-projects/rp1new/">P50DA037844 Y6-10, Project 1</a>: Neurogenetic Substrates of Cocaine Addiction

- Investigator: Paul Meyer, University at Buffalo
- Prior experience/treatment: Exposure to cocaine
- Sex ratio: 74 Female, 75 Male
- Tissues: IC, NAcc, pVTA

##### <a id="P50DA037844-y6-10-proj2" href="https://ratgenes.org/research-projects/rp2/">P50DA037844 Y6-10, Project 2</a>: Socially-acquired nicotine self-administration

- Investigator: Hao Chen, UTHSC
- Prior experience/treatment: Naive
- Sex ratio: 84 Female, 81 Male
- Tissues: NAcc, PL, pVTA

##### <a id="pilot-twas" href="https://ratgenes.org/2020-pilot-grants-awarded/">Pilot: Creating the dataset for TWAS in HS rats</a>

- Investigator: Francesca Telese, UCSD
- Prior experience/treatment: Naive
- Sex ratio: 176 Female, 168 Male
- Tissues: Brain

##### <a id="R01DK106386" href="https://reporter.nih.gov/project-details/8941897">R01DK106386</a>: Systems genetics of adiposity traits in outbred rats

- Publication: [Genetic mapping of multiple traits identifies novel genes for adiposity, lipids, and insulin secretory capacity in outbred rats](https://diabetesjournals.org/diabetes/article-abstract/72/1/135/147723/Genetic-Mapping-of-Multiple-Traits-Identifies)
- Investigator: Leah Solberg Woods, WFU
- Prior experience/treatment: Naive, normal diet
- Sex ratio: All male
- Age: 17 weeks
- Tissues: Adipose, Liver

##### <a id="R01EY021200" href="https://reporter.nih.gov/project-details/10361394">R01EY021200</a>: Genetic Modulators of Glaucoma

- Publication: [Genome-wide association study finds multiple loci associated with intraocular pressure in HS rats](https://doi.org/10.3389/fgene.2022.1029058)
- Investigator: Monica Jablonsky, UTHSC
- Prior experience/treatment: Behavioral testing including brief nicotine exposure
- Sex ratio: 29 Female, 24 Male
- Tissues: Eye

##### <a id="U01DA044468" href="https://ratgenes.org/research-projects-u01da044468/">U01DA044468</a>: Genomic analysis of avoidance learning in addiction

- Investigator: Tom Jhou, University of Maryland
- Prior experience/treatment: Extreme high and low phenotypes in cocaine-avoidance and punishment resistance
- Sex ratio: 48 Female, 45 Male
- Tissues: RMTg

##### <a id="U01DA046077" href="https://ratgenes.org/research-projects-u01da046077/">U01DA046077</a>: Identification of Genetic Features of Delay Discounting Using a Heterogeneous Stock Rat Model

- Investigators: Suzanne Mitchell and Robert Hitzemann, OHSU
- Prior experience/treatment: Behavior testing - delay discounting
- Sex ratio: 100 Female, 100 Male
- Tissues: BLA, NAcc, PL

### Visualizations

- [Expression Heatmap for top expressed genes](/top-expressed/)
- [Expression Heatmap for a query gene list](/query-gene-expression/)
- [eQTL Dashboard](/eqtl-dashboard/)
- [Gene-eQTL Visualizer](/gene-eqtl-visualizer/)

Also try searching for a gene in the toolbar to view its expression and eQTLs.

### Contact {#contact}

For questions or comments related to the data or website, email Daniel Munro at [dmunro@health.ucsd.edu](mailto:dmunro@health.ucsd.edu).

RatGTEx is primarily developed in the labs of [Abraham Palmer at UC San Diego](https://palmerlab.org/) and [Pejman Mohammadi at Seattle Children's Research Institute and University of Washington](https://pejlab.org/).

### Acknowledgments

RatGTEx is produced by the [NIDA Center of Excellence for Genetics, Genomics, and Epigenetics of Substance Use Disorders in Outbred Rats (P30DA060810)](https://ratgenes.org). It is not a project of the GTEx Consortium.

Code for visualizations was adapted from the [GTEx portal viz tools](https://github.com/broadinstitute/gtex-viz), and additional design considerations were inspired by the [GTEx portal](https://gtexportal.org/).

### More info

[RRID:SCR_022145](https://scicrunch.org/resources/data/record/nlx_144509-1/SCR_022145/resolver)

See also [PhenoGen Informatics](https://phenogen.org/index.jsp), the site for quantitative genetics of the transcriptome.
