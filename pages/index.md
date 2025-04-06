---
permalink: /
layout: base
---

<h1 class="text-center mt-4">
    <img src="/assets/images/RatGTExPortal.png" class="img-fluid" width="500px" alt="RatGTEx Portal">
</h1>

This portal provides gene expression, eQTL, and sQTL data for multiple rat tissues. It is managed by the [NIDA Center of Excellence for Genetics, Genomics, and Epigenetics of Substance Use Disorders in Outbred Rats (P30DA060810)](https://ratgenes.org).

### Tissues

The raw data for the tissues in this portal come from multiple studies. The processed data, however, may differ from the results reported in each study for a couple reasons:

1. Different processing decisions could be made for each study, while the data in this portal are processed uniformly to facilitate comparison across tissues.
2. The published study results are immutable, while new versions of data in this portal may be released, for example to improve quality or to use updated reference genome and annotations.

Nevertheless, we also host the original results from individual studies on the [Download](/download/#studies) page when available.

<table class="table table-sm w-auto">
    <caption>Sample sizes reflect the final post-QC datasets.</caption>
    <thead>
        <tr><th>Tissue</th><th>Abbreviation</th><th>Samples</th><th>Originating studies (see below)</th></tr>
    </thead>
    <tbody>
        <tr><td>Adipose</td><td>Adipose</td><td>411</td><td>R01DK106386</td></tr>
        <tr><td>Basolateral amygdala</td><td>BLA</td><td>191</td><td>U01DA046077</td></tr>
        <tr><td>Brain hemisphere</td><td>Brain</td><td>341</td><td>Pilot: Creating the dataset for TWAS in HS rats</td></tr>
        <tr><td>Eye</td><td>Eye</td><td>53</td><td>R01EY021200</td></tr>
        <tr><td>Infralimbic cortex</td><td>IL</td><td>83</td><td>P50DA037844 Y1-5</td></tr>
        <tr><td>Lateral habenula</td><td>LHb</td><td>82</td><td>P50DA037844 Y1-5</td></tr>
        <tr><td>Liver</td><td>Liver</td><td>411</td><td>R01DK106386</td></tr>
        <tr><td>Nucleus accumbens core</td><td>NAcc</td><td>430</td><td>P50DA037844 Y1-5, U01DA046077, P50DA037844 Y5-10</td></tr>
        <tr><td>Orbitofrontal cortex</td><td>OFC</td><td>82</td><td>P50DA037844 Y1-5</td></tr>
        <tr><td>Prelimbic cortex</td><td>PL</td><td>407</td><td>P50DA037844 Y1-5, U01DA046077, P50DA037844 Y5-10</td></tr>
        <tr><td>Posterior ventral tegmental area</td><td>pVTA</td><td>153</td><td>P50DA037844 Y5-10</td></tr>
        <tr><td>Rostromedial tegmental nucleus</td><td>RMTg</td><td>93</td><td>U01DA044468</td></tr>
    </tbody>
</table>

### Data Use

For the tissues above whose originating study has been published (see below), you may use the data in your studies by citing the originating study. For the other tissues, contact us about whether, when, and how you may use them.

### Associated Studies

More info on the [ratgenes.org eQTL page](https://ratgenes.org/research-projects/eqtl/). Additional publications will be posted when available.

Tissues from the same study are from the same set of rats, though the final subsets per tissue differ slightly after quality control filtering. See the [sample info page](/about/samples/) for details.

##### P50DA037844 Y1-5, Project 2: Socially-acquired nicotine self-administration

- Publication: [The regulatory landscape of multiple brain regions in outbred heterogeneous stock rats](https://academic.oup.com/nar/article/50/19/10882/6764417)
- Investigator: Hao Chen, UTHSC
- Prior experience/treatment: Naive
- Sex ratio: 43 Female, 45 Male
- Age: Mean 86.3 days, SD 3.1

##### R01EY021200: Genetic Modulators of Glaucoma

- Publication: [Genome-wide association study finds multiple loci associated with intraocular pressure in HS rats](https://doi.org/10.3389/fgene.2022.1029058)
- Investigator: Monica Jablonsky, UTHSC
- Sex ratio: 29 Female, 24 Male

##### U01DA046077: Identification of Genetic Features of Delay Discounting Using a Heterogeneous Stock Rat Model

- Investigators: Suzanne Mitchell and Robert Hitzemann, OSHU
- Prior experience/treatment: Behavior testing - delay discounting
- Sex ratio: 100 Female, 100 Male

##### R01DK106386: Systems genetics of adiposity traits in outbred rats

- Publication: [Genetic mapping of multiple traits identifies novel genes for adiposity, lipids, and insulin secretory capacity in outbred rats](https://diabetesjournals.org/diabetes/article-abstract/72/1/135/147723/Genetic-Mapping-of-Multiple-Traits-Identifies)
- Investigator: Leah Solberg Woods, WFU
- Prior experience/treatment: Naive, normal diet
- Sex ratio: All male
- Age: 17 weeks

##### Pilot: Creating the dataset for TWAS in HS rats

- Investigator: Francesca Telese, UCSD
- Prior experience/treatment: Naive  
- Sex ratio: 176 Female, 168 Male

##### U01DA044468: Genomic analysis of avoidance learning in addiction

- Investigator: Tom Jhou, University of Maryland
- Prior experience/treatment: Extreme high and low phenotypes in cocaine-avoidance and punishment resistance
- Sex ratio: 48 Female, 45 Male

##### P50DA037844 Y5-10, Project 2: Socially-acquired nicotine self-administration

- Investigator: Hao Chen, UTHSC
- Prior experience/treatment: Naive
- Sex ratio: 84 Female, 81 Male

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
