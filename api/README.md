# RatGTEx API

The API serves RatGTEx v4 data for the GRCr8 assembly. Run commands from the repository root:

```sh
python -m pip install -r requirements.txt
python -m api.app
python -m unittest discover -s tests -p 'test_*.py' -v
```

`app:api` is the Gunicorn entry point. `RATGTEX_DATA_DIR` selects the v4 release directory; its default is this repository's `data/v4`, independently of the working directory. API schema versions describe the interface independently of scientific data releases. The public contract is [openapi.yaml](../openapi.yaml), rendered at `/about/api/`; `/metadata` lists modalities, tissues, cohorts, endpoints and request limits.

## Phenotype and association queries

`/genotypePhenotypes` returns normalized RNA phenotype values and nominal cis associations for one gene, variant, tissue and modality:

```text
/api/v4/genotypePhenotypes?geneId=Ifnlr1&variantId=chr5:153176222&tissueId=IL&modality=splicing
```

`modality` defaults to `expression`. The same request and response structure applies to expression, isoform ratios, intron excision ratios, alternative TSS usage, alternative polyA usage and RNA stability. The `genotypePhenotypes` array contains every available phenotype belonging to the queried gene and modality, including non-significant phenotypes. Each item contains `phenotypeId`, `data`, nominal `pValue`, normalized-phenotype `slope`, `slopeSE` and testing status. Optional `phenotypeId` selects one phenotype from that gene. The shared `genotypes` array and each phenotype's `data` vector use matching sample order. `includeSampleIds=true` includes the aligned sample IDs.

`/phenotypes?geneId=Ifnlr1` provides a lightweight inventory of measured phenotypes across tissues and modalities, without reading sample vectors or genotypes. Optional tissue and modality filters narrow the inventory. Its `datasets` list distinguishes available datasets with no measurements from unavailable datasets; the latter do not prevent other datasets from being listed. The dashboard uses this inventory to align rows before progressively loading plot values. No additional prepared data are required.

`/geneExpression` provides expression-only queries for up to ten genes, with optional tissue filtering. `expressionScale=tpm` (default) returns TPM; `expressionScale=normalized` returns the released quantile-normalized, inverse-normal transformed expression. `/expressionSummary` summarizes the same scales without sending individual values. Displayed vectors receive no additional covariate adjustment. Missing gene/tissue measurements are omitted. `datasetId` identifies the data release; `genotypeDatasetId` identifies the cohort.

`/eqtl` returns published independent expression associations. `/xqtl` returns published independent cross-modality associations and accepts a modality filter. `/significantCisEqtls` returns significant cis-expression associations, with the collection key `significantCisEqtls`. These endpoints report the published mapping results; their permutation and conditional statistics are distinct from the nominal tests calculated for a selected variant.

## Nominal test model

Inference uses the released normalized phenotype, an intercept, alternate-allele dosage and the released **modality-specific covariates**. Both phenotype and genotype are residualized against centered covariates. The two-sided t test has `n - k - 2` degrees of freedom. An eight-entry per-worker cache stores covariate QR decompositions, computed lazily in RAM.

Missing hardcalls are mean-imputed for inference only. Displayed genotypes retain missing values; genotype counts and frequencies use fully called diploid biallelic samples. `maf` is minor allele frequency and `altAlleleFrequency` is alternate allele frequency. Both are null when no calls are available.

Sample identities must match the covariate matrix. Missing covariates, mismatched samples or rank-deficient covariates return `503`. Constant/nonfinite vectors and all-missing genotypes return null statistics with an explicit `associationStatus`. Variants outside the phenotype's 1 Mb cis window, or failing v4's tissue-cohort MAF 0.01 / MAC 2 filters, also return null statistics.

Responses include `pValueSource=computed`, `covariateScope=modality`, association sample count, imputed genotype count and residual degrees of freedom. Float64 CPU arithmetic can differ slightly from tensorQTL's float32 GPU arithmetic and rounded release exports. No request-time permutations or q-values are computed. The precomputed modality-specific gene/group threshold provides significance calibration; grouped modalities use the same gene-group threshold for all that gene's phenotypes. Published cross-modality associations use a distinct covariate model.

## Index preparation and memory

Phenotype SQLite indexes follow the same path and schema across modalities:

```text
data/v4/api/phenotypes/{tissue}/{modality}.{scale}.sqlite3
```

`scale` is `normalized` for all modalities, or `tpm` for expression. Each database stores phenotype IDs, gene IDs, BED coordinates and compressed little-endian float64 sample vectors. Normalized indexes also embed modality-specific mapping thresholds. Sample IDs are stored once per database. Read-only queries retrieve only the requested gene's vectors; they do not load a full tissue/modality matrix.

Build indexes from the existing `ratgtex-server-data` outputs:

```sh
# Normalized phenotypes for all tissues and six modalities
python -m api.prepare_phenotypes --data-dir /path/to/data/v4

# Expression TPM indexes for expression-only views
python -m api.prepare_phenotypes --data-dir /path/to/data/v4 --modalities expression --scales tpm

# Targeted preparation
python -m api.prepare_phenotypes --data-dir /path/to/data/v4 --tissues IL --modalities splicing isoforms

# Exon and independent-association catalog
python -m api.prepare_catalog --data-dir /path/to/data/v4

# Compare computed tests against strong, weak and fixed random released leads
python -m api.validate_associations --data-dir /path/to/data/v4 --report nominal-validation.json
```

Inputs are normalized phenotype BEDs, indexed cohort VCFs, tissue rat IDs, modality-specific covariates, modality-specific `xqtl/top_assoc.*` summaries and the reference/tissue tables. TPM indexes use `expr/expr.tpm.*` when present, with `phenos/*.expression.unnorm.*` as an alternative source. No all-pairs nominal p-value export is required. The API does not read `cis_pvals/*.zip`; archive retention for public downloads is a separate deployment decision.

The exon/association catalog is `data/v4/api/catalog.sqlite3`, prepared from `exon.v4.tsv`, `eqtl/eqtls_indep.v4_rn8.tsv` and `xqtl/xqtls_indep.cross_modality.v4_rn8.tsv`. Its typed tables are indexed by gene, tissue and modality. `/exon`, `/eqtl` and `/xqtl` filter, count and paginate inside SQLite, retrieving only the selected records. Results retain source-file order, and numeric values are stored without additional rounding. Source TSVs remain suitable for downloads.

The catalog builder streams the three inputs, records their sizes/mtimes and schema version, and atomically publishes the complete database. It requires no new scientific exports. Add `api.prepare_catalog` after the three input tables are produced in the HPC release workflow; transfer the database with those matching sources. The three catalog endpoints require a current database and return `503` for missing, invalid or stale catalogs, with no full-table pandas fallback. Source validation is specific to the requested collection. Significant cis associations remain in their per-gene ZIP archive.

Build indexes on the HPC after `ratgtex-server-data` finishes, using the same API code version that will serve them. Keep the builder in this repository and invoke its pinned version from the HPC workflow. Index creation streams one source matrix at a time, checks identifiers and coordinates, and replaces each database atomically. Reserve space for the largest temporary index as well as the final databases. On-server building is also available for repairs.

Transfer indexes and matching source files together with mtimes preserved (`rsync -a` or equivalent). Metadata records the index format, tissue, modality, scale and source sizes/mtimes, including the threshold source. SQLite databases and their explicit little-endian vector encoding are portable between HPC and VPS. Record the serving API commit and validation report with the release. Publish the complete release and restart workers together; per-database replacement is not a transaction across a release's files.

`/genotypePhenotypes` requires current indexes and returns `503` if an index is missing or stale. Expression-only queries can fall back to loading a source matrix. Their shared fallback matrix cache retains at most 128 MiB per worker by default; a larger matrix may be loaded temporarily without retention. This budget does not limit total process memory: reference tables, scientific libraries, response bodies and temporary computation also consume memory.

## Serving and caching

| Environment variable | Default | Meaning |
| --- | --- | --- |
| `RATGTEX_DATA_DIR` | Repository `data/v4` | Release directory |
| `RATGTEX_EXPRESSION_CACHE_MB` | `128` | Retained expression fallback matrix budget per worker, in MiB; `0` disables retention |
| `RATGTEX_RESPONSE_CACHE_MB` | `32` | Serialized response cache budget per worker, in MiB; `0` disables it |
| `RATGTEX_REFERENCE_DATASET` | `ratgtex_v4_round11_2` when available, otherwise the first cohort | Default cohort for variant metadata |

Successful JSON responses are cached for five minutes within a byte budget and a 256-entry limit. Responses support weak ETags, conditional GET/HEAD and negotiated gzip. Errors are not cached. Gene and tissue reference tables load lazily under a lock; exon and independent-association records are queried from the catalog. Each catalog connection uses a 2 MiB SQLite page-cache budget. VCF handles and SQLite connections are request-local and closed after use. ZIP handles used for significant-association lookups are bounded and closed on eviction. Workers have independent caches that assume immutable release files; restart them after data updates. Measure memory before increasing worker count.

## Query conventions

- Gene IDs are rat gene symbols or annotation IDs, matched exactly or unambiguously without regard to case. `/genes?query=Cd300&limit=100&offset=0` searches the gene directory.
- Variant IDs use `chromosome:position`, with 1-based positions. `/variant` accepts `genotypeDatasetId`; genotype/phenotype queries select the cohort from the tissue. Ambiguous positions and multiallelic records cannot identify a supported biallelic variant.
- Use `tissueId` or `tissueSiteDetailId` for tissue filtering, supplying only one. Responses use `tissueSiteDetailId`. List parameters are comma-separated and deduplicated; repeated parameters and empty identifiers return `400`.
- `/eqtl`, `/xqtl`, `/significantCisEqtls` and `/exon` accept `limit` (1–1000) and `offset`. Omitting both returns the full collection. Catalog endpoints apply pagination before records enter Python memory. A significant-association ZIP member is decoded before its pagination.
- `/medianGeneExpression` supports up to 500 genes with `cluster=false`; clustering supports up to 200 genes and uses raw-TPM average-linkage clustering. Unavailable clusters have the value `"Not enough data"`.
- `pValue` is nominal; `pValueBeta` is a beta-approximated permutation p-value; `qValue` is a q-value. `slope` (also exposed as `nes` in expression associations) refers to normalized-phenotype regression, while `log2aFC` is log2 allelic fold change.
- Empty collection lookups return empty arrays. Endpoints requiring one existing resource return `404` when it is absent. Invalid parameters return `400`, unsupported variants `422`, and unavailable release data or indexes `503`. Missing/nonfinite numbers are JSON null.
