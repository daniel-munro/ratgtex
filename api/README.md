# RatGTEx API

The v4 API serves the compiled RatGTEx release files. The existing portal routes and response keys remain supported. `app:api` still works with the current Gunicorn configuration. API schema version 4.1 is an API revision, not a new scientific data release.

## Run and test

From the repository root:

```sh
python -m pip install -r requirements.txt
python -m api.app
python -m unittest discover -s tests -p 'test_*.py' -v
```

Data paths resolve relative to the repository, independently of the working directory. Set `RATGTEX_DATA_DIR` to the absolute path of a different **v4** release directory. This code still uses v4 filenames; changing this setting does not add support for another release layout.

The Flask factory `create_app(config=None, store=None)` supports isolated tests. Importing the application opens no data files or VCF handles. Reference tables load only when requested, under a lock to avoid duplicate cold loads in threaded workers. Gene lookups for exon and association tables use indexes instead of scanning an entire table per request.

## Expression indexes and memory

Previously, startup loaded every tissue's full sample-expression matrix twice, once as `tpm` and again as `iqn`; both loaded the same unnormalized files. There is now one bounded fallback matrix cache, shared by expression and genotype-expression requests. The default retained matrix budget is 128 MiB per worker; an oversized matrix can still be loaded temporarily to serve a cold request but is not retained. The budget is **not** a limit on total process memory: reference tables, ZIP directories, Python/scientific libraries, response bodies, and temporary computation also consume memory.

For regular service, build the optional expression indexes:

```sh
# All tissues, both TPM and inverse-normal transformed expression
python -m api.prepare_expression --data-dir /path/to/data/v4

# Build just normalized expression when TPM indexes already exist
python -m api.prepare_expression --data-dir /path/to/data/v4 --scales normalized

# Rebuild just one tissue and scale
python -m api.prepare_expression --data-dir /path/to/data/v4 --tissues IL --scales tpm
```

This streams each gzip BED into `data/v4/api/expression/{tissue}.{scale}.sqlite3`. Each row contains one gene's zlib-compressed float64 vector; sample IDs are stored once per file. Requests use read-only SQLite queries and decompress only the requested vector. There is no full expression matrix in memory on this path. SQLite connections are closed after use. Building an index uses temporary files and atomic replacement, and duplicate gene/sample IDs or malformed rows abort the build. Index generation adds disk usage (about 1.9 GiB for all current local v4 TPM and normalized indexes); this is a disk-for-memory/latency tradeoff, not an additional required source format.

A source size and modification-time check prevents an index from serving values for a changed source file. Transfer matching BEDs and indexes with their modification times preserved, or rebuild the indexes on the serving host. Stale/missing indexes fall back to the original gzip BEDs. TPM prefers `expr/expr.tpm.*`; the equivalent `phenos/*.expression.unnorm.*` is accepted if the former is absent. `normalized` always uses `phenos/*.expression.norm.*` and never silently falls back to TPM.

The existing `ratgtex-server-data` outputs are sufficient. Indexes can be generated locally, on the web server before deployment, or after compilation on the HPC; no upstream data regeneration is required for the API changes in this revision.

## Serving and cache configuration

| Environment variable | Default | Meaning |
| --- | --- | --- |
| `RATGTEX_DATA_DIR` | Repository `data/v4` | Source release directory |
| `RATGTEX_EXPRESSION_CACHE_MB` | `128` | Maximum retained fallback expression matrix bytes, in MiB, per worker; `0` disables retention |
| `RATGTEX_RESPONSE_CACHE_MB` | `32` | Serialized response cache budget, in MiB, per worker; `0` disables it |
| `RATGTEX_REFERENCE_DATASET` | `ratgtex_v4_round11_2` when present, otherwise first tissue cohort | Default cohort for variant metadata |

Successful JSON responses are cached for five minutes within both a byte budget and a 256-entry cap. Responses too large for the budget are not retained. Responses support weak ETags and conditional GET/HEAD, and gzip when accepted and beneficial. The representation cache stores uncompressed JSON; compression is negotiated per request. Errors are not cached. ZIP directories/handles are cached with an LRU limit of four per worker; archive members are looked up directly, without rebuilding a list of every member on every lookup. Cis p-value lookups scan one gene's member rather than constructing a full DataFrame.

VCFs are opened lazily in request-local contexts and closed after use; no shared mutable iterator or native handle is inherited from an eagerly loaded application. Variant metadata skips sample decoding. Indexed lookups use exact VCF positions, excluding overlapping records that start at another position. Position-only IDs that resolve to multiple records are treated as unavailable, because the identifier cannot disambiguate the alleles.

Caches assume immutable release files. **Publish the matching source files/indexes together and restart all workers after an update.** There is no cross-worker cache invalidation. Avoid multiplying workers before measuring their memory usage. The current Circus command can remain unchanged; this change does not tune worker count, proxy compression, or the existing infinite Gunicorn timeout.

## Query and schema conventions


New integrations and the site use the clearer endpoint names below. The old paths remain supported, with identical parameters, defaults, response bodies and error behavior; OpenAPI marks them deprecated and `/metadata` lists the alias mappings.

| Preferred endpoint | Deprecated alias | Meaning |
| --- | --- | --- |
| `/genotypeExpression` | `/dyneqtl` | Per-sample expression grouped by a variant's genotype |
| `/significantCisEqtls` | `/singleTissueEqtl` | Significant cis-eQTL associations across selected tissues |

The `singleTissueEqtl` response collection key is retained on both association paths to avoid a separate response-schema migration. Other existing field names remain supported. Adding preferred paths is backward compatible within API 4.x; deleting old paths or renaming their required fields would need a breaking API revision. API contract revisions and scientific data releases do not have to happen together. There is no scheduled removal date for these aliases.

The machine-readable contract is [openapi.yaml](../openapi.yaml); the generated public documentation is `pages/about_api.html`. `/api/v4/metadata` describes release/assembly, coordinate conventions, dataset IDs, terminology, endpoints and request limits.

- Gene IDs are the rat gene symbols/annotation IDs in the release, not Ensembl IDs. Matching prefers an exact ID, then an unambiguous case-insensitive match. `/gene` returns `notFound` for unresolved IDs. `/genes?query=Cd300&limit=100&offset=0` provides a searchable, paginated directory; omit `query` to traverse all genes.
- New queries may use `tissueId`; `tissueSiteDetailId` remains an alias. Responses retain `tissueSiteDetailId`. Supplying both aliases, repeating a query parameter, empty list elements, invalid booleans or unknown tissues gives a JSON `400`, rather than a server exception. Lists are comma-separated and deduplicated in first-seen order.
- `/geneExpression` accepts up to 10 genes and optional tissue filtering. TPM is the default. `expressionScale=normalized` selects the quantile-normalized, inverse-normal transformed expression phenotype as released; the API applies no additional covariate adjustment. `includeSampleIds=true` exposes IDs aligned with the data vector. Missing gene/tissue measurements are omitted, not invented or filled with zeros. `datasetId` remains the release identifier; `genotypeDatasetId` identifies the cohort.
- `/expressionSummary` returns sample/nonmissing counts, mean, median, quartiles, range and zero count for the same vectors. It does not send the vectors. Finite observations determine summaries. `/medianGeneExpression` continues to use the precompiled, rounded median matrix and the original raw-TPM average-linkage clustering. Set `cluster=false` for automated extraction without clustering. Clustering is capped at 200 genes; requests without clustering support 500. Empty/singleton/invalid-distance cases use the existing `"Not enough data"` sentinel.
- `/genotypeExpression` returns genotype-stratified expression, **not** a dynamic/temporal eQTL test. It defaults to the released quantile-normalized, inverse-normal transformed phenotype, matching both eQTL plots. Both frontend callers also explicitly request `expressionScale=normalized`. Request `expressionScale=tpm` to obtain unnormalized TPM instead. Ordinary `/geneExpression` requests default to TPM, and the gene-expression plots explicitly select TPM. Genotypes are alternate-allele dosages; partial/missing/non-diploid/non-biallelic calls and expression samples absent from the VCF become `null`. Only fully called diploid samples contribute to frequencies. **`maf` is corrected to minor allele frequency**; `altAlleleFrequency` exposes the quantity that used to be mislabeled `maf`. Both are `null` if there are no called samples. Nominal p-values and thresholds come from compiled mapping results; changing the displayed expression scale does not refit the association.
- `/eqtl` and `/xqtl` retain their original fields and now expose available regression slopes/standard errors, nominal p-values, alternate-allele frequency, rank and position information from the existing source tables. A beta-approximated permutation p-value (`pValueBeta`) is distinct from a nominal p-value (`pValue`) and a q-value. `log2aFC` is distinct from a normalized-phenotype regression slope. Tissue filters are supported; `/xqtl` also accepts `modality`.
- `/significantCisEqtls` contains significant associations, not all tested cis-window variants. It retains `nes` and adds `slope` for the same regression coefficient. Additional runtime p-value rounding was removed; precision already lost upstream cannot be restored here.
- `/eqtl`, `/xqtl`, `/significantCisEqtls` and `/exon` accept optional `limit` (1–1000) and `offset`. Without either, they return the full legacy collection. When either is supplied, `pagination` includes total/limit/offset. These limits bound transfer size; a ZIP member must still be decoded before pagination.
- `/variant` accepts `genotypeDatasetId`, rather than silently assuming that the default cohort contains every release variant. IDs absent/ambiguous/multiallelic in that cohort appear in `notFound`. For genotype-expression requests, the tissue selects the cohort automatically.
- Unknown genes return empty collection lookups for `/gene`, `/geneExpression` and `/medianGeneExpression`; endpoints requiring exactly one existing gene return `404`. Malformed parameters return `400`, unsupported multiallelic genotype-expression requests return `422`, and missing source files return `503`. All numbers are valid JSON: missing/nonfinite numbers are `null`, never `NaN` or `Infinity`.

These are additive schema changes except for the scientific `maf` correction, removal of extra p-value rounding, and intentional validation/error handling changes. The site's existing fields and clustering behavior remain supported. `/geneExpression` defaults to TPM; `/genotypeExpression` now defaults to inverse-normal transformed values, correcting the old copy/paste loading bug and matching the eQTL plot labels. Automated callers that need the former TPM values from `/genotypeExpression` should explicitly set `expressionScale=tpm`.

## Suggested next HPC release changes

These are follow-ups, not prerequisites or changes made to `ratgtex-server-data` in this revision:

1. Add an index-generation step after copying the expression BEDs, or run `api.prepare_expression` after syncing the completed release. Keep it in this repository for now so the index format and reader evolve together.
2. In `scripts/all_gene_esnps.py`, replace the current `float_format="%g"` (six significant digits) with a precision-preserving format such as `"%.17g"`. Likewise review formatting in upstream R exports. Existing archives cannot recover their discarded digits; changing those needs an HPC rebuild.
3. Produce a small release manifest with the assembly/annotation version, pipeline/software versions, normalization method per modality, cohort descriptions, file checksums and exact definitions of each statistic. That would allow richer provenance than the current filename-derived metadata and make future releases easier to validate automatically.
4. If future visualizations need genome-wide or cross-gene association searches, export a queryable index keyed by gene, tissue, modality and genomic position instead of asking the web process to scan bulk ZIP/GZIP files. Existing per-gene archives remain appropriate for current visualizations. Add other phenotype modalities and cross-tissue mash results only with explicit effect-scale and model-provenance fields.
