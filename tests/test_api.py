"""API integration tests with tiny release files, indexed VCFs, and ZIP archives."""
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
import pandas as pd
import pysam

from api.app import create_app, genotype
from api.data_store import DataStore
from api.prepare_expression import build_index


class APITest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.write("tissue_info.v4.tsv", [{"tissueSiteDetailId": "Brain", "dataset": "cohort", "tissueSiteDetail": "Brain", "colorHex": "123456"}])
        genes = [{"geneId": gene, "chromosome": "chr1", "start": 1, "end": 100,
                  "description": "A test gene", "hasEqtl": True, "hasTwas": False,
                  "expr_Brain": True, "testedEqtl_Brain": True, "eqtl_Brain": True}
                 for gene in ["GeneA", "GeneB", "Empty"]]
        self.write("gene.v4.tsv", genes)
        self.write("expr/medianGeneExpression.v4.tsv.gz", [{"geneId": "GeneA", "Brain": 3}, {"geneId": "GeneB", "Brain": 2}])
        self.write("expr/topExpressedGene.v4.tsv", [{"geneId": "GeneA", "tissueSiteDetailId": "Brain", "median": 3, "mtGene": False}, {"geneId": "GeneB", "tissueSiteDetailId": "Brain", "median": 2, "mtGene": True}])
        for scale in ["tpm", "normalized"]:
            name = "expr/expr.tpm.Brain.v4_rn8.bed.gz" if scale == "tpm" else "phenos/phenos.Brain.expression.norm.v4_rn8.bed.gz"
            values = [1, 2, 3, 4, 5] if scale == "tpm" else [-2, -1, 0, 1, 2]
            self.write(name, [{"#chr": "chr1", "start": 0, "end": 1, "phenotype_id": "GeneA", **dict(zip(["s1", "s2", "s3", "s4", "s5"], values))}])
        assoc = [{"tissue": "Brain", "gene_id": "GeneA", "variant_id": "chr1:10", "ref": "A", "alt": "G", "pval_beta": .01, "log2_aFC": .5, "rank": 1, "af": .75, "pval_nominal": 1e-8, "slope": .2, "slope_se": .05}]
        self.write("eqtl/eqtls_indep.v4_rn8.tsv", assoc)
        self.write("xqtl/xqtls_indep.cross_modality.v4_rn8.tsv", [{**assoc[0], "modality": "expression", "phenotype_id": "GeneA"}])
        self.write("eqtl/top_assoc.v4_rn8.tsv", [{"tissue": "Brain", "gene_id": "GeneA", "pval_nominal_threshold": .001}])
        self.write("exon.v4.tsv", [{"geneId": "GeneA", "chromosome": "chr1", "start": 1, "end": 5}])
        with zipfile.ZipFile(self.root / "eqtl/all_gene_esnps.v4.zip", "w") as archive:
            archive.writestr("all_gene_esnps.v4/GeneA.txt", "variantId\ttissueSiteDetailId\tpValue\tnes\nchr1:10\tBrain\t1.23456789e-8\t0.2\n")
        (self.root / "cis_pvals").mkdir()
        with zipfile.ZipFile(self.root / "cis_pvals/Brain.v4.zip", "w") as archive:
            archive.writestr("Brain.v4/GeneA.txt", "variant_id\tpval_nominal\nchr1:10\t1e-8\n")
        (self.root / "geno").mkdir()
        vcf = self.root / "geno/cohort.rn8.vcf"
        header = pysam.VariantHeader()
        header.contigs.add("chr1", length=1000)
        header.formats.add("GT", 1, "String", "Genotype")
        for sample in ["s1", "s2", "s3", "s4"]:
            header.add_sample(sample)
        with pysam.VariantFile(vcf, "w", header=header) as handle:
            # This deletion overlaps chr1:10 but must not hide the exact SNP.
            for position, alleles, calls in [(9, ("AA", "A"), [(0, 0)] * 4),
                    (10, ("A", "G"), [(1, 1), (1, 0), (None, 1), (None, None)]),
                    (20, ("A", "G"), [(None, None)] * 4),
                    (30, ("A", "G", "T"), [(0, 1)] * 4)]:
                rec = handle.new_record(contig="chr1", start=position - 1, alleles=alleles)
                for sample, call in zip(header.samples, calls):
                    rec.samples[sample]["GT"] = call
                handle.write(rec)
        pysam.tabix_compress(str(vcf), str(vcf) + ".gz", force=True)
        pysam.tabix_index(str(vcf) + ".gz", preset="vcf", force=True)
        self.store = DataStore(self.root)
        self.app = create_app({"TESTING": True, "RESPONSE_CACHE_BYTES": 0}, self.store)
        self.client = self.app.test_client()

    def tearDown(self):
        self.store.close()
        self.temp.cleanup()

    def write(self, name, rows):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    def get(self, path):
        result = self.client.get("/api/v4/" + path)
        self.assertEqual(result.status_code, 200, result.json)
        json.loads(result.data, parse_constant=lambda value: self.fail("Nonstandard JSON: " + value))
        return result.json

    def test_startup_is_lazy_and_paths_are_independent_of_cwd(self):
        app = create_app({"DATA_DIR": self.root / "nonexistent", "TESTING": True})
        store = app.extensions["ratgtex_store"]
        self.assertNotIn("genes", store.__dict__)
        self.assertNotIn("tissue_info", store.__dict__)
        self.assertEqual(app.test_client().get("/unknown").status_code, 404)
        self.get("metadata")
        self.assertFalse(self.store._expression)
        self.assertNotIn("genes", self.store.__dict__)

    def test_all_legacy_collections_and_additive_statistics(self):
        self.assertEqual(self.get("tissueInfo")["tissueInfo"][0]["tissueSiteDetailId"], "Brain")
        gene = self.get("gene?geneId=GENEA,missing,GeneA")
        self.assertEqual([row["geneId"] for row in gene["gene"]], ["GeneA"])
        self.assertEqual(gene["notFound"], ["missing"])
        self.assertIs(gene["gene"][0]["statusInTissue"][0]["expressed"], True)
        self.assertNotIn("expr_Brain", gene["gene"][0])
        self.assertEqual(self.get("exon?geneId=GeneA")["exon"][0]["start"], 1)
        eqtl = self.get("eqtl?geneId=GeneA&tissueId=Brain")["eqtl"][0]
        self.assertEqual(eqtl["log2aFC"], .5)
        self.assertEqual(eqtl["slopeSE"], .05)
        self.assertEqual(self.get("xqtl?geneId=GeneA&modality=expression")["xqtl"][0]["phenotypeId"], "GeneA")
        association = self.get("singleTissueEqtl?geneId=GeneA")["singleTissueEqtl"][0]
        self.assertEqual(association["slope"], association["nes"])
        self.assertEqual(association["pValue"], 1.23456789e-8)
        self.assertEqual(self.get("singleTissueEqtl?geneId=Empty")["singleTissueEqtl"], [])

    def test_clear_endpoint_names_preserve_legacy_aliases(self):
        for preferred, legacy, query in [
            ("genotypeExpression", "dyneqtl", "geneId=GeneA&tissueId=Brain&variantId=chr1:10"),
            ("significantCisEqtls", "singleTissueEqtl", "geneId=GeneA&tissueId=Brain&limit=1"),
        ]:
            with self.subTest(endpoint=preferred):
                self.assertEqual(self.get(preferred + "?" + query), self.get(legacy + "?" + query))
                self.assertEqual(self.client.get("/api/v4/" + preferred).status_code, 400)
        meta = self.get("metadata")
        self.assertEqual(meta["endpointAliases"]["/api/v4/dyneqtl"], "/api/v4/genotypeExpression")
        self.assertIn("/api/v4/significantCisEqtls", meta["endpoints"])

    def test_parameter_errors_are_json_and_not_cached(self):
        cases = ["gene", "gene?geneId=", "gene?geneId=GeneA,,GeneB", "gene?geneId=GeneA&geneId=GeneB",
                 "geneExpression?geneId=GeneA&tissueId=invalid", "geneExpression?geneId=GeneA&expressionScale=log",
                 "geneExpression?geneId=GeneA&includeSampleIds=1", "topExpressedGene?tissueId=Brain&filterMtGene=[]",
                 "eqtl?geneId=GeneA,GeneB", "eqtl?geneId=GeneA&limit=0", "eqtl?geneId=GeneA&limit=no",
                 "geneExpression?geneId=GeneA&tissueId=Brain&tissueSiteDetailId=Brain", "variant?variantId=chr1:10&genotypeDatasetId=unknown"]
        for path in cases:
            with self.subTest(path=path):
                response = self.client.get("/api/v4/" + path)
                self.assertEqual(response.status_code, 400, response.json)
                self.assertIn("error", response.json)
                self.assertIn("no-store", response.headers["Cache-Control"])

    def test_expression_scales_samples_and_missing_genes(self):
        row = self.get("geneExpression?geneId=GENEA&tissueId=Brain&includeSampleIds=true")["geneExpression"][0]
        self.assertEqual(row["data"], [1, 2, 3, 4, 5])
        self.assertEqual(row["sampleIds"], ["s1", "s2", "s3", "s4", "s5"])
        norm = self.get("geneExpression?geneId=GeneA&expressionScale=normalized")["geneExpression"][0]
        self.assertEqual(norm["data"], [-2, -1, 0, 1, 2])
        self.assertNotIn("sampleIds", norm)
        self.assertEqual(self.get("geneExpression?geneId=Empty")["geneExpression"], [])
        summary = self.get("expressionSummary?geneId=GeneA")["expressionSummary"][0]
        self.assertEqual([summary[key] for key in ["min", "q1", "median", "q3", "max", "mean"]], [1, 2, 3, 4, 5, 3])

    def test_genotypes_frequency_and_missing_values(self):
        result = self.get("dyneqtl?geneId=GeneA&tissueId=Brain&variantId=chr1:10&includeSampleIds=true")
        self.assertEqual(result["genotypes"], [2, 1, None, None, None])
        self.assertEqual(result["altAlleleFrequency"], .75)
        self.assertEqual(result["maf"], .25)
        self.assertEqual(result["calledSampleCount"], 2)
        self.assertEqual(result["missingGenotypeCount"], 3)
        self.assertEqual(result["pValue"], 1e-8)
        self.assertEqual(len(result["data"]), len(result["sampleIds"]))
        missing = self.get("dyneqtl?geneId=GeneA&tissueId=Brain&variantId=chr1:20")
        self.assertIsNone(missing["maf"])
        self.assertIsNone(missing["altAlleleFrequency"])
        self.assertIsNone(missing["pValue"])
        multi = self.client.get("/api/v4/dyneqtl?geneId=GeneA&tissueId=Brain&variantId=chr1:30")
        self.assertEqual(multi.status_code, 422)
        self.assertEqual(self.client.get("/api/v4/dyneqtl?geneId=GeneA&tissueId=Brain&variantId=chr1:999").status_code, 404)
        self.assertIsNone(genotype((0, None)))
        self.assertIsNone(genotype((1,)))
        self.assertIsNone(genotype((0, 2)))

    def test_genotype_expression_defaults_to_inverse_normal_values_and_can_select_tpm(self):
        query = "dyneqtl?geneId=GeneA&tissueId=Brain&variantId=chr1:10"
        for indexed in [False, True]:
            with self.subTest(indexed=indexed):
                if indexed:
                    for scale, source in [("tpm", "expr/expr.tpm.Brain.v4_rn8.bed.gz"),
                                          ("normalized", "phenos/phenos.Brain.expression.norm.v4_rn8.bed.gz")]:
                        build_index(self.root / source, self.root / f"api/expression/Brain.{scale}.sqlite3")
                    self.store.close()
                default = self.get(query)
                normalized = self.get(query + "&expressionScale=normalized")
                tpm = self.get(query + "&expressionScale=tpm")
                self.assertEqual(default["data"], [-2, -1, 0, 1, 2])
                self.assertEqual(default["data"], normalized["data"])
                self.assertEqual(default["expressionScale"], "normalized")
                self.assertEqual(tpm["data"], [1, 2, 3, 4, 5])
                self.assertEqual(tpm["unit"], "TPM")
                self.assertEqual(default["genotypes"], tpm["genotypes"])
                self.assertEqual(self.get("geneExpression?geneId=GeneA")["geneExpression"][0]["data"], tpm["data"])

    def test_exact_variant_lookup_and_dataset_selection(self):
        result = self.get("variant?variantId=chr1:10,chr1:999,chr1:0,invalid,chr1:30&genotypeDatasetId=cohort")
        self.assertEqual(len(result["variant"]), 1)
        self.assertEqual(result["variant"][0]["ref"], "A")
        self.assertEqual(result["notFound"], ["chr1:999", "chr1:0", "invalid", "chr1:30"])

    def test_pagination_gene_search_and_tissue_filters(self):
        result = self.get("genes?limit=2&offset=1")
        self.assertEqual(result["pagination"], {"total": 3, "limit": 2, "offset": 1})
        self.assertEqual(result["gene"][0]["geneId"], "GeneB")
        self.assertEqual(self.get("genes?query=geneb")["pagination"]["total"], 1)
        self.assertEqual(self.get("eqtl?geneId=GeneA&limit=1&offset=1")["eqtl"], [])
        self.assertEqual(len(self.get("topExpressedGene?tissueId=Brain&filterMtGene=true")["topExpressedGene"]), 1)
        self.assertEqual(self.client.get("/api/v4/xqtl?geneId=GeneA&modality=invalid").status_code, 400)

    def test_empty_and_singleton_clustering(self):
        result = self.get("medianGeneExpression?geneId=missing")
        self.assertEqual(result["medianGeneExpression"], [])
        self.assertEqual(result["clusters"], {"gene": "Not enough data", "tissue": "Not enough data"})
        single = self.get("medianGeneExpression?geneId=GeneA")
        self.assertEqual(single["clusters"]["gene"], "Not enough data")
        multiple = self.get("medianGeneExpression?geneId=GeneA,GeneB")
        self.assertTrue(multiple["clusters"]["gene"].endswith(";"))
        disabled = self.get("medianGeneExpression?geneId=GeneA,GeneB&cluster=false")
        self.assertEqual(disabled["clusters"]["gene"], "Not enough data")

    def test_index_matches_fallback_and_stale_index_is_not_used(self):
        fallback = self.store.expression("Brain", "GeneA")
        source = self.root / "expr/expr.tpm.Brain.v4_rn8.bed.gz"
        index = self.root / "api/expression/Brain.tpm.sqlite3"
        self.assertEqual(build_index(source, index), 1)
        self.store.close()
        with patch("api.data_store.pd.read_csv", side_effect=AssertionError("Index should avoid reading a matrix")):
            indexed = self.store.expression("Brain", "GeneA")
        pd.testing.assert_series_equal(indexed, fallback)
        self.assertIsNone(self.store.expression("Brain", "Empty"))
        previous = index.read_bytes()
        bad = self.root / "bad.gz"
        with gzip.open(bad, "wt") as handle:
            handle.write("bad\n")
        with self.assertRaises(ValueError):
            build_index(bad, index)
        self.assertEqual(index.read_bytes(), previous)
        # An updated source invalidates its index and safely uses the fallback.
        import os
        os.utime(source, ns=(source.stat().st_atime_ns, source.stat().st_mtime_ns + 1000000000))
        pd.testing.assert_series_equal(self.store.expression("Brain", "GeneA"), fallback)
        self.assertGreater(self.store._expression_bytes, 0)

    def test_threaded_cold_metadata_loads_only_once(self):
        from concurrent.futures import ThreadPoolExecutor
        read = pd.read_csv
        with patch("api.data_store.pd.read_csv", wraps=read) as mocked:
            with ThreadPoolExecutor(max_workers=4) as pool:
                frames = list(pool.map(lambda _: self.store.genes, range(8)))
        self.assertEqual(mocked.call_count, 1)
        self.assertTrue(all(frame is frames[0] for frame in frames))

    def test_missing_release_data_is_a_service_error(self):
        source = self.root / "phenos/phenos.Brain.expression.norm.v4_rn8.bed.gz"
        source.unlink()
        with self.assertLogs(self.app.logger, level="ERROR"):
            response = self.client.get("/api/v4/geneExpression?geneId=GeneA&expressionScale=normalized")
        self.assertEqual(response.status_code, 503)
        self.assertIn("error", response.json)

    def test_memory_cache_budget_and_zip_handle_bound(self):
        tiny = DataStore(self.root, expression_cache_bytes=1, archive_cache_size=1)
        try:
            self.assertEqual(len(tiny.expression("Brain", "GeneA")), 5)
            self.assertFalse(tiny._expression)
            tiny.cis_pvalue("Brain", "GeneA", "chr1:10")
            old = next(iter(tiny._archives.values()))
            tiny.significant_associations("GeneA")
            self.assertEqual(len(tiny._archives), 1)
            self.assertIsNone(old.fp)
        finally:
            tiny.close()

    def test_http_cache_gzip_conditional_requests_and_strict_json(self):
        self.app = create_app({"TESTING": True, "RESPONSE_CACHE_BYTES": 100000}, self.store)
        client = self.app.test_client()
        first = client.get("/api/v4/gene?geneId=GeneA", headers={"Accept-Encoding": "gzip"})
        compressed = client.get("/api/v4/metadata", headers={"Accept-Encoding": "gzip"})
        self.assertEqual(compressed.headers["Content-Encoding"], "gzip")
        self.assertIn("genomeAssembly", json.loads(gzip.decompress(compressed.data)))
        identity = client.get("/api/v4/metadata", headers={"Accept-Encoding": "gzip;q=0"})
        self.assertNotIn("Content-Encoding", identity.headers)
        self.assertEqual(identity.headers["ETag"], compressed.headers["ETag"])
        with patch.object(self.store, "resolve_gene", side_effect=AssertionError("Response must be cached")):
            cached = client.get("/api/v4/gene?geneId=GeneA")
            self.assertEqual(cached.status_code, 200)
            self.assertEqual(cached.headers["ETag"], first.headers["ETag"])
            conditional = client.get("/api/v4/gene?geneId=GeneA", headers={"If-None-Match": first.headers["ETag"]})
            self.assertEqual(conditional.status_code, 304)
            self.assertEqual(conditional.data, b"")
        with self.app.app_context():
            self.assertEqual(json.loads(self.app.json.dumps({"value": np.nan, "bool": np.bool_(True)})), {"value": None, "bool": True})


if __name__ == "__main__":
    unittest.main()
