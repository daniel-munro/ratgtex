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
from api.data_store import DataStore, phenotype_index_path
from api.prepare_phenotypes import build_index
from api.prepare_catalog import build_catalog


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
        self.write("xqtl/top_assoc.expression.v4_rn8.tsv", [{"tissue": "Brain", "gene_id": "GeneA", "pval_nominal_threshold": .001}])
        self.write("exon.v4.tsv", [{"geneId": "GeneA", "chromosome": "chr1", "start": 1, "end": 5}])
        self.write("covar/covar.Brain.expression.v4.tsv", [{"ID": "c1", "s1": 0, "s2": 1, "s3": 0, "s4": -1, "s5": 0}])
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
        for sample in ["s1", "s2", "s3", "s4", "s5"]:
            header.add_sample(sample)
        with pysam.VariantFile(vcf, "w", header=header) as handle:
            # This deletion overlaps chr1:10 but must not hide the exact SNP.
            for position, alleles, calls in [(9, ("AA", "A"), [(0, 0)] * 4),
                    (10, ("A", "G"), [(1, 1), (1, 0), (None, 1), (None, None)]),
                    (20, ("A", "G"), [(None, None)] * 4),
                    (30, ("A", "G", "T"), [(0, 1)] * 4),
                    (50, ("A", "G"), [(0, 0), (0, 1), (1, 1), (0, 0), (0, 1)])]:
                rec = handle.new_record(contig="chr1", start=position - 1, alleles=alleles)
                for sample, call in zip(header.samples, calls):
                    rec.samples[sample]["GT"] = call
                handle.write(rec)
        pysam.tabix_compress(str(vcf), str(vcf) + ".gz", force=True)
        pysam.tabix_index(str(vcf) + ".gz", preset="vcf", force=True)
        build_index(self.root / "phenos/phenos.Brain.expression.norm.v4_rn8.bed.gz",
                    phenotype_index_path(self.root, "Brain", "expression"), "expression",
                    self.root / "xqtl/top_assoc.expression.v4_rn8.tsv", "Brain")
        build_catalog(self.root)
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

    def phenotype(self, query):
        response = self.get(query)
        return {**response, **response["genotypePhenotypes"][0]}

    def test_startup_is_lazy_and_paths_are_independent_of_cwd(self):
        app = create_app({"DATA_DIR": self.root / "nonexistent", "TESTING": True})
        store = app.extensions["ratgtex_store"]
        self.assertNotIn("genes", store.__dict__)
        self.assertNotIn("tissue_info", store.__dict__)
        self.assertEqual(app.test_client().get("/unknown").status_code, 404)
        self.get("metadata")
        self.assertFalse(self.store._expression)
        self.assertNotIn("genes", self.store.__dict__)

    def test_collections_and_association_statistics(self):
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
        association = self.get("significantCisEqtls?geneId=GeneA")["significantCisEqtls"][0]
        self.assertEqual(association["slope"], association["nes"])
        self.assertEqual(association["pValue"], 1.23456789e-8)
        self.assertEqual(self.get("significantCisEqtls?geneId=Empty")["significantCisEqtls"], [])

    def test_canonical_endpoints(self):
        meta = self.get("metadata")
        self.assertIn("/api/v4/genotypePhenotypes", meta["endpoints"])
        self.assertIn("/api/v4/significantCisEqtls", meta["endpoints"])
        self.assertNotIn("endpointAliases", meta)
        for path in ["dyneqtl", "genotypeExpression", "singleTissueEqtl"]:
            self.assertEqual(self.client.get("/api/v4/" + path).status_code, 404)

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

    def test_catalog_filters_counts_order_and_numeric_values(self):
        import sqlite3
        rows = [
            {"tissue": tissue, "gene_id": gene, "variant_id": variant, "modality": modality,
             "phenotype_id": gene + "__p", "pval_nominal": pvalue, "slope": slope,
             "rank": rank, "chrom": "chr1", "pos": pos}
            for tissue, gene, variant, modality, pvalue, slope, rank, pos in [
                ("Brain", "GeneA", "chr1:30", "splicing", 1.234567890123456e-30, -.25, 2, 30),
                ("Brain", "GeneB", "chr1:20", "expression", .2, .4, 1, 20),
                ("Liver", "GeneA", "chr1:10", "splicing", .1, .3, 1, 10),
                ("Brain", "GeneA", "chr1:40", "expression", .3, None, 1, 40),
                ("Brain", "GeneA", "chr1:50", "splicing", .4, .7, 3, 50),
            ]]
        self.write("xqtl/xqtls_indep.cross_modality.v4_rn8.tsv", rows)
        self.write("eqtl/eqtls_indep.v4_rn8.tsv", [{k: v for k, v in row.items() if k not in {"modality", "phenotype_id"}} for row in rows])
        self.write("exon.v4.tsv", [{"geneId": gene, "chromosome": "chr1", "start": start, "end": start + 1}
                                    for gene, start in [("GeneA", 30), ("GeneB", 20), ("GeneA", 10)]])
        build_catalog(self.root)
        self.store.genes
        self.store.tissues
        traces = []
        connect = sqlite3.connect
        def traced(*args, **kwargs):
            conn = connect(*args, **kwargs)
            conn.set_trace_callback(traces.append)
            return conn
        with patch("api.data_store.pd.read_csv", side_effect=AssertionError("Catalog must avoid loading source tables")), patch("api.catalog.sqlite3.connect", side_effect=traced):
            result = self.get("xqtl?geneId=GeneA&tissueId=Brain&modality=splicing&limit=1&offset=1")
            self.assertEqual(result["pagination"], {"total": 2, "limit": 1, "offset": 1})
            self.assertEqual(result["xqtl"][0]["variantId"], "chr1:50")
            self.assertEqual(self.get("eqtl?geneId=GeneA&tissueId=Brain&limit=1")["pagination"]["total"], 3)
            self.assertEqual(self.get("exon?geneId=GeneA&limit=1&offset=1")["exon"][0]["start"], 10)
            unpaged = self.get("xqtl?geneId=GeneA")
            self.assertEqual([r["pos"] for r in unpaged["xqtl"]], [30, 10, 40, 50])
            self.assertEqual(unpaged["xqtl"][0]["pValue"], 1.234567890123456e-30)
            self.assertIsNone(unpaged["xqtl"][2]["slope"])
            self.assertNotIn("pagination", unpaged)
            self.assertEqual(self.get("xqtl?geneId=Empty&limit=1")["pagination"]["total"], 0)
            self.assertEqual(self.get("xqtl?geneId=GeneA&modality=isoforms")["xqtl"], [])
            self.assertEqual(self.get("xqtl?geneId=GeneA&limit=1&offset=100")["xqtl"], [])
        self.assertTrue(any("LIMIT 1 OFFSET 1" in sql for sql in traces))
        self.assertTrue(any("count(*)" in sql for sql in traces))
        for name in ["exons", "eqtls", "xqtls"]:
            self.assertNotIn(name, self.store.__dict__)
        self.assertEqual(self.store.catalog_records("exon", "GeneA' OR 1=1 --")[0], [])

    def test_catalog_missing_stale_and_invalid_are_service_errors(self):
        import os
        import sqlite3
        catalog = self.root / "api/catalog.sqlite3"
        catalog.unlink()
        with self.assertLogs(self.app.logger, level="ERROR"):
            self.assertEqual(self.client.get("/api/v4/exon?geneId=GeneA").status_code, 503)
        build_catalog(self.root)
        source = self.root / "exon.v4.tsv"
        os.utime(source, ns=(source.stat().st_atime_ns, source.stat().st_mtime_ns + 1000000000))
        with self.assertLogs(self.app.logger, level="ERROR"):
            self.assertEqual(self.client.get("/api/v4/exon?geneId=GeneA").status_code, 503)
        # Unrelated catalog collections remain usable.
        self.assertEqual(len(self.get("eqtl?geneId=GeneA")["eqtl"]), 1)
        build_catalog(self.root)
        with sqlite3.connect(catalog) as conn:
            conn.execute("UPDATE metadata SET value = 'unknown-format' WHERE key = 'format'")
        with self.assertLogs(self.app.logger, level="ERROR"):
            self.assertEqual(self.client.get("/api/v4/eqtl?geneId=GeneA").status_code, 503)

    def test_failed_catalog_preparation_preserves_published_database(self):
        catalog = self.root / "api/catalog.sqlite3"
        previous = catalog.read_bytes()
        with (self.root / "exon.v4.tsv").open("a") as handle:
            handle.write("bad\trow\n")
        with self.assertRaises(ValueError):
            build_catalog(self.root)
        self.assertEqual(catalog.read_bytes(), previous)
        self.assertEqual(list(catalog.parent.glob("tmp*.sqlite3")), [])

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
        result = self.phenotype("genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:10&includeSampleIds=true")
        self.assertEqual(result["genotypes"], [2, 1, None, None, None])
        self.assertEqual(result["altAlleleFrequency"], .75)
        self.assertEqual(result["maf"], .25)
        self.assertEqual(result["calledSampleCount"], 2)
        self.assertEqual(result["missingGenotypeCount"], 3)
        self.assertIsNone(result["pValue"])
        self.assertEqual(result["associationStatus"], "variant_below_mapping_filter")
        self.assertEqual(result["pValueSource"], "computed")
        self.assertEqual(len(result["data"]), len(result["sampleIds"]))
        missing = self.phenotype("genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:20")
        self.assertIsNone(missing["maf"])
        self.assertIsNone(missing["altAlleleFrequency"])
        self.assertIsNone(missing["pValue"])
        multi = self.client.get("/api/v4/genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:30")
        self.assertEqual(multi.status_code, 422)
        self.assertEqual(self.client.get("/api/v4/genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:999").status_code, 404)
        self.assertIsNone(genotype((0, None)))
        self.assertIsNone(genotype((1,)))
        self.assertIsNone(genotype((0, 2)))

    def test_genotype_phenotypes_use_normalized_values(self):
        result = self.get("genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:50")
        self.assertEqual(result["genotypePhenotypes"][0]["data"], [-2, -1, 0, 1, 2])
        self.assertEqual(result["associationPhenotypeScale"], "normalized")
        self.assertEqual(result["modality"], "expression")
        self.assertEqual(result["unit"], "normalized phenotype")

    def test_exact_variant_lookup_and_dataset_selection(self):
        result = self.get("variant?variantId=chr1:10,chr1:999,chr1:0,invalid,chr1:30&genotypeDatasetId=cohort")
        self.assertEqual(len(result["variant"]), 1)
        self.assertEqual(result["variant"][0]["ref"], "A")
        self.assertEqual(result["notFound"], ["chr1:999", "chr1:0", "invalid", "chr1:30"])

    def test_nominal_expression_uses_normalized_values_and_no_pvalue_archive(self):
        (self.root / "cis_pvals/Brain.v4.zip").unlink()
        normalized = self.phenotype("genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:50")
        self.assertEqual(normalized["associationStatus"], "ok")
        self.assertGreater(normalized["pValue"], 0)
        self.assertEqual(normalized["associationSampleCount"], 5)
        self.assertEqual(normalized["degreesOfFreedom"], 2)

    def test_multimodal_index_and_association_endpoint(self):
        from api.data_store import phenotype_index_path
        samples = [f"s{i}" for i in range(1, 6)]
        source = self.root / "phenos/phenos.Brain.splicing.norm.v4_rn8.bed.gz"
        self.write(str(source.relative_to(self.root)), [
            {"#chr": "chr1", "start": 0, "end": 1, "phenotype_id": name, **dict(zip(samples, data))}
            for name, data in [("GeneA__intron1", [-2, -1, 0, 1, 2]),
                               ("GeneA__intron2", [1, 1, 1, 1, 1]),
                               ("GeneB__intron1", [1, 2, 3, 4, 5])]])
        thresholds = self.root / "xqtl/top_assoc.splicing.v4_rn8.tsv"
        self.write(str(thresholds.relative_to(self.root)), [{"tissue": "Brain", "gene_id": "GeneA", "pval_nominal_threshold": .002}])
        # A different covariate from expression proves modality-specific selection.
        self.write("covar/covar.Brain.splicing.v4.tsv", [{"ID": "c2", **dict(zip(samples, [1, 0, -1, 0, 1]))}])
        index = phenotype_index_path(self.root, "Brain", "splicing")
        build_index(source, index, "splicing", thresholds, "Brain")
        query = "genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:50&modality=splicing&includeSampleIds=true"
        # Prime only the small reference/covariate tables, before the guard.
        self.store.genes
        self.store.tissues
        self.store.datasets
        self.store.covariate_model("Brain", "splicing")
        with patch("api.data_store.pd.read_csv", side_effect=AssertionError("No full matrices should load")):
            result = self.get(query)
        self.assertEqual(result["sampleIds"], samples)
        self.assertEqual(result["pValueThreshold"], .002)
        self.assertEqual([r["phenotypeId"] for r in result["genotypePhenotypes"]], ["GeneA__intron1", "GeneA__intron2"])
        self.assertEqual(result["genotypePhenotypes"][0]["associationStatus"], "ok")
        self.assertEqual(result["genotypePhenotypes"][1]["associationStatus"], "no_residual_phenotype_variation")
        self.assertEqual(len(self.get(query + "&phenotypeId=GeneA__intron1")["genotypePhenotypes"]), 1)
        self.assertEqual(self.get(query + "&phenotypeId=GeneB__intron1")["genotypePhenotypes"], [])
        self.assertEqual(self.client.get("/api/v4/" + query.replace("modality=splicing", "modality=invalid")).status_code, 400)
        import os
        os.utime(thresholds, ns=(thresholds.stat().st_atime_ns, thresholds.stat().st_mtime_ns + 1000000000))
        with self.assertLogs(self.app.logger, level="ERROR"):
            self.assertEqual(self.client.get("/api/v4/" + query).status_code, 503)

    def test_phenotype_inventory_does_not_read_sample_vectors(self):
        with patch("api.data_store.zlib.decompress", side_effect=AssertionError("Inventory must not decompress values")):
            result = self.get("phenotypes?geneId=GeneA&tissueId=Brain")
        self.assertEqual(len(result["phenotypes"]), 1)
        phenotype = result["phenotypes"][0]
        self.assertEqual(phenotype["phenotypeId"], "GeneA")
        self.assertEqual(phenotype["sampleCount"], 5)
        self.assertEqual(phenotype["pValueThreshold"], .001)
        self.assertNotIn("data", phenotype)
        self.assertEqual(result["datasets"][0], {"modality": "expression", "tissueSiteDetailId": "Brain", "status": "ready", "phenotypeCount": 1})
        self.assertTrue(all(d["status"] == "unavailable" for d in result["datasets"][1:]))
        partial = self.client.get("/api/v4/phenotypes?geneId=GeneA&tissueId=Brain")
        self.assertIn("no-store", partial.headers["Cache-Control"])
        available = self.client.get("/api/v4/phenotypes?geneId=GeneA&tissueId=Brain&modality=expression")
        self.assertNotIn("no-store", available.headers["Cache-Control"])
        empty = self.get("phenotypes?geneId=Empty&tissueId=Brain&modality=expression")
        self.assertEqual(empty["phenotypes"], [])
        self.assertEqual(empty["datasets"][0]["status"], "ready")
        self.assertEqual(empty["datasets"][0]["phenotypeCount"], 0)
        self.assertEqual(self.client.get("/api/v4/phenotypes?geneId=GeneA&modality=bad").status_code, 400)

    def test_missing_and_stale_indexes_fail_without_matrix_fallback(self):
        query = "/api/v4/genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:50"
        index = phenotype_index_path(self.root, "Brain", "expression")
        index.unlink()
        with self.assertLogs(self.app.logger, level="ERROR"):
            self.assertEqual(self.client.get(query).status_code, 503)

    def test_association_sample_alignment_and_missing_covariates(self):
        query = "genotypePhenotypes?geneId=GeneA&tissueId=Brain&variantId=chr1:50"
        expected = self.phenotype(query)["pValue"]
        # Different order must preserve sample-identity alignment.
        self.store.close()
        self.write("covar/covar.Brain.expression.v4.tsv", [{"ID": "c1", "s5": 0, "s4": -1, "s3": 0, "s2": 1, "s1": 0}])
        self.assertAlmostEqual(self.phenotype(query)["pValue"], expected)
        self.store.close()
        self.write("covar/covar.Brain.expression.v4.tsv", [{"ID": "c1", "s5": 0, "s4": -1, "s3": 0, "s2": 1, "wrong": 0}])
        with self.assertLogs(self.app.logger, level="ERROR"):
            response = self.client.get("/api/v4/" + query)
        self.assertEqual(response.status_code, 503)
        self.assertIn("samples do not match", response.json["error"])
        self.store.close()
        (self.root / "covar/covar.Brain.expression.v4.tsv").unlink()
        with self.assertLogs(self.app.logger, level="ERROR"):
            self.assertEqual(self.client.get("/api/v4/" + query).status_code, 503)

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
        index = phenotype_index_path(self.root, "Brain", "expression", "tpm")
        self.assertEqual(build_index(source, index, "expression", tissue="Brain", scale="tpm"), 1)
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
            build_index(bad, index, "expression", tissue="Brain", scale="tpm")
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
            tiny.significant_associations("GeneA")
            archive = next(iter(tiny._archives.values()))
            self.assertEqual(len(tiny._archives), 1)
            tiny.close()
            self.assertIsNone(archive.fp)
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
