"""Offline validation against modality-specific released lead associations.

Uses streamed source BEDs, so no new indexes or p-value archives are required.
This is a release check, not a benchmark of production endpoint latency.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import pysam

from .app import genotype, variant_record
from .data_store import DataStore, MODALITIES


def validate(root, tissues, modalities, count, tolerance):
    store = DataStore(root)
    reports = []
    for modality in modalities:
        summary = pd.read_csv(store.root / "xqtl" / f"top_assoc.{modality}.v4_rn8.tsv", sep="\t",
                              usecols=["tissue", "gene_id", "phenotype_id", "variant_id", "pval_nominal"])
        for tissue in tissues:
            candidates = summary.loc[(summary.tissue == tissue) & (summary.pval_nominal > 0)]
            # Include strongest and weakest results as well as a fixed random sample.
            chosen = pd.concat([candidates.nsmallest(2, "pval_nominal"),
                                candidates.nlargest(2, "pval_nominal"),
                                candidates.sample(n=min(count, len(candidates)), random_state=7301)])
            chosen = chosen.drop_duplicates(["phenotype_id", "variant_id"]).head(count)
            wanted = set(chosen.phenotype_id)
            vectors = {}
            source = store.root / "phenos" / f"phenos.{tissue}.{modality}.norm.v4_rn8.bed.gz"
            with gzip.open(source, "rt", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                samples = next(reader)[4:]
                for row in reader:
                    if row[3] in wanted:
                        vectors[row[3]] = np.asarray(row[4:], dtype=float)
                        if len(vectors) == len(wanted):
                            break
            differences, durations, failures = [], [], []
            with pysam.VariantFile(store.root / "geno" / f"{store.datasets[tissue]}.rn8.vcf.gz") as vcf:
                for row in chosen.itertuples():
                    rec = variant_record(row.variant_id, vcf)
                    if rec is None or row.phenotype_id not in vectors:
                        failures.append({"phenotypeId": row.phenotype_id, "reason": "missing_input"})
                        continue
                    genotypes = pd.Series({s: genotype(rec.samples[s]["GT"]) for s in rec.samples}, dtype=object)
                    frame = pd.DataFrame([vectors[row.phenotype_id]], columns=samples)
                    start = perf_counter()
                    actual = store.nominal_associations(tissue, modality, frame, genotypes)[0]
                    durations.append(1000 * (perf_counter() - start))
                    p = actual["pValue"]
                    difference = abs(np.log10(p) - np.log10(row.pval_nominal)) if p and p > 0 else float("inf")
                    if np.isfinite(difference):
                        differences.append(float(difference))
                    if difference > tolerance:
                        failures.append({"phenotypeId": row.phenotype_id, "variantId": row.variant_id,
                                         "expected": row.pval_nominal, "computed": p,
                                         "status": actual["associationStatus"]})
            report = {"tissue": tissue, "modality": modality, "tested": len(chosen),
                      "maxAbsLog10PDelta": max(differences, default=None),
                      "medianComputeMs": float(np.median(durations)) if durations else None,
                      "failures": failures}
            reports.append(report)
            print(f"{tissue} {modality}: {len(chosen)} tested, {len(failures)} failed, max log10 delta {report['maxAbsLog10PDelta']}", flush=True)
    store.close()
    return {"dataRelease": "v4", "covariateScope": "modality", "log10Tolerance": tolerance,
            "tested": sum(r["tested"] for r in reports),
            "failures": sum(len(r["failures"]) for r in reports), "results": reports}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data/v4")
    parser.add_argument("--tissues", nargs="+")
    parser.add_argument("--modalities", nargs="+", choices=MODALITIES, default=list(MODALITIES))
    parser.add_argument("--pairs-per-tissue", type=int, default=12)
    parser.add_argument("--log10-tolerance", type=float, default=.001)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.pairs_per_tissue < 1 or args.log10_tolerance <= 0:
        parser.error("Pair count and tolerance must be positive.")
    store = DataStore(args.data_dir)
    tissues = args.tissues or store.tissues
    if set(tissues) - set(store.tissues):
        parser.error("Unknown tissue ID")
    report = validate(args.data_dir, tissues, args.modalities, args.pairs_per_tissue, args.log10_tolerance)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Total: {report['tested']} tests, {report['failures']} failures")
    raise SystemExit(bool(report["failures"]))


if __name__ == "__main__":
    main()
