"""Build phenotype SQLite indexes from release BEDs and mapping summaries.

Run after ratgtex-server-data compilation. Transfer indexes and source files
with mtimes preserved, and use the same API code version to build and serve.
"""
import argparse
from contextlib import closing
import csv
import gzip
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import zlib

import numpy as np

from .data_store import DataStore, MODALITIES, phenotype_index_path, expression_path, source_signature


def build_index(source, destination, modality, thresholds=None, tissue=None, scale="normalized"):
    source, destination = Path(source), Path(destination)
    signature = source_signature(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=destination.parent, suffix=".sqlite3")
    os.close(fd)
    try:
        with closing(sqlite3.connect(temporary)) as conn, conn, gzip.open(source, "rt", newline="") as handle:
            conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            conn.execute("CREATE TABLE phenotype (phenotype_id TEXT PRIMARY KEY, gene_id TEXT NOT NULL, chromosome TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, data BLOB NOT NULL) WITHOUT ROWID")
            conn.execute("CREATE INDEX phenotype_gene ON phenotype(gene_id)")
            conn.execute("CREATE TABLE threshold (gene_id TEXT PRIMARY KEY, value REAL) WITHOUT ROWID")
            conn.execute("INSERT INTO metadata VALUES ('phenotype_format', 'phenotypes-v1')")
            conn.execute("INSERT INTO metadata VALUES ('modality', ?)", (modality,))
            conn.execute("INSERT INTO metadata VALUES ('tissue', ?)", (tissue,))
            conn.execute("INSERT INTO metadata VALUES ('scale', ?)", (scale,))
            if tissue is None or (scale == "normalized" and thresholds is None):
                raise ValueError("Phenotype indexes require matching tissue and threshold source.")
            if thresholds is not None:
                threshold_signature = source_signature(Path(thresholds))
                conn.execute("INSERT INTO metadata VALUES ('threshold_source', ?)", (json.dumps(threshold_signature),))
                with open(thresholds, newline="") as summary:
                    for entry in csv.DictReader(summary, delimiter="\t"):
                        if entry["tissue"] == tissue:
                            value = float(entry["pval_nominal_threshold"])
                            conn.execute("INSERT INTO threshold VALUES (?, ?)", (entry["gene_id"], value if np.isfinite(value) else None))
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            if header[:4] != ["#chr", "start", "end", "phenotype_id"]:
                raise ValueError(f"Unexpected expression header in {source}")
            samples = header[4:]
            if len(samples) != len(set(samples)):
                raise ValueError("Duplicate sample IDs")
            conn.executemany("INSERT INTO metadata VALUES (?, ?)", [
                ("source", json.dumps(signature)), ("samples", json.dumps(samples)),
                ("format", "float64-little-endian-zlib-v1"),
            ])
            count = 0
            for row in reader:
                if len(row) != len(header):
                    raise ValueError(f"Wrong number of columns on row {count + 2}")
                values = np.asarray(row[4:], dtype="<f8")
                if not row[3] or not samples or (scale == "normalized" and not np.isfinite(values).all()):
                    raise ValueError("Normalized phenotypes require IDs, samples and finite values.")
                gene = row[3] if modality in ("expression", "stability") else row[3].partition("__")[0]
                if modality not in ("expression", "stability") and ("__" not in row[3] or not gene):
                    raise ValueError(f"Expected gene__phenotype ID: {row[3]}")
                start, end = int(row[1]), int(row[2])
                if start < 0 or end <= start:
                    raise ValueError("Invalid BED coordinates")
                conn.execute("INSERT INTO phenotype VALUES (?, ?, ?, ?, ?, ?)",
                             (row[3], gene, row[0], start, end, zlib.compress(values.tobytes(), level=1)))
                count += 1
        if source_signature(source) != signature:
            raise RuntimeError("Expression source changed during indexing")
        if thresholds is not None and source_signature(Path(thresholds)) != threshold_signature:
            raise RuntimeError("Threshold source changed during indexing")
        os.chmod(temporary, 0o644)
        os.replace(temporary, destination)
        return count
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data/v4")
    parser.add_argument("--tissues", nargs="+")
    parser.add_argument("--modalities", nargs="+", choices=MODALITIES, default=list(MODALITIES))
    parser.add_argument("--scales", nargs="+", choices=("normalized", "tpm"), default=["normalized"])
    args = parser.parse_args()
    if "tpm" in args.scales and args.modalities != ["expression"]:
        parser.error("TPM is available for expression; specify --modalities expression --scales tpm.")
    store = DataStore(args.data_dir)
    for tissue in args.tissues or store.tissues:
        if tissue not in store.tissues:
            parser.error(f"Unknown tissue: {tissue}")
        for modality in args.modalities:
            for scale in args.scales:
                source = (expression_path(store.root, tissue, "tpm") if scale == "tpm" else
                          store.root / "phenos" / f"phenos.{tissue}.{modality}.norm.v4_rn8.bed.gz")
                thresholds = store.root / "xqtl" / f"top_assoc.{modality}.v4_rn8.tsv" if scale == "normalized" else None
                destination = phenotype_index_path(store.root, tissue, modality, scale)
                count = build_index(source, destination, modality, thresholds, tissue, scale)
                print(f"{tissue} {modality} {scale}: indexed {count:,} phenotypes", flush=True)


if __name__ == "__main__":
    main()
