"""Build optional per-gene expression indexes without loading whole matrices.

Run locally or on the HPC: python -m api.prepare_expression --help
The existing release files remain the source of truth; publish indexes alongside
those files, preserving their mtimes. Each index is replaced atomically.
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

from .data_store import DataStore, SCALES, expression_path, source_signature


def build_index(source, destination):
    source, destination = Path(source), Path(destination)
    signature = source_signature(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=destination.parent, suffix=".sqlite3")
    os.close(fd)
    try:
        with closing(sqlite3.connect(temporary)) as conn, conn, gzip.open(source, "rt", newline="") as handle:
            conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            conn.execute("CREATE TABLE expression (gene_id TEXT PRIMARY KEY, data BLOB NOT NULL) WITHOUT ROWID")
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
                conn.execute("INSERT INTO expression VALUES (?, ?)", (row[3], zlib.compress(values.tobytes(), level=1)))
                count += 1
        if source_signature(source) != signature:
            raise RuntimeError("Expression source changed during indexing")
        os.chmod(temporary, 0o644)
        os.replace(temporary, destination)
        return count
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data/v4")
    parser.add_argument("--tissues", nargs="+", help="Default: all tissues")
    parser.add_argument("--scales", nargs="+", choices=SCALES, default=list(SCALES), help="Default: both TPM and normalized expression")
    args = parser.parse_args()
    store = DataStore(args.data_dir)
    tissues = args.tissues or store.tissues
    for tissue in tissues:
        if tissue not in store.tissues:
            parser.error(f"Unknown tissue: {tissue}")
        for scale in args.scales:
            source = expression_path(store.root, tissue, scale)
            destination = store.root / "api/expression" / f"{tissue}.{scale}.sqlite3"
            count = build_index(source, destination)
            print(f"{tissue} {scale}: indexed {count:,} genes", flush=True)


if __name__ == "__main__":
    main()
