"""Stream exon and independent-association TSVs into api/catalog.sqlite3.

Run after ratgtex-server-data compilation, using the API code version that
will serve the release. Preserve source mtimes when transferring the files.
"""
import argparse
from contextlib import closing
import csv
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile

from .catalog import (ASSOCIATION_COLUMNS, CATALOG_FORMAT, CATALOG_SOURCES,
                      INTEGER_COLUMNS, REAL_COLUMNS, quote_identifier, source_signature)


def converted(value, column):
    if column in INTEGER_COLUMNS or column in REAL_COLUMNS:
        if value.strip().lower() in {"", "na", "n/a", "nan", "null", "none"}:
            return None
        if column in INTEGER_COLUMNS:
            return int(value)
        number = float(value)
        return number if math.isfinite(number) else None
    return value if value else None


def build_catalog(root, destination=None):
    root = Path(root).resolve()
    destination = Path(destination) if destination is not None else root / "api/catalog.sqlite3"
    signatures = {name: source_signature(root / name) for name in CATALOG_SOURCES.values()}
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=destination.parent, suffix=".sqlite3")
    os.close(fd)
    counts, columns_by_table = {}, {}
    try:
        with closing(sqlite3.connect(temporary)) as conn, conn:
            conn.execute("PRAGMA cache_size = -8192")
            conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            for table, filename in CATALOG_SOURCES.items():
                with (root / filename).open(newline="") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    header = next(reader)
                    if len(set(header)) != len(header):
                        raise ValueError(f"Duplicate columns in {filename}")
                    selected = [(i, name if table == "exon" else ASSOCIATION_COLUMNS[name])
                                for i, name in enumerate(header) if table == "exon" or name in ASSOCIATION_COLUMNS]
                    columns = [name for _, name in selected]
                    required = {"geneId"} if table == "exon" else {"geneId", "tissueSiteDetailId", "variantId"}
                    if table == "xqtl":
                        required.add("modality")
                    if not required.issubset(columns) or len(set(columns)) != len(columns) or "_ordinal" in columns:
                        raise ValueError(f"Invalid catalog columns in {filename}")
                    columns_by_table[table] = columns
                    definitions = [quote_identifier(name) + " " +
                                   ("INTEGER" if name in INTEGER_COLUMNS else "REAL" if name in REAL_COLUMNS else "TEXT")
                                   for name in columns]
                    conn.execute(f"CREATE TABLE {quote_identifier(table)} (_ordinal INTEGER PRIMARY KEY, {', '.join(definitions)})")
                    insert = f"INSERT INTO {quote_identifier(table)} VALUES ({','.join('?' for _ in range(len(columns) + 1))})"
                    count = 0
                    for ordinal, row in enumerate(reader):
                        if len(row) != len(header):
                            raise ValueError(f"Wrong number of columns in {filename}, row {ordinal + 2}")
                        values = [converted(row[i], name) for i, name in selected]
                        if any(not values[columns.index(name)] for name in required):
                            raise ValueError(f"Missing required identifier in {filename}, row {ordinal + 2}")
                        conn.execute(insert, [ordinal, *values])
                        count += 1
                    counts[table] = count
                conn.execute(f'CREATE INDEX {table}_gene ON "{table}" ("geneId", _ordinal)')
                if table != "exon":
                    conn.execute(f'CREATE INDEX {table}_gene_tissue ON "{table}" ("geneId", "tissueSiteDetailId", _ordinal)')
                if table == "xqtl":
                    conn.execute('CREATE INDEX xqtl_gene_modality ON xqtl ("geneId", "modality", _ordinal)')
                    conn.execute('CREATE INDEX xqtl_gene_tissue_modality ON xqtl ("geneId", "tissueSiteDetailId", "modality", _ordinal)')
            conn.executemany("INSERT INTO metadata VALUES (?, ?)", [
                ("format", CATALOG_FORMAT), ("sources", json.dumps(signatures)),
                ("columns", json.dumps(columns_by_table)), ("row_counts", json.dumps(counts)),
            ])
            conn.execute("ANALYZE")
        if any(source_signature(root / name) != signature for name, signature in signatures.items()):
            raise RuntimeError("Catalog sources changed during preparation")
        os.chmod(temporary, 0o644)
        os.replace(temporary, destination)
        return counts
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data/v4")
    args = parser.parse_args()
    counts = build_catalog(args.data_dir)
    for table, count in counts.items():
        print(f"{table}: indexed {count:,} rows")
    print(f"Catalog: {args.data_dir / 'api/catalog.sqlite3'}")


if __name__ == "__main__":
    main()
