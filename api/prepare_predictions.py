"""Build per-tissue and modality ZIP indexes from round 11.2 prediction downloads."""
import argparse
import gzip
from itertools import product
import os
from pathlib import Path
import tempfile
from zipfile import ZIP_DEFLATED, ZipFile

from .data_store import DataStore, MODALITIES
from .predictions import COHORT, RFID_MEMBER, archive_path, gene_member


def build_index(source, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".prediction-", delete=False) as temporary:
        staged = Path(temporary.name)
    try:
        with (
            gzip.open(source, "rb") as handle,
            tempfile.TemporaryFile(dir=destination.parent) as spool,
            ZipFile(staged, "w", ZIP_DEFLATED) as archive,
        ):
            header = handle.readline().rstrip(b"\r\n").split(b"\t")
            if header[0] != b"phenotype_id" or len(header) < 2:
                raise ValueError(f"{source}: expected phenotype_id followed by RFID columns")
            rfids = header[1:]
            if any(not rfid for rfid in rfids) or len(set(rfids)) != len(rfids):
                raise ValueError(f"{source}: empty or duplicate RFID columns")
            archive.writestr(RFID_MEMBER, b"\t".join(rfids) + b"\n")
            offsets = {}
            phenotypes = set()
            # A gene's rows can be interleaved. Keep offsets in RAM and stage
            # only this matrix on disk, preserving the original score strings.
            for number, line in enumerate(handle, start=2):
                phenotype = line.partition(b"\t")[0]
                gene = phenotype.partition(b"__")[0].decode()
                if not gene or line.count(b"\t") != len(rfids) or phenotype in phenotypes:
                    raise ValueError(f"{source}: invalid width, empty gene or duplicate phenotype at row {number}")
                phenotypes.add(phenotype)
                offsets.setdefault(gene, []).append(spool.tell())
                spool.write(line)
            for gene, positions in offsets.items():
                with archive.open(gene_member(gene), "w") as member:
                    for offset in positions:
                        spool.seek(offset)
                        member.write(spool.readline())
        staged.chmod(0o644)
        staged.replace(destination)
        return len(phenotypes)
    finally:
        staged.unlink(missing_ok=True)


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=repo / "data/predict")
    parser.add_argument("--data-dir", type=Path, default=Path(os.environ.get("RATGTEX_DATA_DIR", repo / "data/v4")),
                        help="v4 release directory containing tissue_info.v4.tsv")
    parser.add_argument("--output-dir", type=Path, default=os.environ.get("RATGTEX_PREDICTION_DIR"))
    parser.add_argument("--tissues", nargs="+")
    parser.add_argument("--modalities", nargs="+", choices=MODALITIES)
    parser.add_argument("--allow-partial", action="store_true", help="Skip missing source files for a pilot build")
    args = parser.parse_args()
    tissues = DataStore(args.data_dir).tissues
    for tissue in args.tissues or tissues:
        if tissue not in tissues:
            parser.error(f"Unknown tissue: {tissue}")
    output_dir = args.output_dir or args.data_dir / "api/predictions"
    inputs, missing = [], []
    for tissue, modality in product(args.tissues or tissues, args.modalities or MODALITIES):
        source = args.source_dir / f"pred-{COHORT}-{tissue}-{modality}.tsv.gz"
        if source.is_file():
            inputs.append((source, archive_path(output_dir, tissue, modality)))
        else:
            missing.append(source.name)
    if missing and not args.allow_partial:
        parser.error("Missing prediction files: " + ", ".join(missing))
    if not inputs:
        parser.error("No matching prediction files found.")
    if missing:
        print(f"Skipping {len(missing)} missing prediction files.", flush=True)
    for source, destination in inputs:
        count = build_index(source, destination)
        print(f"{destination}: {count:,} rows, {destination.stat().st_size:,} bytes", flush=True)


if __name__ == "__main__":
    main()
