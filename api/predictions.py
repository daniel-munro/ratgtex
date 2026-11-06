"""Read indexed round 11.2 predictions without converting released scores."""
from pathlib import Path
from zipfile import BadZipFile, ZipFile
from zlib import error as CompressionError


COHORT = "round11_2"
MAX_GENES = 5
MAX_RFIDS = 100
MAX_ROWS = 2000
RFID_MEMBER = "rfids.tsv"


class PredictionsUnavailableError(Exception):
    pass


class PredictionLimitError(Exception):
    pass


def archive_path(root, tissue, modality):
    return Path(root) / f"pred-{COHORT}-{tissue}-{modality}.zip"


def gene_member(gene):
    return f"genes/{gene}.tsv"


def extract_predictions(root, genes, rfids, tissues, modalities):
    rows, unavailable = [], []
    matched = set()
    encoded_rfids = [rfid.encode() for rfid in rfids]
    requested = set(encoded_rfids)
    for tissue in tissues:
        for modality in modalities:
            path = archive_path(root, tissue, modality)
            try:
                archive = ZipFile(path)
            except FileNotFoundError:
                unavailable.append({"tissue": tissue, "modality": modality})
                continue
            except (OSError, BadZipFile) as error:
                raise PredictionsUnavailableError(f"Cannot open the prediction index for {tissue}/{modality}.") from error
            try:
                with archive:
                    header = archive.read(RFID_MEMBER).rstrip(b"\r\n").split(b"\t")
                    positions = {value: index + 1 for index, value in enumerate(header) if value in requested}
                    matched.update(positions)
                    columns = [positions.get(rfid) for rfid in encoded_rfids]
                    for gene in genes:
                        try:
                            handle = archive.open(gene_member(gene))
                        except KeyError:
                            continue
                        with handle:
                            for line in handle:
                                if len(rows) == MAX_ROWS:
                                    raise PredictionLimitError(
                                        f"Result exceeds {MAX_ROWS} phenotype rows. Select fewer genes, tissues or modalities."
                                    )
                                fields = line.rstrip(b"\r\n").split(b"\t")
                                rows.append({
                                    "geneId": gene,
                                    "tissue": tissue,
                                    "modality": modality,
                                    "phenotypeId": fields[0].decode(),
                                    "values": [fields[col].decode() if col is not None else None for col in columns],
                                })
            except (OSError, BadZipFile, CompressionError, EOFError, KeyError, IndexError, UnicodeError) as error:
                raise PredictionsUnavailableError(f"Cannot read the prediction index for {tissue}/{modality}.") from error
    if len(unavailable) == len(tissues) * len(modalities):
        raise PredictionsUnavailableError("No prediction indexes are available for the selected tissues and modalities.")
    represented = {row["geneId"] for row in rows}
    return {
        "cohort": COHORT,
        "rfids": rfids,
        "unmatchedRfids": [rfid for rfid, encoded in zip(rfids, encoded_rfids) if encoded not in matched],
        "genesWithoutPredictions": [gene for gene in genes if gene not in represented],
        "unavailableDatasets": unavailable,
        "predictions": rows,
    }
