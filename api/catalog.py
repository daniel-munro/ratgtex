"""Read-only, gene-scoped queries against the prepared release catalog."""
from contextlib import closing
import json
import sqlite3


CATALOG_FORMAT = "exon-association-catalog-v1"
CATALOG_SOURCES = {
    "exon": "exon.v4.tsv",
    "eqtl": "eqtl/eqtls_indep.v4_rn8.tsv",
    "xqtl": "xqtl/xqtls_indep.cross_modality.v4_rn8.tsv",
}
ASSOCIATION_COLUMNS = {
    "tissue": "tissueSiteDetailId", "gene_id": "geneId",
    "phenotype_id": "phenotypeId", "variant_id": "variantId",
    "pval_beta": "pValueBeta", "log2_aFC": "log2aFC", "chrom": "chromosome",
    "pos": "pos", "af": "altAlleleFrequency", "rank": "rank",
    "tss_distance": "tssDistance", "pval_nominal": "pValue",
    "slope": "slope", "slope_se": "slopeSE", "qval": "qValue",
    "modality": "modality", "ref": "ref", "alt": "alt",
}
INTEGER_COLUMNS = {"start", "end", "exonNumber", "pos", "rank", "tssDistance"}
REAL_COLUMNS = {"pValue", "pValueBeta", "log2aFC", "altAlleleFrequency", "slope", "slopeSE", "qValue"}


class CatalogUnavailableError(ValueError):
    """A current catalog is required to serve the requested collection."""


def source_signature(path):
    stat = path.stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def quote_identifier(name):
    return '"' + name.replace('"', '""') + '"'


def query_catalog(root, collection, gene, tissues=None, modality=None, limit=None, offset=0):
    if collection not in CATALOG_SOURCES:
        raise ValueError("Unknown catalog collection")
    path = root / "api/catalog.sqlite3"
    message = "Missing, invalid or stale API catalog; run python -m api.prepare_catalog."
    if not path.is_file():
        raise CatalogUnavailableError(message)
    try:
        with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as conn:
            conn.row_factory = sqlite3.Row
            # This bounds SQLite's page cache, independently of the response cache.
            conn.execute("PRAGMA cache_size = -2048")
            meta = dict(conn.execute("SELECT key, value FROM metadata"))
            source = CATALOG_SOURCES[collection]
            if (meta.get("format") != CATALOG_FORMAT
                    or json.loads(meta["sources"])[source] != source_signature(root / source)):
                raise CatalogUnavailableError(message)
            columns = json.loads(meta["columns"])[collection]
            if (not isinstance(columns, list) or not columns
                    or any(not isinstance(name, str) or name == "_ordinal" for name in columns)):
                raise CatalogUnavailableError(message)
            where, parameters = ['"geneId" = ?'], [gene]
            if tissues is not None:
                if not tissues:
                    where.append("0")
                else:
                    where.append('"tissueSiteDetailId" IN (' + ",".join("?" for _ in tissues) + ")")
                    parameters.extend(tissues)
            if modality is not None:
                where.append('"modality" = ?')
                parameters.append(modality)
            table = quote_identifier(collection)
            condition = " AND ".join(where)
            sql = f"SELECT {', '.join(map(quote_identifier, columns))} FROM {table} WHERE {condition} ORDER BY _ordinal"
            total = None
            if limit is not None:
                total = conn.execute(f"SELECT count(*) FROM {table} WHERE {condition}", parameters).fetchone()[0]
                sql += " LIMIT ? OFFSET ?"
                parameters.extend([limit, offset])
            rows = [dict(row) for row in conn.execute(sql, parameters)]
            return rows, len(rows) if total is None else total
    except (sqlite3.Error, KeyError, TypeError, json.JSONDecodeError, FileNotFoundError) as error:
        raise CatalogUnavailableError(message) from error
