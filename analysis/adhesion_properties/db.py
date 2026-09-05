"""DuckDB query helpers for the functional annotation database.

functionalDB/function.duckdb is 278GB and opened read-only; every
function here filters against a caller-supplied protein-id list
(registered as a temp view) rather than scanning a full table.
"""

from collections.abc import Iterable
from pathlib import Path

import duckdb
import pandas as pd

FUNCTION_DB_PATH = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb")

_ALLOWED_DOMAIN_TABLES = {"pfam", "cazy_overview", "merops", "signalp", "tmhmm", "targetp"}


def connect(db_path: Path = FUNCTION_DB_PATH) -> "duckdb.DuckDBPyConnection":
    """Open a read-only connection to the functional annotation database."""
    return duckdb.connect(str(db_path), read_only=True)


def fetch_lengths(con, protein_ids: Iterable[str]) -> pd.DataFrame:
    """Return protein_id (without the .protein suffix) + length + peptide
    for the given ids, joined against gene_proteins.

    gene_proteins.protein_id carries a '.protein' suffix (confirmed during
    design review); this strips it so the returned protein_id matches the
    convention every other table and this analysis's own tables use.
    """
    ids_df = pd.DataFrame({"protein_id": list(protein_ids)})
    con.register("_wanted_ids", ids_df)
    query = """
        SELECT
            replace(gp.protein_id, '.protein', '') AS protein_id,
            gp.length AS length,
            gp.peptide AS peptide
        FROM gene_proteins gp
        JOIN _wanted_ids w ON gp.protein_id = w.protein_id || '.protein'
    """
    result = con.execute(query).fetchdf()
    con.unregister("_wanted_ids")
    return result


def fetch_domain_hits(
    con, table: str, protein_ids: Iterable[str], id_col: str = "protein_id"
) -> pd.DataFrame:
    """Return all rows from a functional-annotation table for the given
    protein ids. No dedup — callers decide how to collapse multi-row-
    per-protein tables (pfam) vs. single-row tables (tmhmm/targetp/signalp).

    `table` is restricted to a fixed allow-list (never interpolated from
    external input) since it's spliced directly into the SQL statement.
    """
    if table not in _ALLOWED_DOMAIN_TABLES:
        raise ValueError(f"table must be one of {_ALLOWED_DOMAIN_TABLES}, got {table!r}")

    ids_df = pd.DataFrame({id_col: list(protein_ids)})
    con.register("_wanted_ids", ids_df)
    query = f"SELECT t.* FROM {table} t JOIN _wanted_ids w ON t.{id_col} = w.{id_col}"  # noqa: S608
    result = con.execute(query).fetchdf()
    con.unregister("_wanted_ids")
    return result
