"""Tests for analysis/adhesion_properties/db.py — uses an in-memory DuckDB
with tiny fake tables matching the real functionalDB schema, never the
real 278GB database."""

import duckdb
from db import fetch_domain_hits, fetch_lengths


def _fake_connection():
    con = duckdb.connect(":memory:")
    con.execute("CREATE TABLE gene_proteins (protein_id VARCHAR, length BIGINT, peptide VARCHAR)")
    con.execute(
        "INSERT INTO gene_proteins VALUES "
        "('LOC1_000001-T1.protein', 4, 'SSTT'), "
        "('LOC1_000002-T1.protein', 3, 'PPP')"
    )
    con.execute("CREATE TABLE pfam (protein_id VARCHAR, pfam_id VARCHAR)")
    con.execute(
        "INSERT INTO pfam VALUES ('LOC1_000001-T1', 'DomainA'), ('LOC1_000001-T1', 'DomainB')"
    )
    return con


def test_fetch_lengths_strips_protein_suffix_and_matches_ids():
    con = _fake_connection()
    result = fetch_lengths(con, ["LOC1_000001-T1", "LOC1_000002-T1"])
    result = result.sort_values("protein_id").reset_index(drop=True)
    assert list(result["protein_id"]) == ["LOC1_000001-T1", "LOC1_000002-T1"]
    assert list(result["length"]) == [4, 3]
    assert list(result["peptide"]) == ["SSTT", "PPP"]


def test_fetch_lengths_only_returns_requested_ids():
    con = _fake_connection()
    result = fetch_lengths(con, ["LOC1_000001-T1"])
    assert len(result) == 1
    assert result.iloc[0]["protein_id"] == "LOC1_000001-T1"


def test_fetch_domain_hits_returns_matching_rows():
    con = _fake_connection()
    result = fetch_domain_hits(con, "pfam", ["LOC1_000001-T1"])
    assert len(result) == 2
    assert set(result["pfam_id"]) == {"DomainA", "DomainB"}


def test_fetch_domain_hits_rejects_unknown_table():
    con = _fake_connection()
    try:
        fetch_domain_hits(con, "not_a_real_table", ["LOC1_000001-T1"])
        raise AssertionError("expected ValueError")
    except ValueError:
        pass
