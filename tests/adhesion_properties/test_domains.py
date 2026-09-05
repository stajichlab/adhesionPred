"""Tests for analysis/adhesion_properties/domains.py — uses an in-memory
DuckDB with tiny fake tables, never the real 278GB database."""

import duckdb
from domains import clean_cazy_family, compute_domain_flags, top_domain_table


def _fake_connection():
    con = duckdb.connect(":memory:")
    con.execute("CREATE TABLE pfam (protein_id VARCHAR, pfam_id VARCHAR)")
    con.execute(
        "INSERT INTO pfam VALUES "
        "('A1', 'DomainX'), ('A1', 'DomainX'), ('A2', 'DomainY'), ('B1', 'DomainX'), "
        # Extra distinct proteins so DomainX clears MIN_ADHESION_COUNT_FOR_ENRICHMENT (5)
        # in the enrichment-table tests below, while A1 keeps its duplicate rows so the
        # dedup-by-protein behavior is still exercised (6 raw rows -> 5 distinct proteins).
        "('A2', 'DomainX'), ('A3', 'DomainX'), ('A4', 'DomainX'), ('A5', 'DomainX')"
    )
    con.execute("CREATE TABLE cazy_overview (protein_id VARCHAR, cazyme_fam VARCHAR)")
    con.execute(
        "INSERT INTO cazy_overview VALUES "
        "('A1', 'GH114(59-261)'), ('A2', 'GH114(1-100)'), ('A3', 'GH114(2-200)'), "
        "('A4', 'GH114(3-300)'), ('A5', 'GH114(4-400)')"
    )
    con.execute("CREATE TABLE merops (protein_id VARCHAR, merops_id VARCHAR)")
    con.execute("INSERT INTO merops VALUES ('A2', 'MER0000001')")
    con.execute("CREATE TABLE signalp (protein_id VARCHAR, probability DOUBLE)")
    con.execute("INSERT INTO signalp VALUES ('A1', 0.99)")
    con.execute("CREATE TABLE tmhmm (protein_id VARCHAR, PredHel BIGINT)")
    con.execute("INSERT INTO tmhmm VALUES ('A1', 0), ('A2', 2)")
    con.execute("CREATE TABLE targetp (protein_id VARCHAR, prediction VARCHAR)")
    con.execute("INSERT INTO targetp VALUES ('A1', 'SP')")
    return con


def test_clean_cazy_family_strips_range_suffix():
    assert clean_cazy_family("GH114(59-261)") == "GH114"
    assert clean_cazy_family("GT2_Chitin_synth_2(661-1184)") == "GT2_Chitin_synth_2"
    assert clean_cazy_family("GH20") == "GH20"


def test_compute_domain_flags():
    con = _fake_connection()
    flags = compute_domain_flags(con, ["A1", "A2", "B2"]).set_index("protein_id")
    assert bool(flags.loc["A1", "has_pfam"]) is True
    assert bool(flags.loc["A1", "has_cazy"]) is True
    assert bool(flags.loc["A1", "has_merops"]) is False
    assert bool(flags.loc["A1", "has_signal_peptide"]) is True
    assert bool(flags.loc["A1", "has_tm_helix"]) is False  # PredHel=0
    assert flags.loc["A1", "targetp_category"] == "SP"
    assert bool(flags.loc["A2", "has_tm_helix"]) is True  # PredHel=2
    assert flags.loc["A2", "targetp_category"] == "none"
    assert bool(flags.loc["B2", "has_pfam"]) is False


def test_top_domain_table_counts_distinct_proteins_not_raw_rows():
    con = _fake_connection()
    # A1 has 2 pfam rows for DomainX (6 raw rows total among A1-A5 for DomainX)
    # but distinct-protein count must be 5, not the raw row count of 6.
    top = top_domain_table(
        con,
        "pfam",
        "pfam_id",
        adhesion_ids=["A1", "A2", "A3", "A4", "A5"],
        background_ids=["B1"],
        top_n=5,
    )
    domain_x = top.set_index("domain_id").loc["DomainX"]
    assert domain_x["adhesion_count"] == 5  # distinct proteins, not 6 raw rows
    assert domain_x["background_count"] == 1  # B1


def test_top_domain_table_cleans_cazy_family_names():
    con = _fake_connection()
    top = top_domain_table(
        con,
        "cazy_overview",
        "cazyme_fam",
        adhesion_ids=["A1", "A2", "A3", "A4", "A5"],
        background_ids=[],
        top_n=5,
    )
    assert "GH114" in set(top["domain_id"])
