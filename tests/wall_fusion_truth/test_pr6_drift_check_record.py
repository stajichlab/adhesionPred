"""The committed run record of the PR 6 drift check (job 30282220) and its sensitivity re-reading."""

import csv
import re
import statistics
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_drift_check_2026-10-22"
PIN = "e333a66565fdd21c0fc95743958b08d013a5efed"
NAMES = {"A0ACE4KFX7": "hyd1_like", "P32323": "aga1", "Q07730": "ece1"}
VERDICT = {
    "A0ACE4KFX7": "supported_rule_met",
    "P32323": "supported_rule_met",
    "Q07730": "not_supported",
}
NO_CLADE = "(no clade)"


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _query_cover_supported(row):
    """The same E-value and null thresholds as the freeze, with QUERY coverage instead of max(query, hit)."""
    if not row["best_domain_ievalue"] or float(row["best_domain_ievalue"]) > 1e-5:
        return False
    if row["null_ratio"] == "" or float(row["null_ratio"]) >= 0.5:
        return False
    return float(row["query_cover"]) >= 0.5


def test_no_committed_file_holds_a_sequence():
    seq = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{40,}$")
    files = [f for f in sorted(REC.glob("*")) if f.suffix in {".tsv", ".txt", ".md"}]
    assert {f.name for f in files} >= {
        "drift_hits.tsv",
        "drift_per_clade.tsv",
        "drift_summary.tsv",
        "hits_manifest.tsv",
        "provenance.txt",
        "summary.txt",
        "README.md",
    }
    for f in files:
        assert not any(seq.match(line.strip()) for line in f.read_text().splitlines()), f.name
    assert not any(seq.match(r["hit"]) for r in _rows("drift_hits.tsv"))


def test_provenance_carries_the_pin_and_the_target():
    text = (REC / "provenance.txt").read_text()
    assert f"commit: {PIN} (pinned: {PIN})" in text
    assert "uniprot_fasta: /srv/projects/db/Uniprot/2025_03/uniref50.fasta" in text
    assert "uniprot_fasta_bytes: 28781193344" in text
    assert "queries: P32323 Q07730 A0ACE4KFX7" in text
    assert "hit_sequences_extracted: 233" in text
    assert "phmmer -E 10 --domE 10" in text and "shuffle seed 1" in text


def test_three_queries_and_their_verdicts():
    rows = _rows("drift_summary.tsv")
    assert {r["query"]: r["verdict"] for r in rows} == VERDICT
    assert {r["query"]: int(r["n_jackhmmer_outside"]) for r in rows} == {
        "A0ACE4KFX7": 43,
        "P32323": 64,
        "Q07730": 126,
    }
    assert len(_rows("hits_manifest.tsv")) == len(_rows("drift_hits.tsv")) == 233


def test_counts_recompute_from_the_hit_table():
    hits = _rows("drift_hits.tsv")
    summary = {r["query"]: r for r in _rows("drift_summary.tsv")}
    for q, row in summary.items():
        mine = [h for h in hits if h["query"] == q]
        n = Counter(h["class"] for h in mine)
        assert int(row["n_supported"]) == n["supported_homolog"]
        assert int(row["n_bias_driven"]) == n["bias_driven"] == 0
        assert int(row["n_unsupported"]) == n["unsupported"]
        named = {
            h["clade_fine"]
            for h in mine
            if h["class"] == "supported_homolog" and h["clade_fine"] != NO_CLADE
        }
        assert int(row["supported_named_clades"]) == len(named)
    assert [
        Counter(h["class"] for h in hits if h["query"] == q)["supported_homolog"] for q in NAMES
    ] == [32, 11, 0]


def test_aga1_supported_hits_pass_only_through_the_lenient_coverage_rule():
    hits = _rows("drift_hits.tsv")
    sup = [h for h in hits if h["query"] == "P32323" and h["class"] == "supported_homolog"]
    qc = [float(h["query_cover"]) for h in sup]
    hc = [float(h["hit_cover"]) for h in sup]
    assert len(sup) == 11
    assert (min(qc), max(qc)) == (0.148, 0.31)
    assert round(statistics.median(qc), 3) == 0.217
    assert max(qc) < 0.5 and min(hc) >= 0.5
    assert (min(hc), max(hc)) == (0.549, 0.946)
    assert float(sup[0]["query_stnqgp"]) == 0.644 and float(sup[0]["query_low_complexity"]) == 0.67


def test_hyd1_like_supported_hits_all_have_query_coverage_of_at_least_half():
    sup = [
        h
        for h in _rows("drift_hits.tsv")
        if h["query"] == "A0ACE4KFX7" and h["class"] == "supported_homolog"
    ]
    qc = [float(h["query_cover"]) for h in sup]
    assert len(sup) == 32
    assert (min(qc), max(qc)) == (0.531, 0.99)
    assert round(statistics.median(qc), 3) == 0.724


def test_the_shuffled_null_scored_zero_for_every_hit_with_a_pairwise_hit():
    hits = _rows("drift_hits.tsv")
    with_null = Counter(h["query"] for h in hits if h["null_ratio"] != "")
    assert dict(with_null) == {"A0ACE4KFX7": 32, "P32323": 20}
    assert all(float(h["null_score"]) == 0 for h in hits if h["null_ratio"] != "")


def test_sensitivity_table_with_query_coverage():
    hits = _rows("drift_hits.tsv")
    out = {}
    for q in NAMES:
        sup = [h for h in hits if h["query"] == q and _query_cover_supported(h)]
        named = {h["clade_fine"] for h in sup if h["clade_fine"] != NO_CLADE}
        verdict = "supported_rule_met" if len(sup) >= 5 and len(named) >= 2 else "not_supported"
        out[q] = (len(sup), len(named), verdict)
    assert out == {
        "A0ACE4KFX7": (32, 3, "supported_rule_met"),
        "P32323": (0, 0, "not_supported"),
        "Q07730": (0, 0, "not_supported"),
    }


def test_the_report_states_the_sensitivity_and_changes_nothing():
    text = (ROOT / "docs/reports/2026-10-22-pr6-drift-check.md").read_text()
    assert "not a re-classification and not a new freeze" in text
    assert "owner decision" in text
    assert "drift is confirmed" in text
