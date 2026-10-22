"""The committed run record of the PR 6 gap search, round 2 on UniRef50 (job 30262455)."""

import csv
import gzip
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_gap_uniref50_2026-10-22"
PIN = "e11808b2bf297cc0aa1dfb345dfb3bf568a2c937"
QUERIES = {"A0ACE4KFX7": "KS0017", "P32323": "KS0080", "Q07730": "KS0107", "Q5AGC4": "KS0086"}
STATUS = {
    "KS0017": "buildable",
    "KS0080": "buildable_iterated_only",
    "KS0107": "buildable_iterated_only",
    "KS0086": "buildable",
}


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _hits():
    with gzip.open(REC / "gap_homologs.tsv.gz", "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_committed_table_holds_a_sequence():
    seq = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{40,}$")
    files = [f for f in sorted(REC.glob("*")) if f.suffix in {".tsv", ".txt", ".md"}]
    assert files
    for f in files:
        assert not any(seq.match(line.strip()) for line in f.read_text().splitlines()), f.name
    lines = _hits()
    assert lines and set(lines[0]) == {
        "query",
        "method",
        "hit_accession",
        "evalue",
        "score",
        "taxid",
        "clade_fine",
        "where",
    }
    assert not any(seq.match(h["hit_accession"]) for h in lines)


def test_provenance_carries_the_pin_and_the_uniref50_target():
    text = (REC / "provenance.txt").read_text()
    assert f"commit: {PIN} (pinned: {PIN})" in text
    assert "pfam_release: 38.2" in text
    assert "pfam_sha256_same_as_earlier_jobs: yes" in text
    assert "uniprot_fasta: /srv/projects/db/Uniprot/2025_03/uniref50.fasta" in text
    assert "uniprot_fasta_sequences: 70198728" in text
    assert "uniprot_fasta_bytes: 28781193344" in text
    assert "settings: phmmer -E 1e-5; jackhmmer -N 3 --incE 1e-5 --incdomE 1e-5" in text


def test_four_queries_and_their_statuses():
    rows = _rows("gap_per_query.tsv")
    assert {r["query"]: r["cluster_id"] for r in rows} == QUERIES
    assert {r["cluster_id"]: r["status"] for r in rows} == STATUS
    assert "buildable 2; buildable_iterated_only 2" in (REC / "summary.txt").read_text()


def test_counts_recompute_from_the_hit_table():
    n = Counter((h["query"], h["method"], h["where"]) for h in _hits())
    per_query = {r["query"]: r for r in _rows("gap_per_query.tsv")}
    for q, r in per_query.items():
        for m in ("phmmer", "jackhmmer"):
            assert int(r[f"{m}_outside_seq"]) == n.get((q, m, "outside_clade"), 0)
            assert int(r[f"{m}_inside_clade_seq"]) == n.get((q, m, "inside_clade"), 0)
            assert int(r[f"{m}_non_fungal_seq"]) == n.get((q, m, "non_fungal"), 0)
    assert n[("Q5AGC4", "phmmer", "outside_clade")] == 379
    assert n[("Q5AGC4", "jackhmmer", "outside_clade")] == 953
    assert n[("A0ACE4KFX7", "phmmer", "outside_clade")] == 5
    assert n[("A0ACE4KFX7", "jackhmmer", "outside_clade")] == 43
    assert n[("P32323", "phmmer", "outside_clade")] == 4
    assert n[("P32323", "jackhmmer", "outside_clade")] == 64
    assert ("Q07730", "phmmer", "outside_clade") not in n
    assert n[("Q07730", "jackhmmer", "outside_clade")] == 126


def test_the_drift_flag_rests_on_thousands_of_non_fungal_hits():
    n = Counter((h["query"], h["method"], h["where"]) for h in _hits())
    assert n[("P32323", "jackhmmer", "non_fungal")] == 2688
    assert n[("Q07730", "jackhmmer", "non_fungal")] == 6896
    for q in ("P32323", "Q07730"):
        assert n[(q, "jackhmmer", "non_fungal")] > 20 * n[(q, "jackhmmer", "outside_clade")]


def test_hyd1_like_one_pass_rule_rests_on_the_no_clade_bucket():
    out = [
        h
        for h in _hits()
        if h["query"] == "A0ACE4KFX7" and h["method"] == "phmmer" and h["where"] == "outside_clade"
    ]
    clades = Counter(h["clade_fine"] for h in out)
    assert len(out) == 5
    assert clades == {"Sordariomycetes": 3, "": 2}
    assert len([c for c in clades if c]) == 1
    # the committed clade table shows the bucket the script counts as a clade
    bucket = [
        r
        for r in _rows("gap_clades.tsv")
        if r["query"] == "A0ACE4KFX7"
        and r["method"] == "phmmer"
        and r["clade_fine"] == "(no clade)"
    ]
    assert [int(r["n_outside_seq"]) for r in bucket] == [2]


def test_ecm33_outside_hits_and_the_below_ga_check():
    out = [
        h
        for h in _hits()
        if h["query"] == "Q5AGC4" and h["method"] == "phmmer" and h["where"] == "outside_clade"
    ]
    named = Counter(h["clade_fine"] for h in out if h["clade_fine"])
    assert len(out) == 379
    assert len(named) == 18
    assert sum(1 for v in named.values() if v >= 5) == 10
    assert named["Dothideomycetes"] == 106
    rows = {r["sequence"]: r for r in _rows("ecm33_check.tsv")}
    assert set(rows) == {"P38248", "Q5AGC4"}
    assert rows["Q5AGC4"]["score_seq"] == "5.5" and rows["P38248"]["score_seq"] == "0.6"
    for r in rows.values():
        assert r["model"] == "PF12454" and r["ga_seq"] == "23.2" and r["verdict"] == "below_GA"
