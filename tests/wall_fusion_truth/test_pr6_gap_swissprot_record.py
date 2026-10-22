"""The committed run record of the PR 6 gap search, round 1 on Swiss-Prot (job 30260197)."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_gap_swissprot_2026-10-22"
PIN = "6f8aa839d10159549c7240bb0e28a65e4afac2cb"
QUERIES = {"A0ACE4KFX7": "KS0017", "P32323": "KS0080", "Q07730": "KS0107", "Q5AGC4": "KS0086"}


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_committed_table_holds_a_sequence():
    seq = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{40,}$")
    files = [f for f in sorted(REC.glob("*")) if f.suffix in {".tsv", ".txt", ".md"}]
    assert files
    for f in files:
        assert not any(seq.match(line.strip()) for line in f.read_text().splitlines()), f.name


def test_provenance_carries_the_pin_and_the_inputs():
    text = (REC / "provenance.txt").read_text()
    assert f"commit: {PIN} (pinned: {PIN})" in text
    assert "pfam_release: 38.2" in text
    assert "pfam_sha256_same_as_earlier_jobs: yes" in text
    assert "uniprot_fasta: /srv/projects/db/Uniprot/2025_03/uniprot_sprot.fasta" in text
    assert "uniprot_fasta_sequences: 573661" in text


def test_four_queries_all_not_buildable_with_the_expected_reasons():
    rows = _rows("gap_per_query.tsv")
    assert {r["query"]: r["cluster_id"] for r in rows} == QUERIES
    assert {r["status"] for r in rows} == {"not_buildable"}
    reason = {r["cluster_id"]: r["reason"] for r in rows}
    assert reason["KS0017"] == "only 1 outside sequences (minimum 5)"
    assert reason["KS0080"] == "no outside fungal homolog"
    assert reason["KS0107"] == "no outside fungal homolog"
    assert reason["KS0086"] == "only 4 outside sequences (minimum 5)"


def test_counts_recompute_from_the_hit_table():
    hits = _rows("gap_homologs.tsv")
    n = {}
    for h in hits:
        key = (h["query"], h["method"], h["where"])
        n[key] = n.get(key, 0) + 1
    assert n[("Q5AGC4", "phmmer", "outside_clade")] == 4
    assert n[("Q5AGC4", "jackhmmer", "outside_clade")] == 4
    assert n[("Q5AGC4", "jackhmmer", "non_fungal")] == 37
    assert n[("A0ACE4KFX7", "jackhmmer", "outside_clade")] == 1
    assert ("A0ACE4KFX7", "phmmer", "outside_clade") not in n
    assert not any(k[0] in {"P32323", "Q07730"} and k[2] == "outside_clade" for k in n)
    per_query = {r["query"]: r for r in _rows("gap_per_query.tsv")}
    for q, r in per_query.items():
        for m in ("phmmer", "jackhmmer"):
            assert int(r[f"{m}_outside_seq"]) == n.get((q, m, "outside_clade"), 0)
            assert int(r[f"{m}_inside_clade_seq"]) == n.get((q, m, "inside_clade"), 0)


def test_ecm33_outside_hits_are_two_clades_in_three_species():
    out = [
        h
        for h in _rows("gap_homologs.tsv")
        if h["query"] == "Q5AGC4" and h["method"] == "phmmer" and h["where"] == "outside_clade"
    ]
    assert len(out) == 4
    assert {h["clade_fine"] for h in out} == {"Eurotiomycetes", "Taphrinomycotina"}
    assert len({h["taxid"] for h in out}) == 3


def test_ecm33_check_is_below_ga_for_both_proteins():
    rows = {r["sequence"]: r for r in _rows("ecm33_check.tsv")}
    assert set(rows) == {"P38248", "Q5AGC4"}
    for r in rows.values():
        assert r["model"] == "PF12454" and r["verdict"] == "below_GA"
        assert float(r["ga_seq"]) == 23.2 and float(r["score_seq"]) < float(r["ga_seq"])
    assert float(rows["Q5AGC4"]["score_seq"]) == 5.5
    assert float(rows["P38248"]["score_seq"]) == 0.6
