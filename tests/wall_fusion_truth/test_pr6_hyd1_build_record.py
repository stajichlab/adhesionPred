"""The committed run record of the PR 6 hyd1_like build (jobs 30321490, 30325384 and 30335767)."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_hyd1_like_2026-10-23"
PIN = "9f4843b1f7d3111c1f2c5db8544076b78d2b8a67"
FREEZE_SHA = "65ed9bb6dbcc6035a4093537a5afc0f92a179d3f4356902c1aba83a12e552ff2"


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_committed_file_holds_a_sequence():
    seq = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{40,}")
    suffixes = {".tsv", ".txt", ".md", ".json", ".log"}
    files = [f for f in sorted(REC.rglob("*")) if f.suffix in suffixes]
    assert len(files) > 100
    for f in files:
        for line in f.read_text().splitlines():
            assert not seq.search(line), f
    assert not list(REC.rglob("*.aln.fa")) and not list(REC.rglob("*.hmm"))
    assert not list(REC.rglob("*.faa")) and not list(REC.rglob("*.fa"))


def test_provenance_carries_the_pin_and_the_freeze_hash():
    paths = list(REC.glob("build/*/*/provenance.txt"))
    assert len(paths) == 9
    for p in paths:
        text = p.read_text()
        assert f"commit {PIN} (pinned {PIN})" in text, p
        assert f"freeze sha256: {FREEZE_SHA}" in text, p
    assert f"commit {PIN}" in (REC / "extract/provenance.txt").read_text()
    for p in REC.glob("build/*/*/cutoff.json"):
        assert FREEZE_SHA in p.read_text(), p


def test_nine_unit_mode_rows_and_sordariomycetes_cluster30_not_buildable():
    rows = _rows("report/recall_by_unit.tsv")
    assert len(rows) == 9
    assert {r["independence_mode"] for r in rows if r["independence_mode"]} == {"cluster30", "none"}
    bad = [
        r
        for r in rows
        if r["unit"] == "clade:Sordariomycetes"
        and r["label"] == "descriptive"
        and r["rule_ok"] == "not_buildable"
    ]
    assert len(bad) == 1
    assert (REC / "build/clade_Sordariomycetes/cluster30/manifest.json").exists()
    assert (
        '"status": "not_buildable"'
        in (REC / "build/clade_Sordariomycetes/cluster30/manifest.json").read_text()
    )
    assert not (REC / "build/clade_Sordariomycetes/cluster30/cutoff.json").exists()


def test_no_unit_is_testable_and_u_is_at_most_two():
    built = [r for r in _rows("report/recall_by_unit.tsv") if r["rule_ok"] != "not_buildable"]
    assert len(built) == 8
    for r in built:
        assert int(r["u_indep"]) <= 2, r["unit"]
        assert int(r["min_u_freeze"]) == 3
        assert r["testable_u"] == "no" and r["rule_ok"] == "not_testable"
        assert r["label"] in {"descriptive", "partial"}
    assert sum(1 for r in built if r["testable_u"] == "yes") == 0


def test_the_kiss_unit_is_recovered():
    kiss = [r for r in _rows("report/recall_by_unit.tsv") if r["unit"] == "clade:Eurotiomycetes"]
    assert len(kiss) == 2
    for r in kiss:
        assert r["n_held_out"] == "1" and r["custom_hit"] == "1"
    scores = {}
    for mode in ("cluster30", "none"):
        text = (REC / f"build/clade_Eurotiomycetes/{mode}/hits_custom.tsv").read_text()
        scores[mode] = float(text.splitlines()[1].split("\t")[1])
    assert scores["cluster30"] > 50 and scores["none"] > 70


def test_panel_ignored_is_nineteen_ids():
    rows = _rows("report/recall_by_unit.tsv")
    assert {r["panel_ignored_n"] for r in rows if r["panel_ignored_n"]} == {"19"}
    import json

    cut = json.loads((REC / "build/all/none/cutoff.json").read_text())
    assert len(cut["panel_ignored"]) == 19
    assert "B6QUR1" in cut["panel_ignored"]
    assert (REC / "negcheck/proteome_panel_ga.tsv").read_text().count("\n") == 20
