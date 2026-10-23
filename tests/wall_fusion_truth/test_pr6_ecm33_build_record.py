"""The committed run record of the PR 6 ecm33 build (jobs 30285974 and 30287890)."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_ecm33_build_2026-10-23"
PIN = "ba15e0c83ed52b2636ef6837ad51e0edf3d4c698"
FREEZE_SHA = "a7553076d4963e4daedf3d7fdf0eba28b8e8349847b2644482fc2e854d7ce915"


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
            assert not seq.fullmatch(line.strip()), f
    assert not list(REC.rglob("*.aln.fa")) and not list(REC.rglob("*.hmm"))


def test_provenance_carries_the_pin_and_the_freeze_hash():
    paths = list(REC.glob("ecm33_build/*/*/provenance.txt"))
    assert len(paths) == 23
    for p in paths:
        text = p.read_text()
        assert f"commit: {PIN} (pinned: {PIN})" in text, p
        assert f"freeze sha256: {FREEZE_SHA}" in text, p
    assert f"commit: {PIN}" in (REC / "extract_provenance.txt").read_text()
    assert FREEZE_SHA in (REC / "report/summary.txt").read_text()


def test_eleven_units_and_twenty_two_rows():
    rows = _rows("report/recall_by_unit.tsv")
    assert len(rows) == 22
    assert len({r["unit"] for r in rows}) == 11
    assert {r["independence_mode"] for r in rows} == {"cluster30", "none"}
    assert len(_rows("report/side_by_side.tsv")) == 11


def test_custom_model_hits_every_held_out_member_and_u_is_zero():
    rows = _rows("report/recall_by_unit.tsv")
    for r in rows:
        assert r["custom_hit"] == r["n_held_out"], r["unit"]
        assert r["iteration_hit"] == r["n_held_out"], r["unit"]
        assert r["u_missed_by_every_comparator"] == "0", r["unit"]
        assert r["testable_u"] == "no" and r["rule_ok"] == "not_testable"
        assert r["status"] == "descriptive"
    assert sum(1 for r in rows if r["testable_u"] == "yes") == 0
    assert sum(int(r["n_held_out"]) for r in rows if r["independence_mode"] == "cluster30") == 469


def test_extra_calls_and_hard_negatives_are_within_the_frozen_limits():
    extra = _rows("report/extra_calls.tsv")
    assert len(extra) == 22
    assert all(r["within_limit"] == "yes" for r in extra)
    assert max(int(r["extra_calls_total"]) for r in extra) == 5
    hard = _rows("report/hard_negative_calls.tsv")
    assert len(hard) == 198
    assert all(float(r["frozen_limit"]) == 0.05 for r in hard)
    assert all(float(r["call_rate"]) <= 0.05 for r in hard)
    assert all(r["within_limit"] == "yes" for r in hard)


def test_ignored_negatives_are_hit_far_above_the_cutoff():
    rows = _rows("report/ignored_negatives_line.tsv")
    assert len(rows) == 44
    assert {r["accession"] for r in rows} == {"P38248", "Q4WNS8"}
    for r in rows:
        assert r["custom_hit"] == "yes"
        assert float(r["custom_score"]) > 7 * float(r["custom_cutoff"])
        if r["accession"] == "P38248":
            assert r["pfam_relaxed_score"] == "" and r["pfam_relaxed_hit"] == "no"
        else:
            assert r["pfam_relaxed_hit"] == "yes"
