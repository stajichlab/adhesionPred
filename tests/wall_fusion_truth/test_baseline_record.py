"""The committed run record of the Task 5.4 baseline job (30233470): no sequence, the pin, no pooled row, no tier change."""

import csv
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/kiss_baseline_2026-10-22"
PIN = "9472e32e319a213211a69051f87208c8ff2a3f5f"


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _wilson(k, n, z=1.959964):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def test_no_committed_table_holds_a_sequence():
    seq = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{40,}$")
    for f in sorted(REC.glob("*")):
        if f.suffix in {".tsv", ".txt", ".md"}:
            assert not any(seq.match(line.strip()) for line in f.read_text().splitlines()), f.name


def test_provenance_carries_the_pin_and_the_pfam_release():
    text = (REC / "provenance.txt").read_text()
    assert f"commit: {PIN} (pinned: {PIN})" in text
    assert "pfam_release: 38.2" in text
    assert "models: 39 (with GA 39, without GA 0)" in text
    assert _rows("no_ga_models.tsv") == []


def test_recall_by_slice_has_no_pooled_row_and_the_wilson_interval_recomputes():
    rows = _rows("recall_by_slice.tsv")
    assert {r["slice"] for r in rows} == {"adhesion", "wall", "fusion"}
    for r in rows:
        n, k = int(r["n_clusters"]), int(r["k_clusters_any"])
        if n:
            lo, hi = _wilson(k, n)
            assert (
                abs(lo - float(r["wilson95_low"])) < 6e-4
                and abs(hi - float(r["wilson95_high"])) < 6e-4
            )
        if r["slice"] != "fusion":
            assert r["status"] == ("estimated" if n >= 5 else "descriptive" if n else "none")
        if r["slice"] == "fusion":
            assert r["status"] in {"descriptive", "none"}


def test_missed_clusters_count_and_family_reps_are_never_counted():
    missed = _rows("clusters_missed.tsv")
    by_slice = {s: sum(r["slice"] == s for r in missed) for s in ("adhesion", "wall", "fusion")}
    assert by_slice == {"adhesion": 12, "wall": 28, "fusion": 4}
    for r in _rows("family_rep_descriptive.tsv"):
        assert r["counts_toward_floor"] == "no" and r["may_be_test_control"] == "no"
    assert all(r["note"].startswith("descriptive") for r in _rows("wall_rows_false_hits.tsv"))
