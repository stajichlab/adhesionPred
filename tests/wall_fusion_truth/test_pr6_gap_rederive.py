"""The re-derivation of the round 2 verdicts under the fixed clade rule (Piece 1 of the drift-check PR)."""

import csv
import gzip
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
RECORD = ROOT / "docs/reports/data/pr6_gap_uniref50_2026-10-22"
REDERIVED = ROOT / "docs/reports/data/pr6_gap_uniref50_2026-10-22_rederived"


def load(name):
    spec = importlib.util.spec_from_file_location(name, WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


rd = load("pr6_gap_rederive")


def table(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_committed_rederived_table_equals_a_fresh_run(tmp_path):
    assert rd.main([str(RECORD), str(tmp_path / "o")]) == 0
    fresh = (tmp_path / "o/rederived_verdicts.tsv").read_text()
    assert fresh == (REDERIVED / "rederived_verdicts.tsv").read_text()


def test_only_hyd1_like_changes_and_it_becomes_iterated_only():
    rows = table(REDERIVED / "rederived_verdicts.tsv")
    assert len(rows) == 8
    changed = {r["genes"] for r in rows if r["changed"] == "yes"}
    assert changed == {"hyd1_like"}
    by = {(r["genes"], r["method"]): r for r in rows}
    assert by[("hyd1_like", "phmmer")]["status_committed"] == "buildable"
    assert by[("hyd1_like", "phmmer")]["status_fixed"] == "buildable_iterated_only"
    # the one-pass set: 5 sequences, 2 of them without a clade, so 1 named clade: not met
    p = by[("hyd1_like", "phmmer")]
    assert (p["n_outside_seq"], p["n_unassigned_seq"], p["n_named_clades_fixed"]) == ("5", "2", "1")
    assert p["rule_met_committed"] == "yes" and p["rule_met_fixed"] == "no"
    # the iterated set meets it in 3 named clades
    j = by[("hyd1_like", "jackhmmer")]
    assert (j["n_outside_seq"], j["n_unassigned_seq"], j["n_named_clades_fixed"]) == (
        "43",
        "3",
        "3",
    )
    for g, st in (
        ("aga1", "buildable_iterated_only"),
        ("ece1", "buildable_iterated_only"),
        ("ecm33", "buildable"),
    ):
        assert by[(g, "phmmer")]["status_fixed"] == st


def test_the_strict_column_is_unchanged_in_every_row():
    for r in table(REDERIVED / "rederived_verdicts.tsv"):
        assert r["strict_committed"] == r["strict_fixed"]


def test_counts_match_the_committed_hit_table():
    with gzip.open(RECORD / "gap_homologs.tsv.gz", "rt", newline="") as fh:
        hits = list(csv.DictReader(fh, delimiter="\t"))
    for r in table(REDERIVED / "rederived_verdicts.tsv"):
        mine = [
            h
            for h in hits
            if h["query"] == r["query"]
            and h["method"] == r["method"]
            and h["where"] == "outside_clade"
        ]
        assert len(mine) == int(r["n_outside_seq"])
        assert sum(1 for h in mine if not h["clade_fine"]) == int(r["n_unassigned_seq"])
        assert len({h["clade_fine"] for h in mine if h["clade_fine"]}) == int(
            r["n_named_clades_fixed"]
        )


def test_the_run_record_is_not_touched_and_refuses_a_filled_directory(tmp_path):
    out = tmp_path / "o"
    out.mkdir()
    (out / "x").write_text("keep")
    assert rd.main([str(RECORD), str(out)]) == 1
    assert (out / "x").read_text() == "keep"
    assert (
        "buildable\t" in (RECORD / "gap_per_query.tsv").read_text()
    )  # the record still says what it said


def test_no_sequence_lines_in_the_rederived_files():
    for p in REDERIVED.iterdir():
        for line in p.read_text().splitlines():
            assert not (len(line) >= 60 and line.isalpha() and line.isupper())
