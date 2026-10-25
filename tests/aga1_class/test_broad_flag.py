"""Tests of analysis/aga1_class/broad_flag.py (synthetic rows) and of the real-record counts quoted in the report."""

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/aga1_class/broad_flag.py"
spec = importlib.util.spec_from_file_location("broad_flag", SCRIPT)
bf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bf)


def row(acc, sp, gpi, run, group="hard_negative", tool="aiupred"):
    return {
        "accession": acc, "role_group": group, "signal_peptide_pred": sp, "gpi_pred": gpi,
        "disorder_longest_run": str(run), "disorder_tool": tool,
    }  # fmt: skip


def flags(rows, n=200):
    return {r["accession"]: r["broad_flag"] for r in bf.flag_rows(rows, n)}


def test_all_three_parts_needed_and_threshold_is_inclusive():
    rows = [
        row("A", "yes", "yes", 200), row("B", "yes", "yes", 199), row("C", "no", "yes", 500), row("D", "yes", "no", 500),
    ]  # fmt: skip
    assert flags(rows) == {"A": "yes", "B": "no", "C": "no", "D": "no"}


def test_the_threshold_is_a_parameter():
    rows = [row("A", "yes", "yes", 150)]
    assert flags(rows, 200) == {"A": "no"} and flags(rows, 100) == {"A": "yes"}
    assert {r["broad_min_run"] for r in bf.flag_rows(rows, 100)} == {100}


def test_three_valued_logic_no_beats_na_and_na_beats_yes():
    rows = [
        row("A", "yes", "NA", 300), row("B", "no", "NA", 300), row("C", "yes", "NA", 10), row("D", "yes", "yes", "NA"),
    ]  # fmt: skip
    assert flags(rows) == {"A": "NA", "B": "no", "C": "no", "D": "NA"}


@pytest.mark.parametrize(
    "rows, n, tool, msg",
    [
        ([row("A", "yes", "yes", 5, tool="iupred3")], 200, "aiupred", "expected only 'aiupred'"),
        ([row("A", "maybe", "yes", 5)], 200, "aiupred", "expected yes, no or NA"),
        ([row("A", "yes", "yes", "x")], 200, "aiupred", "non-negative integer"),
        ([row("A", "yes", "yes", 5)], 0, "aiupred", "at least 1"),
        ([], 200, "aiupred", "no rows"),
    ],
)
def test_bad_input_stops(rows, n, tool, msg):
    with pytest.raises(bf.FlagError, match=msg):
        bf.flag_rows(rows, n, tool)


def test_missing_disorder_columns_stop():
    with pytest.raises(bf.FlagError, match="lacks columns"):
        bf.flag_rows([{"accession": "A"}], 200)


def _run(*args):
    return subprocess.run(
        [sys.executable, "-I", str(SCRIPT), *args], capture_output=True, text=True
    )


def test_cli_writes_default_200_and_refuses_overwrite(tmp_path):
    rec = tmp_path / "rec.tsv"
    cols = list(row("A", "yes", "yes", 1))
    with rec.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows([row("A", "yes", "yes", 250), row("B", "yes", "yes", 120)])
    out = tmp_path / "o.tsv"
    r = _run("--record", str(rec), "--out", str(out), "--summary-out", str(tmp_path / "s.tsv"))
    assert r.returncode == 0 and "minimum disordered run 200" in r.stdout
    got = {x["accession"]: x["broad_flag"] for x in csv.DictReader(out.open(), delimiter="\t")}
    assert got == {"A": "yes", "B": "no"}
    again = _run("--record", str(rec), "--out", str(out))
    assert again.returncode == 1 and "exists" in again.stderr and "FATAL" in again.stderr
    bad = _run("--record", str(rec), "--out", str(tmp_path / "x.tsv"), "--disorder-tool", "iupred3")
    assert bad.returncode == 1 and "FATAL" in bad.stderr and not (tmp_path / "x.tsv").exists()


REAL = ROOT / "docs/reports/data/aga1_disorder_2026-10-24/aga1_signalp_predgpi_aiupred.tsv"


def _real(n):
    with REAL.open(newline="") as fh:
        flagged = bf.flag_rows(list(csv.DictReader(fh, delimiter="\t")), n)
    return {s["role_group"]: s for s in bf.summarise(flagged)}


def test_real_counts_at_200_and_100():
    s200, s100 = _real(200), _real(100)
    assert (s200["class_positive_A"]["yes"], s200["hard_negative"]["yes"], s200["background_assumed_negative"]["yes"]) == (
        23, 8, 7,
    )  # fmt: skip
    assert (s100["class_positive_A"]["yes"], s100["hard_negative"]["yes"], s100["background_assumed_negative"]["yes"]) == (
        26, 23, 31,
    )  # fmt: skip


def test_record_hashes_match_and_tables_carry_their_setting():
    import hashlib
    import re

    rec = ROOT / "docs/reports/data/aga1_broad_flag_2026-10-24"
    found = re.findall(r"^([0-9a-f]{64})  (\S+)$", (rec / "README.md").read_text(), flags=re.M)
    assert {p for _, p in found} == {
        "broad_flag_run200.tsv", "broad_flag_run100.tsv", "summary_run200.tsv", "summary_run100.tsv",
    }  # fmt: skip
    for sha, p in found:
        assert hashlib.sha256((rec / p).read_bytes()).hexdigest() == sha, p
    for n in (200, 100):
        with (rec / f"broad_flag_run{n}.tsv").open(newline="") as fh:
            assert {r["broad_min_run"] for r in csv.DictReader(fh, delimiter="\t")} == {str(n)}
