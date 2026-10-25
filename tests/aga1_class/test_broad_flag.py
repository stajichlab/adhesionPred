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
    r = _run(
        "--record",
        str(rec),
        "--out",
        str(out),
        "--no-stalk-rule",
        "--summary-out",
        str(tmp_path / "s.tsv"),
    )
    assert r.returncode == 0 and "minimum disordered run 200" in r.stdout
    got = {x["accession"]: x["broad_flag"] for x in csv.DictReader(out.open(), delimiter="\t")}
    assert got == {"A": "yes", "B": "no"}
    again = _run("--record", str(rec), "--out", str(out), "--no-stalk-rule")
    assert again.returncode == 1 and "exists" in again.stderr and "FATAL" in again.stderr
    bad = _run(
        "--record",
        str(rec),
        "--out",
        str(tmp_path / "x.tsv"),
        "--no-stalk-rule",
        "--disorder-tool",
        "iupred3",
    )
    assert bad.returncode == 1 and "FATAL" in bad.stderr and not (tmp_path / "x.tsv").exists()


REC = ROOT / "docs/reports/data/aga1_broad_flag_2026-10-24"
RECORD = ROOT / "docs/reports/data/aga1_disorder_2026-10-24/aga1_signalp_predgpi_aiupred.tsv"


def _table(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _real(min_run, stalk=True, rescue=None):
    st = {r["accession"]: r for r in _table(REC / "stalk_features_w30_f0.5.tsv")} if stalk else None
    flagged = bf.flag_rows(_table(RECORD), min_run, "aiupred", st, rescue)
    return {s["role_group"]: s for s in bf.summarise(flagged)}


def _yes(summary, *groups):
    return tuple(summary[g]["yes"] for g in groups)


GROUPS = ("class_positive_A", "class_positive_B", "hard_negative", "background_assumed_negative")


def test_real_counts_without_and_with_the_stalk_rule():
    assert _yes(_real(200, stalk=False), *GROUPS) == (23, 18, 8, 7)
    assert _yes(_real(100, stalk=False), *GROUPS) == (26, 23, 23, 31)
    assert _yes(_real(200), *GROUPS) == (
        18,
        18,
        8,
        7,
    )  # the rule costs five tier A proteins, removes no negative
    assert _yes(_real(100), *GROUPS) == (22, 22, 20, 22)


def test_real_counts_with_the_exploratory_tail_rescue():
    s = _real(200, rescue=1.4)
    assert _yes(s, *GROUPS) == (19, 18, 8, 7) and s["adhesin_other_mechanism"]["yes"] == 8
    assert _yes(_real(200, rescue=2.0), *GROUPS) == (19, 18, 8, 7)


def test_the_grid_agrees_with_the_flag_tables():
    g = _table(REC / "stalk_grid.tsv")

    def yes(n, w, f, group):
        return next(
            int(r["yes"]) for r in g
            if (r["min_run"], r["stalk_window"], r["stalk_min_frac"], r["role_group"]) == (str(n), w, f, group)
        )  # fmt: skip

    assert tuple(yes(200, "none", "none", x) for x in GROUPS) == (23, 18, 8, 7)
    assert tuple(yes(200, "30", "0.5", x) for x in GROUPS) == (18, 18, 8, 7)
    assert tuple(yes(100, "30", "0.5", x) for x in GROUPS) == (22, 22, 20, 22)
    assert len(g) == 2 * 10 * len({r["role_group"] for r in g})


def test_every_committed_table_carries_its_setting():
    for p in REC.glob("broad_flag_*.tsv"):
        rows = _table(p)
        assert len(rows) == 850
        assert (
            len({r["broad_min_run"] for r in rows}) == 1
            and len({r["stalk_rule"] for r in rows}) == 1
        )
        assert ("nostalk" in p.name) == (rows[0]["stalk_rule"] == "none")


def test_record_hashes_match_every_committed_file():
    import hashlib
    import re

    found = {
        p: s
        for s, p in re.findall(
            r"^([0-9a-f]{64})  (\S+)$", (REC / "README.md").read_text(), flags=re.M
        )
    }
    assert set(found) == {p.name for p in REC.iterdir() if p.name != "README.md"}
    for p, sha in found.items():
        assert hashlib.sha256((REC / p).read_bytes()).hexdigest() == sha, p


# ---- stalk_features.py and the stalk-aware flag (synthetic) ----
sspec = importlib.util.spec_from_file_location(
    "stalk_features", ROOT / "analysis/aga1_class/stalk_features.py"
)
sf = importlib.util.module_from_spec(sspec)
sspec.loader.exec_module(sf)


def test_stalk_window_must_lie_wholly_inside_one_run():
    res = "A" * 5 + "ST" * 15 + "A" * 5  # 30 Ser/Thr-rich residues, with a flanking A on each side
    on = [0.1] * 5 + [0.9] * 30 + [0.1] * 5
    assert sf.stalk(res, on, 30, 0.5, 0.5) == (30, 1.0)
    split = (
        on[:20] + [0.1] + on[21:]
    )  # one low-scoring residue splits the run in two parts shorter than 30
    assert sf.stalk(res, split, 30, 0.5, 0.5) == (0, 0.0)
    assert sf.stalk(res, on, 31, 0.5, 0.5) == (0, 0.0)


def test_stalk_needs_the_ser_thr_fraction_and_reports_the_longest_holding_run():
    sc = [0.9] * 70 + [0.1] + [0.9] * 40  # a run of 70, one low residue, a run of 40
    # first run: 15 S then A (a 30-window of exactly 15/30 = 0.5 counts); second run all S
    assert sf.stalk("S" * 15 + "A" * 55 + "G" + "S" * 40, sc, 30, 0.5, 0.5)[0] == 70
    assert (
        sf.stalk("A" * 70 + "G" + "S" * 40, sc, 30, 0.5, 0.5)[0] == 40
    )  # only the short run holds a window
    assert sf.stalk("S" * 14 + "A" * 56 + "G" + "A" * 40, sc, 30, 0.5, 0.5) == (0, 14 / 30)


def test_ctail_is_the_best_15_window_in_the_last_40_and_na_when_short():
    assert sf.ctail("A" * 20 + "I" * 15 + "S" * 5) == ("4.500", 5)
    assert sf.ctail("ACDE") == ("NA", "NA")


def test_read_raw_reads_hash_headers_and_stops_on_bad_rows(tmp_path):
    ok = tmp_path / "r.txt"
    ok.write_text("#  banner\n#>P1 x\n1\tA\t0.1\n2\tS\t0.9\n\n#>P2\n1\tT\t0.5\n")
    assert {k: "".join(v[0]) for k, v in sf.read_raw(ok).items()} == {"P1": "AS", "P2": "T"}
    for text, msg in (
        ("#>P1\n2\tA\t0.1\n", "out of order"),
        ("#>P1\n1\tA\n", "expected 'pos aa score'"),
        ("#>P1\n1\tA\t0.1\n#>P1\n1\tA\t0.1\n", "duplicate"),
        ("1\tA\t0.1\n", "before any record"),
    ):
        ok.write_text(text)
        with pytest.raises(sf.StalkError, match=msg):
            sf.read_raw(ok)


def test_features_stop_when_accessions_differ_from_the_manifest(tmp_path):
    man = tmp_path / "m.txt"
    man.write_text("P1\nP2\n")
    with pytest.raises(sf.StalkError, match="differ from the manifest"):
        sf.features({"P1": (["A"], [0.1])}, man, 30, 0.5, 0.5)
    with pytest.raises(sf.StalkError, match="out of range"):
        sf.features({"P1": (["A"], [0.1]), "P2": (["A"], [0.1])}, man, 30, 1.5, 0.5)


def _stalk_row(acc, run, hydro="0.0"):
    return {
        "accession": acc, "stalk_window": "30", "stalk_min_frac": "0.5", "stalk_threshold": "0.5", "stalk_run_len": str(run),
        "ctail_hydro": hydro,
    }  # fmt: skip


def test_flag_uses_the_stalk_run_instead_of_the_longest_run():
    rows = [row("A", "yes", "yes", 500), row("B", "yes", "yes", 100)]
    st = {"A": _stalk_row("A", 150), "B": _stalk_row("B", 250)}
    out = {r["accession"]: r for r in bf.flag_rows(rows, 200, stalk=st)}
    assert out["A"]["broad_flag"] == "no" and out["B"]["broad_flag"] == "yes"
    assert out["A"]["stalk_rule"] == "window=30;min_st_frac=0.5;disorder_threshold=0.5"


def test_tail_rescue_only_for_signal_peptide_plus_gpi_no_and_needs_the_stalk_table():
    rows = [
        row("A", "yes", "no", 300), row("B", "no", "no", 300), row("C", "yes", "NA", 300), row("D", "yes", "no", 300),
    ]  # fmt: skip
    st = {
        k: _stalk_row(k, 300, h)
        for k, h in (("A", "2.0"), ("B", "2.0"), ("C", "2.0"), ("D", "1.0"))
    }
    out = {r["accession"]: r for r in bf.flag_rows(rows, 200, stalk=st, rescue=1.4)}
    assert [out[k]["broad_flag"] for k in "ABCD"] == ["yes", "no", "NA", "no"]
    assert [out[k]["gpi_rescued"] for k in "ABCD"] == ["yes", "no", "no", "no"]
    assert out["A"]["gpi_pred"] == "no"  # the original call stays visible
    with pytest.raises(bf.FlagError, match="needs a stalk table"):
        bf.flag_rows(rows, 200, rescue=1.4)


def test_stalk_table_must_match_the_record_accessions_and_settings():
    rows = [row("A", "yes", "yes", 300)]
    with pytest.raises(bf.FlagError, match="same accessions"):
        bf.flag_rows(rows, 200, stalk={"B": _stalk_row("B", 1)})
    mixed = {"A": _stalk_row("A", 1), "B": {**_stalk_row("B", 1), "stalk_window": "20"}}
    with pytest.raises(bf.FlagError, match="mixes several"):
        bf.flag_rows([row("A", "yes", "yes", 1), row("B", "yes", "yes", 1)], 200, stalk=mixed)


def test_cli_requires_an_explicit_stalk_choice(tmp_path):
    r = _run("--record", str(RECORD), "--out", str(tmp_path / "o.tsv"))
    assert r.returncode == 2 and "stalk" in r.stderr
