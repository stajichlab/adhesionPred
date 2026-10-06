import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.repeats import parse_repeat_table, repeat_rows, run_detector


def prot(pid, state="ok", length=100):
    return Protein(pid, "M" * length, "x" * 64, state, "", 0.0)


REPEATS = (
    "strain\tprotein\tlength\trep_period\trep_score\trep_start\trep_end\trep_n_copies\trep_coverage\n"
    "S\tA\t300\t47\t0.8\t80\t270\t4.0\t0.62\n"
    "S\tB\t300\t0\t0.1\t0\t0\t0\t0.0\n"
    "S\tC\t300\t12\t0.5\t10\t60\t4.1\t0.17\n"
    "S\tD\t300\t12\t0.5\t10\t60\t2.4\t0.40\n"
)


def test_repeat_call_needs_period_coverage_and_copies(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows([prot(x) for x in "ABCD"] + [prot("E")], parse_repeat_table(path))
    assert [r.get("call") for r in rows] == [
        "called",
        "not_called",
        "not_called",
        "not_called",
        None,
    ]
    assert rows[-1]["state"] == "error"  # not in the detector output


def test_repeat_thresholds_are_parameters(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows([prot("C")], parse_repeat_table(path), min_coverage=0.10)
    assert rows[0]["call"] == "called"


def test_repeat_table_needs_its_columns_and_unique_proteins(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text("protein\trep_period\nA\t3\n")
    with pytest.raises(ValueError, match="missing column"):
        parse_repeat_table(path)
    path.write_text("protein\trep_period\trep_n_copies\trep_coverage\nA\t3\t3\t0.5\nA\t3\t3\t0.5\n")
    with pytest.raises(ValueError, match="duplicate protein"):
        parse_repeat_table(path)


def test_run_detector_builds_the_command(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr("subprocess.run", lambda cmd, check: calls.append((cmd, check)))
    cmd = run_detector(
        tmp_path, "repeat14", "in.faa", "out.tsv", python="py", extra=["--mode", "exact"]
    )
    assert calls == [(cmd, True)]
    assert cmd[0] == "py" and cmd[1].endswith("analysis/cocci_repeats/14_repeat_detect_general.py")
    assert cmd[2:] == ["in.faa", "--out", "out.tsv", "--mode", "exact"]


def test_a_protein_shorter_than_the_detector_minimum_is_not_called_not_an_error(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows(
        [prot("SHORT", length=50), prot("LONG_MISSING", length=200)], parse_repeat_table(path)
    )
    assert (rows[0]["state"], rows[0]["call"]) == ("ok", "not_called")
    assert rows[1]["state"] == "error"  # long enough to be profiled, but absent from the table
