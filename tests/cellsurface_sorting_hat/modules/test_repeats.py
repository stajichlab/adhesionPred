import subprocess

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.base import invalid_row
from cellsurface_sorting_hat.modules.repeats import (
    DETECTOR_MIN_LEN,
    parse_repeat_table,
    repeat_rows,
    run_detector,
)


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
    rows = repeat_rows([prot(x) for x in "ABCD"], parse_repeat_table(path), min_coverage=0.10)
    assert rows[2]["call"] == "called"


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
    monkeypatch.setattr("subprocess.run", lambda cmd, check, **kw: calls.append((cmd, check, kw)))
    (tmp_path / "out.tsv").write_text("x")
    cmd = run_detector(
        tmp_path, "repeat14", "in.faa", tmp_path / "out.tsv", python="py", extra=["--mode", "exact"]
    )
    assert [(c, ch) for c, ch, _ in calls] == [(cmd, True)]
    assert calls[0][2].get("shell") is not True
    assert cmd[0] == "py" and cmd[1].endswith("analysis/cocci_repeats/14_repeat_detect_general.py")
    assert cmd[2:] == [
        "in.faa",
        "--out",
        str(tmp_path / "out.tsv"),
        "--min-len",
        str(DETECTOR_MIN_LEN),
        "--mode",
        "exact",
    ]


def test_run_detector_errors_name_module_and_fasta(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="repeat99"):
        run_detector(tmp_path, "repeat99", "in.faa", "o.tsv")

    def fail(cmd, check, **kw):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr("subprocess.run", fail)
    with pytest.raises(RuntimeError, match=r"repeat02 failed on in\.faa.*IndexError"):
        run_detector(tmp_path, "repeat02", "in.faa", tmp_path / "o.tsv")

    monkeypatch.setattr("subprocess.run", lambda cmd, check, **kw: None)
    with pytest.raises(RuntimeError, match=r"repeat02 ran on in\.faa but wrote no output"):
        run_detector(tmp_path, "repeat02", "in.faa", tmp_path / "missing.tsv")


def test_a_protein_shorter_than_the_detector_minimum_is_not_called_not_an_error(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    parsed = parse_repeat_table(path)
    rows = repeat_rows(
        [prot(x) for x in "ABCD"] + [prot("SHORT", length=50), prot("LONG_MISSING", length=200)],
        parsed,
    )
    assert (rows[4]["state"], rows[4]["call"]) == ("ok", "not_called")
    assert rows[5]["state"] == "error"  # long enough to be profiled, but absent from the table


def table(*rows):
    head = "protein\trep_period\trep_n_copies\trep_coverage\n"
    return head + "".join("\t".join(map(str, r)) + "\n" for r in rows)


def test_thresholds_are_inclusive_and_period_zero_is_never_called(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(table(("COV", 5, 3.0, 0.25), ("COP", 5, 2.5, 0.30), ("P0", 0, 4, 0.5)))
    rows = repeat_rows([prot(x) for x in ("COV", "COP", "P0")], parse_repeat_table(path))
    assert [r["call"] for r in rows] == ["called", "called", "not_called"]


def test_length_boundary_is_79_not_called_80_uses_the_rule(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(table(("L80", 5, 4.0, 0.5)))
    rows = repeat_rows([prot("L79", length=79), prot("L80", length=80)], parse_repeat_table(path))
    assert (rows[0]["state"], rows[0]["call"], rows[0]["period"]) == ("ok", "not_called", 0)
    assert (rows[1]["state"], rows[1]["call"]) == ("ok", "called")
    # 79 aa and absent from the table is not an error, 80 aa and absent is
    assert repeat_rows([prot("L80", length=80)], {})[0]["state"] == "error"


def test_trailing_stars_do_not_count_towards_the_length(tmp_path):
    p = Protein("S", "M" * 79 + "*", "x" * 64, "ok", "", 0.0)
    assert repeat_rows([p], {})[0]["call"] == "not_called"


def test_a_table_id_not_in_the_fasta_is_an_error(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(table(("A", 5, 4.0, 0.5), ("OTHER", 5, 4.0, 0.5)))
    with pytest.raises(ValueError, match=r"1 repeat-table ID\(s\) are not in the FASTA.*'OTHER'"):
        repeat_rows([prot("A")], parse_repeat_table(path))


def test_non_numeric_field_names_file_and_line(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(table(("A", 5, 4.0, 0.5), ("B", 5, "abc", 0.5)))
    with pytest.raises(ValueError, match=r"r\.tsv:3: bad numeric field"):
        parse_repeat_table(path)
    path.write_text(table(("A", "inf", 4.0, 0.5)))
    with pytest.raises(ValueError, match=r"r\.tsv:2: bad numeric field"):
        parse_repeat_table(path)


def test_invalid_protein_gets_invalid_row(tmp_path):
    p = prot("BAD", state="invalid")
    assert repeat_rows([p], {}) == [invalid_row(p)]


@pytest.mark.parametrize(
    ("row", "what"),
    [
        ("S\tA\t300\t47\t0.8\t80\t270\t4.0", "short row"),  # rep_coverage is None
        ("S\tA\t300\t47\t0.8\t80\t270\tnan\t0.5", "not finite"),
        ("S\tA\t300\t47\t0.8\t80\t270\tinf\t0.5", "not finite"),
        ("S\tA\t300\t47\t0.8\t80\t270\t4.0\t-0.5", "negative"),
        ("S\tA\t300\t-47\t0.8\t80\t270\t4.0\t0.5", "negative"),
    ],
)
def test_parse_repeat_table_refuses_a_short_row_and_bad_numbers_with_path_and_line(
    tmp_path, row, what
):
    path = tmp_path / "r.tsv"
    head = "strain\tprotein\tlength\trep_period\trep_score\trep_start\trep_end\trep_n_copies\trep_coverage\n"
    path.write_text(head + row + "\n")
    with pytest.raises(ValueError, match=rf"{path}:2.*{what}"):
        parse_repeat_table(path)


def test_parse_repeat_table_still_reads_empty_values_as_zero(tmp_path):
    path = tmp_path / "r.tsv"
    head = "strain\tprotein\tlength\trep_period\trep_score\trep_start\trep_end\trep_n_copies\trep_coverage\n"
    path.write_text(head + "S\tA\t300\t\t0.1\t0\t0\t\t\n")
    assert parse_repeat_table(path)["A"] == {"period": 0, "copies": 0.0, "coverage": 0.0}
