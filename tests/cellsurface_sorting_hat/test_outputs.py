"""Unit tests of the writers, without the command line."""

import csv
import gzip

from cellsurface_sorting_hat.engine import CallRecord
from cellsurface_sorting_hat.outputs import (
    LONG_COLUMNS,
    RunInfo,
    render_report,
    write_long,
    write_wide,
)


def info(**over):
    base = {
        "version": "x",
        "taxa": {40: 2},
        "taxonomy_sha256": "t",
        "config_sha256": "c",
        "default_gate": "step1_rule@R0",
        "thresholds": {"antigen_percentile_max": 15},
        "n_proteins": 2,
        "n_invalid": 0,
        "invalid": [],
        "n_trailing_stop": 0,
        "module_states": {"step1_rule@R0": "ok"},
    }
    base.update(over)
    return RunInfo(**base)


def rec(protein, call, variant, value, status="unvalidated", basis="", other=""):
    return CallRecord(protein, call, variant, value, status, basis, other)


def read(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_long_table_has_the_documented_columns(tmp_path):
    write_long(
        tmp_path / "l.tsv.gz", [rec("A", "tandem_repeat_protein", "", "called", "smoke", "m:t")]
    )
    rows = read(tmp_path / "l.tsv.gz")
    assert rows[0] == LONG_COLUMNS
    assert rows[1] == ["A", "tandem_repeat_protein", "", "called", "smoke", "m:t", ""]


def test_wide_table_names_gated_columns_with_the_variant(tmp_path):
    records = [
        rec("A", "signal_peptide_protein", "R0", "called", "estimated"),
        rec(
            "A",
            "other_not_surface",
            "R0",
            "called",
            "estimated",
            other="cocci_specificity_rank_top15",
        ),
        rec("B", "signal_peptide_protein", "R0", "not_called"),
        rec("B", "other_not_surface", "R0", "not_assessable"),
    ]
    write_wide(tmp_path / "w.tsv.gz", records, ["A", "B"])
    rows = read(tmp_path / "w.tsv.gz")
    assert rows[0] == [
        "protein",
        "signal_peptide_protein[R0]",
        "signal_peptide_protein[R0]_status",
        "other_not_surface[R0]",
        "other_not_surface[R0]_status",
        "other_not_surface[R0]_basis",
    ]
    assert rows[1][1:3] == ["called", "estimated"] and rows[1][-1] == "cocci_specificity_rank_top15"


def test_report_warns_about_error_partial_and_inconsistent_modules():
    text = render_report(
        info(
            module_states={"a": "error", "b": "partial", "c": "ok"},
            module_notes={"a": "no result table"},
            inconsistent={"c": 2},
        ),
        [rec("A", "tandem_repeat_protein", "", "called")],
    )
    assert "WARNING: Module a is in state error. no result table" in text
    assert "WARNING: Module b is in state partial." in text
    assert "WARNING: Module c gives different rows to 2 group(s)" in text


def test_report_truncates_a_long_invalid_list():
    invalid = [[f"P{i}", "internal stop codon"] for i in range(25)]
    text = render_report(info(n_invalid=25, invalid=invalid), [rec("A", "x", "", "called")])
    assert "- P19: internal stop codon" in text and "- P20:" not in text
    assert "... and 5 more (see proteins.tsv.gz)" in text


def test_report_counts_values_per_call_and_variant():
    records = [
        rec("A", "signal_peptide_protein", "R0", v) for v in ("called", "called", "not_called")
    ]
    text = render_report(info(), records)
    assert "| signal_peptide_protein | R0 | 2 | 1 | 0 |" in text


def test_known_limits_print_the_configured_antigen_percentile():
    report = render_report(info(thresholds={"antigen_percentile_max": 10}), [])
    assert "is the top 10% of a fixed Coccidioides" in report
    assert "top 15%" not in report
