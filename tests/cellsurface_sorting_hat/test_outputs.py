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


def test_run_json_carries_the_module_identities(tmp_path):
    import json

    from cellsurface_sorting_hat.outputs import write_run_json

    identities = [{"name": "repeat02", "version": "1", "params_hash": "p", "artefact_hash": "a"}]
    write_run_json(tmp_path / "run.json", info(module_identities=identities))
    assert json.loads((tmp_path / "run.json").read_text())["module_identities"] == identities
    write_run_json(tmp_path / "run2.json", info())
    assert json.loads((tmp_path / "run2.json").read_text())["module_identities"] == []


def _spec_limits():
    """The numbered known-limits items of the orchestrator spec, joined to one line each."""
    import re
    from pathlib import Path

    text = (
        Path(__file__).parents[2] / "docs/superpowers/specs/2026-10-04-orchestrator-design.md"
    ).read_text()
    block = text.split("Known limits, printed in the report header", 1)[1].split("\n### ", 1)[0]
    items = []
    for line in block.splitlines()[1:]:
        if re.match(r"^\d+\. ", line):
            items.append(re.sub(r"^\d+\. ", "", line).strip())
        elif line.startswith("   ") and items:
            items[-1] += " " + line.strip()
    return items


def test_the_spec_and_the_report_list_the_same_known_limits():
    from cellsurface_sorting_hat.outputs import _known_limits

    code = [" ".join(t.split()) for t in _known_limits({"antigen_percentile_max": 15})]
    spec = _spec_limits()
    assert len(spec) == len(code) == 9
    for i, (a, b) in enumerate(zip(spec, code, strict=True), 1):
        if i != 5:  # the spec writes the percentile as P%; the code prints the configured value
            assert a == b, f"item {i}"


def test_the_last_known_limit_says_the_signalp_leakage_was_not_measured():
    report = render_report(info(), [])
    assert (
        "9. Leakage: overlap between the Phase C positives and the SignalP 6 training data was "
        "not measured." in report
    )


def test_report_counts_the_surface_attachment_basis_by_sub_label():
    records = [
        rec("a", "surface_attachment_candidate", "R0", "called", other="hydrophobin_domain"),
        rec(
            "b",
            "surface_attachment_candidate",
            "R0",
            "called",
            other="tandem_repeat_protein,hsba_domain",
        ),
        rec("c", "surface_attachment_candidate", "R0", "called", other="hydrophobin_domain"),
        rec("d", "surface_attachment_candidate", "R0", "not_called"),
    ]
    text = render_report(info(), records)
    assert "## Surface attachment: which evidence held" in text
    assert "| tandem repeat | 1 |" in text
    assert "| hydrophobin | 2 |" in text
    assert "| HsbA | 1 |" in text
    assert "| wall family domain (PA14, CFEM) | 0 |" in text
    assert "| proteins called | 3 |" in text
    assert "adhesin repeat" not in text


def test_report_counts_a_protein_once_when_two_variants_call_it():
    records = [
        rec("a", "surface_attachment_candidate", "R0", "called", other="hydrophobin_domain"),
        rec("a", "surface_attachment_candidate", "R1", "called", other="hydrophobin_domain"),
    ]
    text = render_report(info(), records)
    assert "| proteins called | 1 |" in text
    assert "| hydrophobin | 1 |" in text


def test_wide_table_has_a_basis_column_for_surface_attachment_candidate(tmp_path):
    records = [
        rec(
            "A",
            "surface_attachment_candidate",
            "R0",
            "called",
            other="tandem_repeat_protein,hsba_domain",
        )
    ]
    write_wide(tmp_path / "w.tsv.gz", records, ["A"])
    rows = read(tmp_path / "w.tsv.gz")
    assert rows[0][-1] == "surface_attachment_candidate[R0]_basis"
    assert rows[1][-1] == "tandem_repeat_protein,hsba_domain"


def test_report_defines_the_three_statuses_and_says_composites_are_derived():
    text = render_report(info(), [])
    assert "## What the statuses mean" in text
    for word in ("`unvalidated`", "`smoke`", "`estimated`"):
        assert word in text
    assert "weakest status of the leaf calls that decided" in text


def test_report_has_no_attachment_section_without_called_proteins():
    text = render_report(info(), [rec("a", "surface_attachment_candidate", "R0", "not_called")])
    assert "Surface attachment: which evidence held" not in text
