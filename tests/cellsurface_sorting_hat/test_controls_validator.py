"""L1 control-set validator (D15, D22): rules of ``docs/agent-tasks/COMMON-RULES.md``.

The good fixture passes. Each defect is made by a small edit of a copy of it and must fail with a
message that names the file, the row (accession) and the rule. No network, no ``names.dmp`` and no
HPCC data are needed. The check of the real committed sets is at the end.
"""

import shutil
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "controls"
REPO = Path(__file__).resolve().parents[2]


def _validate(path, **kw):
    from cellsurface_sorting_hat.controls import validate_controls

    return validate_controls(path, **kw)


@pytest.fixture
def good(tmp_path):
    dest = tmp_path / "set"
    shutil.copytree(FIXTURES / "good", dest)
    return dest


def _edit(path, fn):
    path.write_text(fn(path.read_text()))


def _messages(report):
    return [str(d) for d in report.defects]


def _assert_defect(report, file, row, rule):
    hits = [
        str(d)
        for d in report.defects
        if file in str(d) and rule in str(d) and (row or "") in str(d)
    ]
    assert hits, f"no defect with {file!r}, {row!r}, {rule!r} in {_messages(report)}"
    if row:
        assert any(f"row {row}" in h for h in hits)


def test_good_fixture_passes(good):
    rep = _validate(good)
    assert rep.defects == []
    assert rep.ok


def test_committed_good_fixture_passes():
    assert _validate(FIXTURES / "good").ok


def test_size_floor_is_reported_not_failed(good):
    rep = _validate(good)
    assert rep.ok
    size = rep.size
    assert size["positive_rows"] == 2 and size["negative_rows"] == 2
    assert size["positive_clusters"] == 1 and size["negative_clusters"] == 2
    assert size["positive_rows_to_floor"] == 18
    assert size["negative_clusters_to_floor"] == 18
    assert any("floor" in n for n in rep.notes)


def test_missing_column(good):
    _edit(good / "controls.tsv", lambda t: t.replace("\tleakage\t", "\tleak\t", 1))
    rep = _validate(good)
    assert not rep.ok
    _assert_defect(rep, "controls.tsv", None, "columns")
    assert any("leakage" in m for m in _messages(rep))


def test_bad_leakage_value(good):
    _edit(good / "controls.tsv", lambda t: t.replace("\tpartial\t", "\tsome\t", 1))
    rep = _validate(good)
    _assert_defect(rep, "controls.tsv", "P2", "leakage")


def test_bad_label(good):
    _edit(good / "controls.tsv", lambda t: t.replace("\tnegative\t", "\tneg\t", 1))
    rep = _validate(good)
    _assert_defect(rep, "controls.tsv", "N1", "label")


def test_bad_tuned_value(good):
    _edit(good / "controls.tsv", lambda t: t.replace("\tyes\t", "\tmaybe\t", 1))
    _assert_defect(_validate(good), "controls.tsv", "N1", "tuned")


def test_missing_sequence(good):
    _edit(good / "sequences.faa", lambda t: t.replace(">N2 test\nMKTAYIAKQR\n", ""))
    rep = _validate(good)
    _assert_defect(rep, "sequences.faa", "N2", "sequence")


def test_row_missing_from_clusters(good):
    _edit(good / "clusters.tsv", lambda t: t.replace("N2\tC3\n", ""))
    rep = _validate(good)
    _assert_defect(rep, "clusters.tsv", "N2", "cluster")


def test_count_mismatch_in_counts_md(good):
    _edit(
        good / "counts.md",
        lambda t: t.replace("| Rows (proteins) | 4 |", "| Rows (proteins) | 5 |"),
    )
    rep = _validate(good)
    _assert_defect(rep, "counts.md", None, "counts")
    assert any("5" in m and "4" in m for m in _messages(rep))


def test_cluster_count_mismatch_per_class(good):
    _edit(
        good / "counts.md",
        lambda t: t.replace("| Clusters negative | 2 |", "| Clusters negative | 3 |"),
    )
    _assert_defect(_validate(good), "counts.md", None, "counts")


def test_duplicate_accession(good):
    def dup(text):
        lines = text.splitlines(True)
        return "".join(lines + [lines[1]])

    _edit(good / "controls.tsv", dup)
    rep = _validate(good)
    _assert_defect(rep, "controls.tsv", "P1", "duplicate")


@pytest.mark.parametrize(
    "old,new,acc,rule",
    [
        ("\tPMID:1\t", "\t\t", "P1", "source"),
        ("\tBinds host cells.\t", "\t\t", "P1", "quote"),
    ],
)
def test_empty_source_or_quote(good, old, new, acc, rule):
    _edit(good / "controls.tsv", lambda t: t.replace(old, new, 1))
    _assert_defect(_validate(good), "controls.tsv", acc, rule)


def test_bad_date(good):
    _edit(good / "controls.tsv", lambda t: t.replace("2026-10-02", "2026-13-45", 1))
    _assert_defect(_validate(good), "controls.tsv", "N1", "date")


def test_date_in_wrong_format(good):
    _edit(good / "controls.tsv", lambda t: t.replace("2026-10-01", "01/10/2026", 1))
    _assert_defect(_validate(good), "controls.tsv", "P1", "date")


def test_non_integer_taxon_id(good):
    _edit(good / "controls.tsv", lambda t: t.replace("\t746128\t", "\tabc\t", 1))
    _assert_defect(_validate(good), "controls.tsv", "P1", "taxon_id")


def test_empty_species(good):
    _edit(good / "controls.tsv", lambda t: t.replace("\tAspergillus fumigatus\t", "\t\t", 1))
    _assert_defect(_validate(good), "controls.tsv", "P1", "species")


def test_missing_file(good):
    (good / "clusters.tsv").unlink()
    _assert_defect(_validate(good), "clusters.tsv", None, "file")


def test_mixed_cluster_needs_flag(good):
    _edit(good / "clusters.tsv", lambda t: t.replace("N1\tC2", "N1\tC1"))
    rep = _validate(good)
    _assert_defect(rep, "clusters.tsv", "N1", "mixed")


def test_mixed_cluster_is_allowed_when_the_row_says_so(good):
    _edit(good / "clusters.tsv", lambda t: t.replace("N1\tC2", "N1\tC1"))
    flags = {"accession": "mixed_cluster_with_positive", "N1": "yes"}

    def add_column(text):
        lines = text.splitlines()
        return "".join(f"{line}\t{flags.get(line.split(chr(9))[0], 'no')}\n" for line in lines)

    _edit(good / "controls.tsv", add_column)
    _edit(
        good / "counts.md",
        lambda t: t.replace(
            "clusters (30% identity, 50% coverage) | 3", "clusters (30% identity, 50% coverage) | 2"
        ),
    )
    rep = _validate(good)
    assert not any("mixed" in m for m in _messages(rep)), _messages(rep)
    assert not any("controls.tsv" in m and "N1" in m for m in _messages(rep))


NAMES = "746128\t|\tAspergillus fumigatus\t|\t\t|\tscientific name\t|\n"


def test_names_dmp_is_optional_and_checked_when_given(good, tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(NAMES)
    assert _validate(good, names_dmp=names).ok
    _edit(
        good / "controls.tsv",
        lambda t: t.replace("Aspergillus fumigatus", "Aspergillus fumigata", 1),
    )
    assert _validate(good).ok
    rep = _validate(good, names_dmp=names)
    _assert_defect(rep, "controls.tsv", "P1", "names.dmp")


def test_cli_exit_codes(good, capsys):
    from cellsurface_sorting_hat.controls import main

    assert main([str(good)]) == 0
    _edit(good / "controls.tsv", lambda t: t.replace("2026-10-02", "2026-13-45", 1))
    assert main([str(good)]) == 1
    assert "row N1" in capsys.readouterr().out


def test_cli_on_committed_bad_fixture(capsys):
    from cellsurface_sorting_hat.controls import main

    assert main([str(FIXTURES / "bad")]) == 1
    assert "row N1" in capsys.readouterr().out


def test_calibrate_cli_has_the_controls_subcommand(good, capsys):
    from cellsurface_sorting_hat.calibration.cli import main

    assert main(["controls", str(good)]) == 0


# Committed sets. Known state, pinned so that a new defect or a fixed one is noticed. The data are
# not edited here: the owner decides whether to fix them.
# analysis/lectin_validation/controls: two negative chitinases (Q12713, Q4U4T0, n1_shared) are in the
# 30% cluster of the positive C1G842 and the table has no mixed_cluster_with_positive column.
# The owner decides whether to fix the data or to add the flag. Sets without controls.tsv
# (data/controls/repeat-mechanism, repeat-adhesins-other-species) use another layout and are not
# checked here. data/controls/hard-negatives-shared passes. This list must only get shorter.
KNOWN_FAILURES = {"analysis/lectin_validation/controls": 2}


def _committed_sets():
    found = [p.parent for p in (REPO / "data").rglob("controls.tsv")]
    found += [p.parent for p in (REPO / "analysis").rglob("controls.tsv")]
    return sorted(found)


def test_committed_sets_match_the_known_state():
    state = {}
    for path in _committed_sets():
        rep = _validate(path)
        if rep.defects:
            state[str(path.relative_to(REPO))] = len(rep.defects)
    assert state == KNOWN_FAILURES
