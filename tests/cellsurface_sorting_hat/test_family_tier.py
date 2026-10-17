"""Task 3.3 (D3, D11, D19): the family tier is derived by code from a committed measurement file.

The tier describes the family model's own evidence. The call status describes the call's
measurement. Neither implies the other.
"""

import csv
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "data" / "sorting_hat" / "family_table.tsv"
MEASURE = ROOT / "data" / "sorting_hat" / "family_measure.tsv"
XFAIL_33 = pytest.mark.xfail(strict=True, reason="Task 3.3: family_tier module")
XFAIL_34 = pytest.mark.xfail(strict=True, reason="Task 3.4: tier column in family_table.tsv")


def row(**over):
    """A per-family measurement row that derives ``estimated``: 25/25, 22/22, narrow intervals."""
    base = {
        "pfam_acc": "PF00001",
        "name": "X",
        "measure_basis": "per_family_reviewed",
        "denominator": "family_members",
        "n_pos": "25",
        "pos_hits": "",
        "pos_denominator": "",
        "coverage": "",
        "n_neg": "25",
        "n_clusters_pos": "22",
        "n_clusters_neg": "22",
        "neg_hits": "1",
        "neg_clusters_hit": "1",
        "negative_basis": "curated",
        "sensitivity": "0.9",
        "sens_lo": "0.85",
        "sens_hi": "0.95",
        "specificity": "0.9",
        "spec_lo": "0.85",
        "spec_hi": "0.95",
        "n_pos_in_reference": "0",
        "n_pos_unknown": "0",
    }
    base.update(over)
    return base


# ---- derivation ---------------------------------------------------------------------------------


@XFAIL_33
def test_a_reviewed_per_family_measure_with_enough_data_is_estimated():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row()) == "estimated"


@XFAIL_33
def test_the_same_measure_on_shared_descriptive_sets_is_capped_at_smoke():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row(measure_basis="descriptive_shared_sets")) == "smoke"


@XFAIL_33
def test_five_and_five_derive_smoke():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row(n_pos="5", n_neg="5", n_clusters_pos="5", n_clusters_neg="5")) == "smoke"


@XFAIL_33
def test_no_row_derives_untested():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(None) == "untested"


@XFAIL_33
def test_no_positive_or_no_negative_derives_untested():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row(n_pos="0")) == "untested"
    assert derive_tier(row(n_neg="0")) == "untested"


@XFAIL_33
def test_a_row_without_n_clusters_pos_derives_at_most_smoke():
    """FD-A3: the E1 positives are not clustered, so the floor must not be skipped."""
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row(n_clusters_pos="")) == "smoke"
    assert derive_tier(row(n_clusters_neg="")) == "smoke"


@XFAIL_33
def test_negatives_that_are_not_curated_derive_at_most_smoke():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row(negative_basis="assumed")) == "smoke"


@XFAIL_33
def test_positives_that_are_mostly_in_reference_or_unknown_cap_the_tier_at_smoke():
    from cellsurface_sorting_hat.family_tier import derive_tier

    assert derive_tier(row(n_pos_in_reference="20")) == "smoke"
    assert derive_tier(row(n_pos_unknown="20")) == "smoke"
    assert derive_tier(row(n_pos_in_reference="6", n_pos_unknown="6")) == "estimated"


@XFAIL_33
def test_an_unknown_measure_basis_is_refused():
    from cellsurface_sorting_hat.family_tier import derive_tier

    with pytest.raises(ValueError, match="measure_basis"):
        derive_tier(row(measure_basis="vibes"))


@XFAIL_33
def test_a_descriptive_row_counts_the_hits_as_positives():
    """A descriptive row has no member list of its own: the E1 hits are its positives (D19)."""
    from cellsurface_sorting_hat.family_tier import derive_tier

    descriptive = row(
        measure_basis="descriptive_shared_sets",
        denominator="e1_adhesins_all_classes",
        n_pos="",
        pos_hits="9",
        pos_denominator="54",
        coverage="0.167",
        sensitivity="",
        sens_lo="",
        sens_hi="",
        specificity="",
        spec_lo="",
        spec_hi="",
        n_neg="137",
        n_clusters_pos="",
        n_clusters_neg="76",
    )
    assert derive_tier(descriptive) == "smoke"
    assert derive_tier({**descriptive, "pos_hits": "0"}) == "untested"  # Hyr1: no E1 hit


# ---- the committed measurement file -------------------------------------------------------------


def _read(path):
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader), reader.fieldnames


@XFAIL_33
def test_the_measurement_file_uses_coverage_not_sensitivity_for_the_shared_sets():
    rows, header = _read(MEASURE)
    assert "coverage" in header and "measure_basis" in header
    assert {r["measure_basis"] for r in rows} == {"descriptive_shared_sets"}
    assert all(r["sensitivity"] == "" for r in rows)  # D19: coverage is never sensitivity
    assert len(rows) == 24


@XFAIL_33
def test_the_committed_measures_derive_no_tier_above_smoke_and_hyr1_is_untested():
    from cellsurface_sorting_hat.family_tier import derive_tiers

    tiers = derive_tiers(TABLE, MEASURE)
    assert set(tiers.values()) <= {"untested", "smoke"}
    by_name = {r["name"]: tiers[r["pfam_acc"]] for r in _read(TABLE)[0]}
    assert by_name["Hyr1"] == "untested"
    assert by_name["PIR1-like_C"] == "untested"  # 0 of 54 E1 positives
    assert by_name["CFEM"] == "smoke"
    assert by_name["PA14"] == "smoke"


@XFAIL_34
def test_the_committed_table_equals_the_derived_tiers():
    from cellsurface_sorting_hat.family_tier import check

    assert check(TABLE, MEASURE) == []


@XFAIL_34
def test_a_table_with_one_stronger_tier_fails(tmp_path):
    from cellsurface_sorting_hat.family_tier import check

    rows, header = _read(TABLE)
    for r in rows:
        if r["name"] == "CFEM":
            r["specificity_tier"] = "estimated"
    path = tmp_path / "family_table.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    problems = check(path, MEASURE)
    assert len(problems) == 1 and "PF05730" in problems[0] and "estimated" in problems[0]


@XFAIL_33
def test_write_changes_only_the_tier_column(tmp_path):
    from cellsurface_sorting_hat.family_tier import main

    header = ["pfam_acc", "name", "specificity_tier", "note"]
    table = tmp_path / "t.tsv"
    table.write_text(
        "\t".join(header) + "\nPF05730\tCFEM\testimated\tkeep me\nPF99999\tOther\tsmoke\tand me\n"
    )
    measure = tmp_path / "m.tsv"
    with open(MEASURE, encoding="utf-8", newline="") as src, open(measure, "w", newline="") as dst:
        dst.write(src.read())
    assert main(["--check", "--table", str(table), "--measure", str(measure)]) == 1
    assert main(["--write", "--table", str(table), "--measure", str(measure)]) == 0
    rows = _read(table)[0]
    assert [(r["pfam_acc"], r["specificity_tier"], r["note"]) for r in rows] == [
        ("PF05730", "smoke", "keep me"),
        ("PF99999", "untested", "and me"),
    ]
    assert main(["--check", "--table", str(table), "--measure", str(measure)]) == 0


# ---- code location (D15) ------------------------------------------------------------------------


def test_no_test_of_the_blocking_suite_imports_from_analysis():
    pattern = re.compile(r"^\s*(from|import)\s+analysis\b|sys\.path[^\n]*analysis", re.MULTILINE)
    offenders = [
        str(p.relative_to(ROOT))
        for p in Path(__file__).parent.rglob("*.py")
        if pattern.search(p.read_text(encoding="utf-8"))
    ]
    assert offenders == []


def test_the_tier_code_is_in_src_not_in_analysis():
    assert not (ROOT / "analysis" / "family_tests" / "derive_tier.py").exists()
