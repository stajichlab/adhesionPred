"""Task 3.8: the new rules are written down once in the rules file, and CLASSES.md says what a reader needs."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
RULES = (ROOT / "docs" / "paper" / "03-status-and-validation-rules.md").read_text(encoding="utf-8")
CLASSES = (ROOT / "docs" / "CLASSES.md").read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "term",
    [
        "active",
        "specificity_tier",
        "narrow_ok",
        "not_applicable",
        "derived_from",
        "negative_basis",
        "leakage_exception",
        "primary_label",
        "measure_basis",
        "mechanism_without_class",
    ],
)
def test_the_rules_file_names_each_new_term(term):
    section = RULES.split("## 8. Decisions of 2026-10-17")[1].split("## 7.")[0]
    assert f"`{term}" in section


def test_classes_md_separates_the_family_tier_from_the_call_status():
    assert "neither implies the other" in CLASSES
    assert "`active` means enabled by the owner, not tested" in CLASSES


def test_classes_md_has_the_interim_sentence_until_the_status_files_are_current():
    assert "Interim (2026-10-17)" in CLASSES and "stale" in CLASSES
