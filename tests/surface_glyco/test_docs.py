"""Docs and changelog agree with the code."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_changelog_lists_the_breaking_rename_under_unreleased():
    text = (ROOT / "CHANGELOG.md").read_text()
    unreleased = text.split("## [Unreleased]")[1].split("\n## [")[0]
    assert "BREAKING" in unreleased and "surface_glyco" in unreleased


def test_readme_uses_current_names_only():
    text = (ROOT / "README.md").read_text()
    assert "adhesion_train" not in text and "adhesion_predict " not in text
    assert "surface_glyco_train" in text


def test_agents_does_not_claim_a_thread_safe_cache():
    assert "thread-safe" not in (ROOT / "AGENTS.md").read_text().lower()


def test_decisions_doc_records_the_2026_09_30_changes():
    text = (ROOT / "docs" / "PLAN-2026-09-30-pipeline-and-decisions.md").read_text()
    assert "Changes made later on 2026-09-30" in text
