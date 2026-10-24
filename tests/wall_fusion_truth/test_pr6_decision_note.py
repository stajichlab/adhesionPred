"""The PR 6 decision note exists, has its sections, and the Ecm33 model is untouched."""

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NOTE = ROOT / "docs/reports/2026-10-24-pr6-decision-note.md"
MODEL_DIR = ROOT / "data/custom_models/ecm33_all_2026-10-23"
MODEL_SHA = "9a6e2c812540b0b1892f081611bf4044ec674bc524d68f5c6fc3d16001c50937"

SECTIONS = [
    "## Words used in this note",
    "## 1. Where PR 6 stands",
    "## 2. Decided versus open",
    "## 3. Every family in play",
    "## 4. Options per open question",
    "## 5. Recommended order of work",
    "## 6. What we do NOT know",
]


def test_note_exists_with_sections():
    text = NOTE.read_text()
    for heading in SECTIONS:
        assert heading in text, heading


def test_note_recommends_but_does_not_decide():
    assert "DECIDES NOTHING" in NOTE.read_text()


def test_readme_owner_decision_paragraph():
    text = (MODEL_DIR / "README.md").read_text()
    assert "## Owner decision 2026-10-24" in text
    assert "report-only" in text and "additional development" in text


def test_model_hash_unchanged():
    digest = hashlib.sha256((MODEL_DIR / "model.hmm").read_bytes()).hexdigest()
    assert digest == MODEL_SHA
    assert MODEL_SHA in (MODEL_DIR / "SHA256SUMS").read_text()
