"""Make analysis/kingdom_survey importable as plain modules for tests."""

import sys
from pathlib import Path

KINGDOM_SURVEY_DIR = Path(__file__).resolve().parents[2] / "analysis" / "kingdom_survey"
sys.path.insert(0, str(KINGDOM_SURVEY_DIR))
