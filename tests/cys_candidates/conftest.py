"""Make analysis/cys_candidates importable as a plain module for tests."""

import sys
from pathlib import Path

CYS_DIR = Path(__file__).resolve().parents[2] / "analysis" / "cys_candidates"
sys.path.insert(0, str(CYS_DIR))
