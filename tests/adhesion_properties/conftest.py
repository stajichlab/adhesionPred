"""Make analysis/adhesion_properties importable as plain modules for tests."""

import sys
from pathlib import Path

ADHESION_PROPERTIES_DIR = Path(__file__).resolve().parents[2] / "analysis" / "adhesion_properties"
sys.path.insert(0, str(ADHESION_PROPERTIES_DIR))
