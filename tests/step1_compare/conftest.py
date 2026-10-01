"""Make analysis/step1_compare importable as plain modules for tests."""

import importlib.util
import sys
from pathlib import Path

import pytest

STEP1_DIR = Path(__file__).resolve().parents[2] / "analysis" / "step1_compare"
TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(STEP1_DIR))
sys.path.insert(0, str(STEP1_DIR / "jobs"))
sys.path.insert(0, str(TESTS_DIR))


def load_script(name: str):
    """Import a numbered script such as 01_extract_go_truth.py as a module."""
    spec = importlib.util.spec_from_file_location(f"step1_{name}", STEP1_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def fixtures_dir() -> Path:
    return TESTS_DIR / "fixtures"


@pytest.fixture
def mini_ontology(fixtures_dir):
    import go_obo

    return go_obo.parse_obo(fixtures_dir / "mini.obo")
