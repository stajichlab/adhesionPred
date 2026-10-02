"""Make analysis/step1_compare importable as plain modules for tests."""

import importlib.util
import sys
from pathlib import Path

import pytest

STEP1_DIR = Path(__file__).resolve().parents[2] / "analysis" / "step1_compare"
TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(STEP1_DIR))
sys.path.insert(0, str(STEP1_DIR / "jobs"))
sys.path.insert(0, str(STEP1_DIR / "phasec"))
sys.path.insert(0, str(TESTS_DIR))


def load_script(name: str):
    """Import a numbered script such as 01_extract_go_truth.py as a module."""
    spec = importlib.util.spec_from_file_location(f"step1_{name}", STEP1_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_phasec(name: str):
    """Import a numbered Phase C script such as 08_build_eval_tables.py as a module."""
    path = STEP1_DIR / "phasec" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"phasec_{name}", path)
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


@pytest.fixture(scope="session")
def phasec_chain(tmp_path_factory):
    """Phase C fixture chain 08 -> 09 (stub MMseqs2) -> 10 -> 11, built once per session.

    Step 11 runs only when phasec/11_evaluate.py exists (it is written after 10). Tests that
    change a file must copy chain["work"] first (phasec_fixture.copy_work)."""
    pytest.importorskip("sklearn")
    import phasec_fixture

    upto = 11 if (STEP1_DIR / "phasec" / "11_evaluate.py").exists() else 10
    return phasec_fixture.run_chain(tmp_path_factory.mktemp("phasec_chain"), upto, load_phasec)
