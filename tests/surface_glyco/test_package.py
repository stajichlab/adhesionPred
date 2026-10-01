"""Package identity: names, version fallback, entry points."""

import importlib
import sys
from pathlib import Path

import pytest

tomllib = pytest.importorskip("tomllib")  # Python >= 3.11

ROOT = Path(__file__).resolve().parents[2]


def test_new_package_imports_and_has_a_version_string():
    import surface_glyco

    assert isinstance(surface_glyco.__version__, str) and surface_glyco.__version__


def test_old_package_is_gone_from_this_checkout():
    assert not (ROOT / "src" / "adhesion_predict").exists()
    try:
        old = importlib.import_module("adhesion_predict")
    except ModuleNotFoundError:
        return
    # An editable install from another checkout may exist in this environment; ours must not.
    assert not Path(old.__file__).resolve().is_relative_to(ROOT)
    sys.modules.pop("adhesion_predict", None)


def test_pyproject_names_match_the_decision():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert cfg["project"]["name"] == "surface_glyco"
    assert cfg["project"]["scripts"] == {
        "surface_glyco_predict": "surface_glyco.scripts.predict:cli",
        "surface_glyco_train": "surface_glyco.scripts.train:cli",
        "surface_glyco_evaluate": "surface_glyco.scripts.evaluate:cli",
    }
    assert cfg["tool"]["setuptools"]["package-data"] == {"surface_glyco": ["models/*"]}
    assert cfg["project"]["requires-python"] == ">=3.11"
