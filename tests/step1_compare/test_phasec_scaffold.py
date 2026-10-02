"""Phase C package rules: allowed imports, no BASH_SOURCE, the class mapping (spec 3.1)."""

import ast
import sys

import labelmap
import paths
import pytest

PHASEC = paths.STEP1_DIR / "phasec"
ALLOWED_THIRD_PARTY = {"numpy", "scipy", "sklearn"}


def imported_names(source: str) -> list[str]:
    names = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names += [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names.append((node.module or "").split(".")[0])
    return names


def forbidden_imports(source: str, local: set[str]) -> list[str]:
    return [
        n
        for n in imported_names(source)
        if n not in sys.stdlib_module_names and n not in local and n not in ALLOWED_THIRD_PARTY
    ]


def _local() -> set[str]:
    return {p.stem for p in PHASEC.glob("*.py")} | {p.stem for p in paths.STEP1_DIR.glob("*.py")}


def test_phasec_imports_only_allowed():
    modules = sorted(PHASEC.glob("*.py"))
    assert len(modules) >= 2, modules  # evalio and labelmap from Task 1 on
    for path in modules:
        assert forbidden_imports(path.read_text(), _local()) == [], path.name


def test_phasec_import_check_catches_a_forbidden_import():
    path = PHASEC / "labelmap.py"
    copy = path.read_text() + "\nimport pandas\nfrom torch import nn\n"
    assert forbidden_imports(copy, _local()) == ["pandas", "torch"]


def test_no_phasec_file_uses_bash_source():
    files = sorted(PHASEC.rglob("*.py")) + sorted(PHASEC.rglob("*.sh"))
    assert files
    for f in files:
        assert "BASH_SOURCE" not in f.read_text(), f


CASES = [
    ("P-ext", "", "pos", "wall"),
    ("P-ext", "P-gpi", "pos", "wall"),
    ("P-ext", "PM-TM", "neg", "PM-TM"),
    ("P-ext", "pm-unresolved", "excluded", "pm-unresolved"),
]


@pytest.mark.parametrize("label,d8,cls,stratum", CASES)
def test_class_of_truth_table_p_ext(label, d8, cls, stratum):
    assert labelmap.class_of(label, d8) == cls
    assert labelmap.stratum_of(label, "wall", d8) == stratum


def test_class_of_truth_table():
    # every label x d8_class pair of spec 3.1
    expected = {
        ("P-ext", ""): "pos",
        ("P-ext", "P-gpi"): "pos",
        ("P-ext", "PM-TM"): "neg",
        ("P-ext", "pm-unresolved"): "excluded",
    }
    for label in labelmap.LABELS:
        for d8 in labelmap.D8_CLASSES:
            if label == "P-ext":
                want = expected[(label, d8)]
            else:
                want = {"N-int": "neg", "N-sec": "neg"}.get(label, "excluded")
            assert labelmap.class_of(label, d8) == want, (label, d8)
    assert labelmap.stratum_of("P-ext", "extracellular-only", "") == "extracellular-only"
    assert labelmap.stratum_of("N-sec", "", "") == "N-sec"
    assert labelmap.stratum_of("ambiguous", "", "") == "ambiguous"


def test_class_of_refuses_unknown_values():
    with pytest.raises(ValueError, match="unknown label"):
        labelmap.class_of("P-wall", "")
    with pytest.raises(ValueError, match="unknown d8_class"):
        labelmap.class_of("P-ext", "TM")
    with pytest.raises(ValueError, match="subset"):
        labelmap.stratum_of("P-ext", "", "")
