import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/apply_rubric.py"
spec = importlib.util.spec_from_file_location("apply_rubric", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def d(A="no", B="no", C="no", D="no", E="no"):
    return m.decide({"line_A": A, "line_B": B, "line_C": C, "line_D": D, "line_E": E})


def test_c_alone_is_hydrophobin():
    assert d(C="yes")[0] == "hydrophobin"


def test_a_needs_b_or_d():
    assert d(A="yes")[0] == "unresolved"
    assert d(A="yes", B="yes")[0] == "hydrophobin"
    assert d(A="yes", D="yes")[0] == "hydrophobin"
    assert d(B="yes")[0] == "unresolved"  # only the Pfam hit
    assert d(D="yes")[0] == "unresolved"


def test_e_alone_is_not_hydrophobin():
    assert d(E="yes")[0] == "not_hydrophobin"


def test_conflict_between_hydrophobin_evidence_and_e_is_unresolved():
    out = d(A="yes", B="yes", E="yes")
    assert out[0] == "unresolved" and "conflict" in out[1]


def test_c_overrides_e():
    assert d(C="yes", E="yes")[0] == "hydrophobin"


def test_nothing_is_unresolved():
    assert d()[0] == "unresolved"


def test_not_checked_does_not_count_as_yes():
    assert d(A="yes", B="not checked")[0] == "unresolved"
