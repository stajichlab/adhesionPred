import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/prefreeze_data.py"
spec = importlib.util.spec_from_file_location("prefreeze_data", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_candidates_from_scores_apply_the_conditions():
    scores = {"a": 12.0, "b": 9.0, "c": 7.0, "d": 6.0, "e": 5.0}
    meta = {
        "a": {"r0": "called", "n_cys": 8, "strict": False, "labelled": False},
        "b": {"r0": "not_called", "n_cys": 8, "strict": False, "labelled": False},
        "c": {"r0": "called", "n_cys": 6, "strict": False, "labelled": False},
        "d": {"r0": "called", "n_cys": 9, "strict": True, "labelled": False},
        "e": {"r0": "called", "n_cys": 9, "strict": False, "labelled": True},
    }
    cands = m.candidates("P", scores, meta)
    assert cands == [{"id": "a", "proteome": "P", "score": 12.0}]


def test_hard_negative_scores_use_a_floor_value_for_proteins_that_are_not_called_by_the_conditions():
    scores = {"x-G": 20.0, "y-G": 30.0}
    meta = {
        "x-G": {"r0": "called", "n_cys": 8},
        "y-G": {"r0": "not_called", "n_cys": 8},
        "z-G": {"r0": "called", "n_cys": 8},  # no hmmsearch hit
    }
    groups = {"x-G": "G", "y-G": "G", "z-G": "G"}
    out = m.hard_negative_scores(scores, meta, groups)
    assert sorted(out["G"]) == sorted([20.0, m.NOT_CALLED, m.NOT_CALLED])


def test_n_cys_counts_cysteines():
    assert m.n_cys("MKCACC") == 3


def test_exact_sequence_labels():
    proteome = {"p1": "MKAC", "p2": "MKAA"}
    assert m.exact_matches({"lp1": "MKAC"}, proteome) == {"p1"}
