import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/measure.py"
spec = importlib.util.spec_from_file_location("measure", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_wilson_matches_reference_values():
    # reference values computed independently (Wilson score interval, z = 1.959964)
    lo, hi = m.wilson(5, 5)
    assert round(lo, 3) == 0.566
    assert round(hi, 3) == 1.0
    assert round(m.wilson(4, 4)[0], 3) == 0.510
    assert round(m.wilson(4, 5)[0], 3) == 0.376
    assert m.wilson(0, 0) == (0.0, 1.0)


def test_rule_6_3_cases():
    assert m.keep_rescue([]) == ("no evidence", 0, 0)
    assert m.keep_rescue(["hydrophobin"] * 4)[0] == "no evidence"  # below the 5-cluster floor
    assert m.keep_rescue(["hydrophobin"] * 5)[0] == "keep"
    assert m.keep_rescue(["hydrophobin"] * 4 + ["not_hydrophobin"])[0] == "no evidence"
    assert (
        m.keep_rescue(["hydrophobin"] * 5 + ["unresolved"] * 10)[0] == "keep"
    )  # unresolved are not counted


def test_confusion_counts_and_exclusions():
    labels = {"a": 1, "b": 1, "c": 0, "d": 0, "e": 0}
    called = {"a", "c"}
    r = m.confusion(labels, called)
    assert (r["tp"], r["fn"], r["fp"], r["tn"]) == (1, 1, 1, 2 + 0)
    assert round(r["sens"], 3) == 0.5 and round(r["spec"], 3) == round(2 / 3, 3)
