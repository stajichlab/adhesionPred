import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/ship_rule.py"
spec = importlib.util.spec_from_file_location("ship_rule", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_cluster_outcomes_mixed_cluster_is_unresolved():
    clusters = {"a": "c1", "b": "c1", "c": "c2", "d": "c3"}
    dec = {"a": "hydrophobin", "b": "not_hydrophobin", "c": "hydrophobin", "d": "unresolved"}
    out = m.cluster_outcomes(clusters, dec)
    assert out == {"c1": "unresolved", "c2": "hydrophobin", "c3": "unresolved"}


def test_precision_rule_branches():
    assert m.precision_rule({}) == ("no evidence", 0, 0, 0.0)
    four = {f"c{i}": "hydrophobin" for i in range(4)}
    assert m.precision_rule(four)[0] == "no evidence"  # below the floor of 5
    five = {f"c{i}": "hydrophobin" for i in range(5)}
    assert m.precision_rule(five)[0] == "pass"
    low = {
        **{f"c{i}": "hydrophobin" for i in range(2)},
        **{f"n{i}": "not_hydrophobin" for i in range(4)},
    }
    assert m.precision_rule(low)[0] == "fail"
    unresolved_many = {**five, **{f"u{i}": "unresolved" for i in range(30)}}
    assert m.precision_rule(unresolved_many)[0] == "pass"


def test_ship_rule_needs_every_condition():
    ok = {"recall": True, "cost": True, "hard_negatives": True, "precision": "pass"}
    assert m.ship(ok) is True
    for k, v in (
        ("recall", False),
        ("cost", False),
        ("hard_negatives", False),
        ("precision", "fail"),
        ("precision", "no evidence"),
    ):
        assert m.ship({**ok, k: v}) is False
