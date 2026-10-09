import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/evaluate.py"
spec = importlib.util.spec_from_file_location("evaluate", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

CLUSTERS = {"c1": {"a", "b"}, "c2": {"c"}, "c3": {"d", "e"}}
MISSED = {"b", "c", "d"}


def test_a_cluster_is_recovered_when_a_pfam_missed_member_is_called():
    assert m.recovered_clusters({"b"}, CLUSTERS, MISSED) == {"c1"}
    assert m.recovered_clusters({"a"}, CLUSTERS, MISSED) == set()  # a is not Pfam-missed
    assert m.recovered_clusters({"b", "c", "d", "e"}, CLUSTERS, MISSED) == {"c1", "c2", "c3"}


def test_unrecovered_count_u():
    assert m.unrecovered({"c1"}, ["c1", "c2", "c3"]) == ["c2", "c3"]


def test_relaxed_pass_rule_needs_three_of_six():
    assert m.relaxed_recall_pass({"a", "b", "c"}, 6) is True
    assert m.relaxed_recall_pass({"a", "b"}, 6) is False


def test_hmm_decision_branches():
    assert m.hmm_decision(u=0, recovered_by_hmm=set(), unrecovered_by_relaxed=[]) == "not testable"
    assert (
        m.hmm_decision(u=1, recovered_by_hmm={"c2"}, unrecovered_by_relaxed=["c2"])
        == "not testable"
    )
    assert (
        m.hmm_decision(u=2, recovered_by_hmm={"c2"}, unrecovered_by_relaxed=["c2", "c3"])
        == "recall rule met"
    )
    assert (
        m.hmm_decision(u=2, recovered_by_hmm=set(), unrecovered_by_relaxed=["c2", "c3"])
        == "recall rule not met"
    )
    assert (
        m.hmm_decision(u=3, recovered_by_hmm={"c2"}, unrecovered_by_relaxed=["c2", "c3", "c4"])
        == "recall rule not met"
    )  # ceil(3/2) = 2 needed, only 1 recovered
    assert (
        m.hmm_decision(
            u=3, recovered_by_hmm={"c2", "c3"}, unrecovered_by_relaxed=["c2", "c3", "c4"]
        )
        == "recall rule met"
    )
    assert (
        m.hmm_decision(u=3, recovered_by_hmm={"c9"}, unrecovered_by_relaxed=["c2", "c3", "c4"])
        == "recall rule not met"
    )


def test_wilson_interval_of_three_of_six():
    lo, hi = m.wilson(3, 6)
    assert round(lo, 2) == 0.19 and round(hi, 2) == 0.81


def test_call_requires_score_r0_and_cysteines():
    scores = {"x": 5.0, "y": 5.0, "z": 5.0, "w": 1.0}
    r0 = {"x": "called", "y": "not_called", "z": "called", "w": "called"}
    cys = {"x": 8, "y": 8, "z": 7, "w": 9}
    assert m.called_set(scores, 2.6, r0, cys) == {"x"}
