import importlib.util
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/literature/lp_clusters.py"
)
spec = importlib.util.spec_from_file_location("lp_clusters", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_reserve_keeps_only_clean_lp_clusters():
    clusters = {"lp1": "c1", "lp2": "c1", "lp3": "c2", "lp4": "c3", "t1": "c2"}
    lp = {"lp1", "lp2", "lp3", "lp4"}
    t12 = {"t1"}
    unverified = {"lp4"}
    assert m.reserve(clusters, lp, t12, unverified) == {"c1"}


def test_reserve_is_empty_when_every_lp_cluster_has_a_t2_member():
    assert m.reserve({"lp1": "c1", "t": "c1"}, {"lp1"}, {"t"}, set()) == set()


def test_gap_report_finds_actual_gaps():
    seq = (
        "MKC"
        + "A" * 7
        + "CC"
        + "A" * 39
        + "C"
        + "A" * 21
        + "C"
        + "A" * 5
        + "CC"
        + "A" * 17
        + "C"
        + "GG"
    )
    assert m.actual_gaps(seq) == [7, 0, 39, 21, 5, 0, 17]
    assert m.stated_gaps("CN{7}CCN{39}CN{21}CN{5}CCN{17}C") == [7, 39, 21, 5, 17]


def test_stated_gaps_of_a_missing_pattern_is_none():
    assert m.stated_gaps("") is None
