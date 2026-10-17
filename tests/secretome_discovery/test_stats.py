import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/secretome_discovery/stats.py"
spec = importlib.util.spec_from_file_location("sd_stats", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def rows(t_present, t_n, b_present, b_n, pool=400, bpool=400):
    out = []
    for i in range(t_n):
        out.append(
            {"genus": f"T{i}", "group": "target", "pool": pool, "hits": 3 if i < t_present else 0}
        )
    for i in range(b_n):
        out.append(
            {
                "genus": f"B{i}",
                "group": "background",
                "pool": bpool,
                "hits": 2 if i < b_present else 0,
            }
        )
    return out


def test_a_planted_family_is_found_with_a_small_p():
    s = m.family_stats(rows(10, 10, 0, 40))
    assert s["t_present"] == 10 and s["b_present"] == 0
    assert s["fisher_p"] < 1e-8
    assert s["odds_ratio"] > 1


def test_fisher_matches_a_hand_computed_value():
    # 3 of 3 target, 0 of 3 background: p(two-sided) = 2 / C(6,3) = 0.1
    s = m.family_stats(rows(3, 3, 0, 3))
    assert s["fisher_p"] == pytest.approx(0.1)


def test_copies_per_1000_use_the_pool_size():
    s = m.family_stats(rows(2, 2, 2, 2, pool=500, bpool=1000))
    assert s["t_copies_median"] == pytest.approx(1000 * 3 / 500)
    assert s["b_copies_median"] == pytest.approx(1000 * 2 / 1000)


def test_a_family_in_two_target_genera_is_not_reported():
    fams = {"PF_A": rows(2, 10, 0, 40), "PF_B": rows(8, 10, 0, 40)}
    out = m.report(fams, min_genera=3)
    assert [r["family"] for r in out["enriched"]] == ["PF_B"]


def test_depletion_is_listed_apart_from_enrichment():
    fams = {"PF_C": rows(0, 10, 30, 40), "PF_D": rows(8, 10, 0, 40)}
    out = m.report(fams, min_genera=3)
    assert [r["family"] for r in out["depleted"]] == ["PF_C"]
    assert [r["family"] for r in out["enriched"]] == ["PF_D"]


def test_bh_is_monotone_and_not_below_p():
    q = m.bh([0.001, 0.04, 0.03, 0.5])
    assert all(a >= b for a, b in zip(q, [0.001, 0.04, 0.03, 0.5], strict=True))
    assert max(q) <= 1.0


def test_logistic_p_is_not_significant_when_presence_follows_pool_size():
    # larger pool -> present; groups differ only in pool size
    rs = []
    for i in range(20):
        rs.append(
            {"genus": f"T{i}", "group": "target", "pool": 200 + 10 * i, "hits": 1 if i >= 10 else 0}
        )
        rs.append(
            {
                "genus": f"B{i}",
                "group": "background",
                "pool": 200 + 10 * i,
                "hits": 1 if i >= 10 else 0,
            }
        )
    s = m.family_stats(rs)
    assert s["logit_p"] > 0.5


def test_fitch_counts_a_single_gain_and_the_losses_inside_the_clade():
    # ((A,B),(C,D)),((E,F),(G,H)): present in A,B,C only -> one gain, one loss (D)
    nwk = "(((A,B),(C,D)),((E,F),(G,H)));"
    g, lo = m.fitch_gains_losses(nwk, {"A", "B", "C"}, all_tips=list("ABCDEFGH"))
    assert (g, lo) == (1, 1)


def test_fitch_counts_two_gains_for_two_separate_clades():
    nwk = "(((A,B),(C,D)),((E,F),(G,H)));"
    g, lo = m.fitch_gains_losses(nwk, {"A", "B", "G", "H"}, all_tips=list("ABCDEFGH"))
    assert (g, lo) == (2, 0)


def test_gains_inside_and_outside_the_target_clade_are_counted_apart():
    nwk = "(((A,B),(C,D)),((E,F),(G,H)));"
    d = m.fitch_detail(nwk, {"A", "B", "G", "H"}, target={"A", "B", "C", "D"})
    assert d == {"gains_in_target": 1, "gains_elsewhere": 1, "losses": 0}


def test_mann_whitney_gets_its_own_bh_over_families_with_three_target_genera():
    # PF_E: presence equal in both groups, copy number much higher in the target
    def r(genus, group, hits, pool=400):
        return {"genus": genus, "group": group, "pool": pool, "hits": hits}

    rows_e = [r(f"T{i}", "target", 8) for i in range(10)] + [
        r(f"B{i}", "background", 1) for i in range(40)
    ]
    rows_f = (
        [r(f"T{i}", "target", 1) for i in range(2)]
        + [r(f"T{i}", "target", 0) for i in range(2, 10)]
        + [r(f"B{i}", "background", 1) for i in range(40)]
    )
    out = m.report({"PF_E": rows_e, "PF_F": rows_f}, min_genera=3)
    allrows = {s["family"]: s for s in out["all"]}
    assert allrows["PF_E"]["mwu_q"] < 1e-5
    assert allrows["PF_F"]["mwu_q"] != allrows["PF_F"]["mwu_q"]  # nan: fewer than 3 target genera
