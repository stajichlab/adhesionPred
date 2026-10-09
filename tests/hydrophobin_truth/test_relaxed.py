import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/relaxed.py"
spec = importlib.util.spec_from_file_location("relaxed", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

TBLOUT = """#                                                               --- full sequence ---
# target name        accession  query name           accession    E-value  score  bias
p1 - Hydrophobin PF01185.24 1e-5 25.5 0.1 x
p1 - Eas PF22354.3 1e-3 30.2 0.0 x
p2 - Hydrophobin PF01185.24 1e-2 8.6 0.0 x
# [ok]
"""


def test_parse_tblout_keeps_the_best_score_per_target():
    assert m.parse_tblout(TBLOUT.splitlines()) == {"p1": 30.2, "p2": 8.6}


def cand(score, pid="x", proteome="A"):
    return {"id": pid, "proteome": proteome, "score": score}


SIZES = {"A": 10000, "B": 10000}


def test_lowest_cutoff_is_the_floor_when_nothing_scores_above_it():
    assert m.lowest_cutoff([], SIZES, {}, 5, 0.05, floor=0.0) == 0.0


def test_lowest_cutoff_respects_the_per_proteome_limit_and_returns_the_lowest_included_score():
    # 7 candidates in proteome A at scores 1..7: limit 5 per 10,000 means at most 5 calls, so the two lowest are excluded
    cands = [cand(float(s), pid=f"p{s}") for s in range(1, 8)]
    c = m.lowest_cutoff(cands, SIZES, {}, 5, 0.05, floor=0.0)
    assert c == 3.0  # calls p3..p7 = 5 calls; the cutoff is the lowest included score


def test_lowest_cutoff_applies_the_limit_in_every_proteome():
    cands = [cand(float(s), pid=f"a{s}", proteome="A") for s in range(1, 4)] + [
        cand(float(s) + 10, pid=f"b{s}", proteome="B") for s in range(1, 8)
    ]
    c = m.lowest_cutoff(cands, SIZES, {}, 5, 0.05, floor=0.0)
    assert c == 13.0  # B needs its two lowest (11, 12) excluded; the lowest included score is 13.0


def test_hard_negative_rate_can_bind():
    # no proteome candidates; a hard-negative group of 20 proteins with scores 1..20 and a rate limit of 10%
    hn = {"CFEM": [float(s) for s in range(1, 21)]}
    c = m.lowest_cutoff([], SIZES, hn, 5, 0.10, floor=0.0)
    assert c == 19.0  # 2 of 20 proteins at or above 19.0 is the 10% limit


def test_hsba_group_is_ignored_by_the_rate_rule():
    hn = {"HsbA": [float(s) for s in range(1, 21)]}
    assert m.lowest_cutoff([], SIZES, hn, 5, 0.0, floor=0.0, ignore_groups={"HsbA"}) == 0.0


def test_no_positive_score_is_an_input():
    import inspect

    assert "positive" not in " ".join(inspect.signature(m.lowest_cutoff).parameters)


def test_option_rule_picks_the_option_with_more_recovered_clusters():
    recovered = {"default": {"c1", "c2"}, "nobias": {"c1", "c2", "c3"}}
    assert m.choose_option(recovered) == "nobias"
    assert (
        m.choose_option({"default": {"c1"}, "nobias": {"c1"}}) == "default"
    )  # tie: default filters


def test_conditioned_extra_calls_require_r0_cys_and_exclude_strict_and_labelled():
    items = {
        "ok": {"score": 9.0, "r0": "called", "n_cys": 8, "strict": False, "labelled": False},
        "no_sp": {"score": 9.0, "r0": "not_called", "n_cys": 8, "strict": False, "labelled": False},
        "few_cys": {"score": 9.0, "r0": "called", "n_cys": 7, "strict": False, "labelled": False},
        "strict": {"score": 9.0, "r0": "called", "n_cys": 8, "strict": True, "labelled": False},
        "labelled": {"score": 9.0, "r0": "called", "n_cys": 8, "strict": False, "labelled": True},
    }
    assert m.conditioned_extra(items) == ["ok"]
