"""Tests of the metric functions on toy data."""

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis" / "lectin_validation"))

import lectin_metrics as lm  # noqa: E402


def test_wilson_zero_successes():
    lo, hi = lm.wilson(0, 10)
    assert lo == 0.0
    assert hi == pytest.approx(0.2775, abs=1e-3)


def test_wilson_all_successes():
    lo, hi = lm.wilson(10, 10)
    assert hi == pytest.approx(1.0, abs=1e-9)
    assert lo == pytest.approx(0.7225, abs=1e-3)


def test_wilson_half():
    lo, hi = lm.wilson(5, 10)
    assert lo == pytest.approx(0.2366, abs=1e-3)
    assert hi == pytest.approx(0.7634, abs=1e-3)


def test_wilson_empty_is_nan():
    lo, hi = lm.wilson(0, 0)
    assert math.isnan(lo)
    assert math.isnan(hi)


def test_wilson_rejects_bad_k():
    with pytest.raises(ValueError):
        lm.wilson(3, 2)


def test_hit_accessions_models_and_floors():
    rows = [("a", "PF1", 30.0), ("b", "PF1", 120.0), ("c", "PF2", 50.0), ("d", "PF3", 10.0)]
    assert lm.hit_accessions(rows, ["PF1"]) == {"a", "b"}
    assert lm.hit_accessions(rows, ["PF1"], {"PF1": 100.0}) == {"b"}
    assert lm.hit_accessions(rows, ["PF1", "PF2"], {"PF1": 100.0}) == {"b", "c"}
    assert lm.hit_accessions(rows, []) == set()


def test_summarize_cluster_hit_if_any_member_hits():
    members = [("a", "c1"), ("b", "c1"), ("c", "c2"), ("d", "c3")]
    s = lm.summarize(members, {"b", "d"})
    assert (s["n_seq"], s["k_seq"]) == (4, 2)
    assert (s["n_clu"], s["k_clu"]) == (3, 2)


def test_summarize_empty_group():
    s = lm.summarize([], set())
    assert s["n_seq"] == 0
    assert s["n_clu"] == 0
    assert math.isnan(s["clu_lo"])


def test_summarize_no_hits():
    s = lm.summarize([("a", "c1"), ("b", "c2")], set())
    assert s["k_clu"] == 0
    assert s["clu_lo"] == 0.0


def test_fmt():
    assert lm.fmt(0, 0, float("nan"), float("nan")) == "n=0"
    assert lm.fmt(3, 4, 0.3, 0.95) == "3/4 (0.30-0.95)"


def test_meets_floor():
    assert lm.meets_floor(20, 20)
    assert not lm.meets_floor(19, 40)
    assert not lm.meets_floor(40, 3)
