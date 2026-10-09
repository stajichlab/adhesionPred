import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/cost.py"
spec = importlib.util.spec_from_file_location("cost", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_classify_extra_calls():
    relaxed = {"a", "b", "c", "d", "e"}
    strict = {"a"}
    hsba = {"b"}
    labelled = {"c"}
    out = m.classify(relaxed, strict, hsba, labelled)
    assert out["strict_also"] == ["a"]
    assert out["hsba_overlap"] == ["b"]
    assert out["labelled_extra"] == ["c"]
    assert out["unlabelled_extra"] == ["d", "e"]


def test_strict_takes_precedence_over_hsba_and_labels():
    out = m.classify({"a"}, {"a"}, {"a"}, {"a"})
    assert (
        out["strict_also"] == ["a"] and out["hsba_overlap"] == [] and out["unlabelled_extra"] == []
    )


def test_rate_per_10000():
    assert m.per_10000(3, 10000) == 3.0
    assert round(m.per_10000(4, 8107), 2) == 4.93
    assert m.per_10000(0, 0) == 0.0
