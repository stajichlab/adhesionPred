import importlib.util
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/make_calibration_inputs.py"
)
spec = importlib.util.spec_from_file_location("make_calibration_inputs", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_rows_label_positives_negatives_and_exclude_t3():
    proteome = ["p1", "p2", "p3", "p4", "p5"]
    rows = m.build_rows(proteome, positives={"p1", "p2"}, excluded={"p3"})
    assert rows == [("p1", 1), ("p2", 1), ("p4", 0), ("p5", 0)]


def test_positive_missing_from_proteome_is_an_error():
    try:
        m.build_rows(["p1"], positives={"p9"}, excluded=set())
    except ValueError as e:
        assert "p9" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_every_row_gets_a_cluster_or_the_table_is_refused():
    rows = [("p1", 1), ("p2", 0)]
    out = m.attach_clusters(rows, {"p1": "c1", "p2": "c2"})
    assert out == [("p1", 1, "c1"), ("p2", 0, "c2")]
    try:
        m.attach_clusters(rows, {"p1": "c1"})
    except ValueError as e:
        assert "p2" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_split_part_of_a_cluster_is_one_value():
    clusters = {"p1": "c1", "p2": "c1", "p3": "c2", "p4": "c3", "p5": "c3"}
    part = m.split_parts(clusters, seed=3, dev_fraction=0.34)
    assert part["p1"] == part["p2"] and part["p4"] == part["p5"]
