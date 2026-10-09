import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/evidence_sheets.py"
spec = importlib.util.spec_from_file_location("evidence_sheets", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_cys_gaps():
    assert m.cys_gaps("MKCAACCAAAC") == [2, 0, 3]
    assert m.cys_gaps("MKAAA") == []


def test_parse_blast_keeps_best_hits_in_order():
    text = "q1\tA1\t80.0\t90\t1e-30\t120\tprot A\nq1\tA2\t60.0\t80\t1e-10\t60\tprot B\nq2\tA3\t99.0\t100\t1e-50\t200\tprot C\n"
    out = m.parse_blast(text.splitlines(), top=1)
    assert out["q1"] == ["A1 80.0% cov90 1e-30 prot A"] and out["q2"] == [
        "A3 99.0% cov100 1e-50 prot C"
    ]


def test_kind_labels():
    assert m.kinds(strict=True, relaxed=True, labelled=False) == "strict_and_relaxed_unlabelled"
    assert m.kinds(strict=False, relaxed=True, labelled=False) == "relaxed_only_unlabelled"
    assert m.kinds(strict=True, relaxed=False, labelled=False) == "strict_only_unlabelled"
    assert m.kinds(strict=True, relaxed=True, labelled=True) is None
