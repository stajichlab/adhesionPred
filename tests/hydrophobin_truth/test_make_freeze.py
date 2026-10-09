import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/make_freeze.py"
spec = importlib.util.spec_from_file_location("make_freeze", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_allowed_score_files_are_tuning_only():
    tuning = ["Scer_S288C", "Calb_SC5314"]
    ok = [
        "Scer_S288C.default.tblout",
        "hn_tuning.nobias.tblout",
        "t2_positives.default.tblout",
        "fold_03.Calb_SC5314.default.tblout",
        "fold_03.hn_tuning.nobias.tblout",
    ]
    assert m.disallowed(ok, tuning) == []


def test_a_test_proteome_or_test_part_file_is_disallowed():
    tuning = ["Scer_S288C"]
    bad = [
        "Afum_Af293.default.tblout",
        "hn_test.default.tblout",
        "lp.default.tblout",
        "fold_03.Bbas_ARSEF2860.default.tblout",
    ]
    assert sorted(m.disallowed(bad, tuning)) == sorted(bad)


def test_sha256_of_files_is_listed_sorted(tmp_path):
    (tmp_path / "b.txt").write_text("b")
    (tmp_path / "a.txt").write_text("a")
    out = m.hashes(tmp_path, ["b.txt", "a.txt"])
    assert list(out) == ["a.txt", "b.txt"] and all(len(v) == 64 for v in out.values())
