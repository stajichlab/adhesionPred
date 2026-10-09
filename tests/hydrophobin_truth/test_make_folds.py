import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/make_folds.py"
spec = importlib.util.spec_from_file_location("make_folds", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_one_fold_per_cluster_and_pfam_missed_flag():
    clusters = {"a": "c1", "b": "c1", "c": "c2", "d": "c3"}
    folds = m.loco_folds(clusters, pfam_missed={"c", "b"})
    assert [f["held_out_cluster"] for f in folds] == ["c1", "c2", "c3"]
    assert [f["pfam_missed"] for f in folds] == [1, 1, 0]
    assert folds[0]["train"] == ["c", "d"] and folds[0]["test"] == ["a", "b"]


GROUPS = {
    "Afum": ["Afum_Af293", "Afum_A1163", "Afum_W72310"],
    "Bbas": ["Bbas"],
    "Fful": ["Fful"],
    "Fgra": ["Fgra"],
    "Pexp": ["Pexp"],
    "Post": ["Post"],
    "Scer": ["Scer"],
    "Calb": ["Calb"],
    "Cimm": ["Cimm"],
    "Bder": ["Bder"],
}
TRUTH = {"Afum", "Bbas", "Fful", "Fgra", "Pexp", "Post"}


def test_species_split_is_reproducible_balanced_and_keeps_strains_together():
    a = m.split_species_groups(GROUPS, seed=1, truth_groups=TRUTH)
    b = m.split_species_groups(GROUPS, seed=1, truth_groups=TRUTH)
    assert a == b
    assert sorted(a.values()).count("tuning") == 5 and sorted(a.values()).count("test") == 5
    assert len({a[g] for g in ["Afum"]}) == 1
    for part in ("tuning", "test"):
        assert sum(1 for g, p in a.items() if p == part and g in TRUTH) >= 2


def test_every_seed_satisfies_the_truth_group_constraint():
    for seed in range(50):
        s = m.split_species_groups(GROUPS, seed=seed, truth_groups=TRUTH)
        for part in ("tuning", "test"):
            assert sum(1 for g, p in s.items() if p == part and g in TRUTH) >= 2
