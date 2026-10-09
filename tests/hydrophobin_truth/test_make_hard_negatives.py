import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/make_hard_negatives.py"
spec = importlib.util.spec_from_file_location("make_hard_negatives", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def row(acc, name="Some protein", seq="MKC" * 20):
    return {"accession": acc, "name": name, "sequence": seq, "length": len(seq)}


def test_small_secreted_rule():
    ok = row("a", seq="M" + "AC" * 10 + "A" * 40)  # 10 Cys, length 61
    assert m.small_secreted_ok(ok)
    assert not m.small_secreted_ok(row("b", seq="M" + "A" * 100))  # no Cys
    assert not m.small_secreted_ok(row("c", seq="MC" * 400))  # too long
    assert not m.small_secreted_ok(
        row("d", name="Class I hydrophobin x", seq="M" + "AC" * 10 + "A" * 40)
    )
    assert not m.small_secreted_ok(row("e", seq="MCC" + "A" * 20))  # too short or too few Cys


def test_hydrophobin_like_entries_are_removed_from_every_group():
    groups = {"CFEM": [row("a"), row("b")], "HsbA": [row("b"), row("c")]}
    kept, removed = m.remove_hydrophobin_like(groups, hits={"b"}, names={"c": "Rodlet protein"})
    assert [r["accession"] for r in kept["CFEM"]] == ["a"]
    assert kept["HsbA"] == []
    assert sorted(removed) == ["b", "c"]


def test_unreviewed_sample_is_seeded_and_capped():
    pool = [row(f"u{i}") for i in range(1000)]
    a = m.sample_unreviewed(pool, 300, seed=5)
    b = m.sample_unreviewed(pool, 300, seed=5)
    assert a == b and len(a) == 300
    assert len(m.sample_unreviewed(pool[:10], 300, seed=5)) == 10


def test_cluster_parts_keep_clusters_whole_and_split_roughly_in_half():
    cl = {f"p{i}": f"c{i // 2}" for i in range(40)}
    part = m.split_parts(cl, seed=3)
    for c in set(cl.values()):
        assert len({part[p] for p, x in cl.items() if x == c}) == 1
    n_test = sum(1 for v in part.values() if v == "test")
    assert 10 <= n_test <= 30
