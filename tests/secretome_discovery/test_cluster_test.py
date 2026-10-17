import gzip
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/secretome_discovery/cluster_test.py"
spec = importlib.util.spec_from_file_location("cluster_test", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_scan_clusters_collects_pool_genera_for_pfam_free_clusters_only(tmp_path):
    f = tmp_path / "c.tsv.gz"
    lines = [
        "AAA_1\tAAA_1",
        "AAA_1\tBBB_1",
        "AAA_1\tCCC_1",  # cluster 1: pool proteins in genera gA and gB; no Pfam
        "AAA_2\tAAA_2",
        "AAA_2\tBBB_2",  # cluster 2: BBB_2 has a Pfam hit -> excluded
        "ZZZ_1\tZZZ_1",
        "ZZZ_1\tAAA_3",  # cluster 3: only one chosen member, not in the pool
    ]
    with gzip.open(f, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    tag_genus = {"AAA": "gA", "BBB": "gB", "CCC": "gC"}
    pool = {"S": {"AAA_1", "BBB_1", "AAA_2"}, "U": {"AAA_1", "BBB_1"}}
    pfam = {"BBB_2"}
    res = m.scan_clusters(f, tag_genus, pool, pfam)
    assert set(res["S"]) == {"AAA_1"}
    assert res["S"]["AAA_1"] == {"gA": 1, "gB": 1}
    assert res["U"]["AAA_1"] == {"gA": 1, "gB": 1}


def test_cluster_family_test_finds_a_target_restricted_cluster():
    genera = {f"T{i}": "target" for i in range(10)} | {f"B{i}": "background" for i in range(40)}
    pool_size = {g: 400 for g in genera}
    clusters = {
        "C1": {f"T{i}": 1 for i in range(8)},  # in 8 of 10 target genera, no background
        "C2": {f"T{i}": 1 for i in range(2)},  # only 2 target genera
        "C3": {f"B{i}": 1 for i in range(20)},
    }
    out = m.test_clusters(clusters, genera, pool_size, min_genera=3)
    assert [r["family"] for r in out["enriched"]] == ["C1"]
    assert out["n_tested"] == 3
