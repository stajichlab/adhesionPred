import gzip
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/secretome_discovery/run_discovery.py"
spec = importlib.util.spec_from_file_location("run_discovery", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def write_pools(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("protein\tpool\tlength\tpct_cys\tpct_pro\tpct_st\tsequons_per_100aa\tgpi\tpfam\n")
        for p, pool, pf in rows:
            fh.write(f"{p}\t{pool}\t100\t1\t1\t1\t1\t0\t{pf}\n")


def test_genus_table_counts_pool_size_and_proteins_per_family(tmp_path):
    f = tmp_path / "x.pools.tsv.gz"
    write_pools(f, [("a", "U", "PF1;PF2"), ("b", "S", "PF1"), ("c", "U", "")])
    pool_s, fam_s = m.genus_counts(f, pool="S")
    pool_u, fam_u = m.genus_counts(f, pool="U")
    assert pool_s == 3 and fam_s == {"PF1": 2, "PF2": 1}
    assert pool_u == 2 and fam_u == {"PF1": 1, "PF2": 1}


def test_family_rows_include_genera_without_a_hit_with_zero_hits():
    genus = {
        "G1": ("target", 100, {"PF1": 2}),
        "G2": ("background", 200, {}),
        "G3": ("background", 150, {"PF1": 1}),
    }
    fams = m.family_rows(genus)
    rows = {r["genus"]: r for r in fams["PF1"]}
    assert rows["G2"]["hits"] == 0 and rows["G2"]["pool"] == 200
    assert rows["G1"]["group"] == "target" and rows["G3"]["hits"] == 1


def test_pruned_tree_keeps_only_the_chosen_tips():
    nwk = "(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);"
    out = m.prune_newick(nwk, {"A", "C", "E"})
    import io

    from Bio import Phylo

    tips = sorted(t.name for t in Phylo.read(io.StringIO(out), "newick").get_terminals())
    assert tips == ["A", "C", "E"]
