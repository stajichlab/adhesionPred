"""clade_floors.py: clusters per clade and stratum, floor status (fixture, no data)."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "clade_floors", ROOT / "analysis/wall_fusion_truth/clade_floors.py"
)
cf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cf)


def cl(i, classes, sl, fine, coarse=""):
    return {
        "cluster_id": f"KS{i}",
        "classes": classes,
        "slices": sl,
        "clades_fine": fine,
        "clades_coarse": coarse or fine,
    }


def table(rows):
    return {(r[0], r[1], r[2]): r[3:] for r in cf.clade_floors(rows)}


def test_counts_clade_only_and_shared_and_status_at_the_floor():
    rows = [cl(i, "positive", "wall", "A") for i in range(4)]  # 4 clade-only in A
    rows.append(cl(10, "positive", "wall", "A,B"))  # shared A and B
    t = table(rows)
    assert t[("fine", "wall", "A")] == [5, 4, 1, "estimated", "descriptive", "descriptive"]
    assert t[("fine", "wall", "B")] == [1, 0, 1, "descriptive", "none", "descriptive"]


def test_negative_only_clusters_and_the_negative_member_do_not_count():
    rows = [cl(1, "negative", "", "A"), cl(2, "negative,positive", "adhesion", "A")]
    t = table(rows)
    assert t[("fine", "adhesion", "A")][0] == 1
    assert ("fine", "", "A") not in t


def test_fusion_is_always_descriptive_even_with_many_clusters():
    rows = [cl(i, "positive", "fusion", "A") for i in range(9)]
    t = table(rows)
    assert t[("fine", "fusion", "A")][:6] == [9, 9, 0, "descriptive", "descriptive", "none"]


def test_levels_are_separate():
    rows = [cl(1, "positive", "adhesion", "A", "X"), cl(2, "positive", "adhesion", "B", "X")]
    t = table(rows)
    assert t[("coarse", "adhesion", "X")][:3] == [2, 2, 0]
    assert t[("fine", "adhesion", "A")][:3] == [1, 1, 0]


def test_committed_table_is_what_the_script_makes():
    rec = ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22"
    import csv

    with open(rec / "clusters_split.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    want = (
        "\t".join(cf.HEADER)
        + "\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in cf.clade_floors(rows))
    )
    assert (rec / "clade_floors.tsv").read_text() == want
