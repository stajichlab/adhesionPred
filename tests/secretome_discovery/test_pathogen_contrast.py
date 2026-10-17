import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/secretome_discovery/pathogen_contrast.py"
spec = importlib.util.spec_from_file_location("pathogen_contrast", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_contrast_counts_strict_broad_and_others():
    presence = {"A": {"PF1"}, "B": {"PF1", "PF2"}, "C": {"PF2"}, "D": set()}
    groups = {
        "A": "mammal_pathogen",
        "B": "mammal_pathogen_mixed",
        "C": "other_nonpathogen",
        "D": "unclassified",
    }
    rows = {r["family"]: r for r in m.contrast(presence, groups, ["PF1", "PF2"])}
    r1 = rows["PF1"]
    assert (r1["strict_present"], r1["strict_n"]) == (1, 1)
    assert (r1["others_present"], r1["others_n"]) == (
        1,
        3,
    )  # B counts among the others in the strict contrast
    assert (r1["broad_present"], r1["broad_n"]) == (2, 2)
    assert (r1["others_broad_present"], r1["others_broad_n"]) == (0, 2)


def test_the_shipped_group_file_covers_every_onygenales_genus():
    import csv

    root = Path(__file__).resolve().parents[2]
    have = {
        r["genus"]
        for r in csv.DictReader(
            open(root / "analysis/secretome_discovery/onygenales_host_groups.tsv"), delimiter="\t"
        )
    }
    sel = {
        r["genus"]
        for r in csv.DictReader(
            open(root / "analysis/fungi_scan/selection_discovery.tsv"), delimiter="\t"
        )
        if r["group"] == "onygenales"
    }
    assert sel <= have
