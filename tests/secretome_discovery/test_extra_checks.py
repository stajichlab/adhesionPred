import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/secretome_discovery/extra_checks.py"
spec = importlib.util.spec_from_file_location("extra_checks", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_clade_table_counts_species_with_the_family_per_clade():
    species = {"s1": 1, "s2": 1, "s3": 2}
    fams = {"s1": {"PF1": 2}, "s2": {}, "s3": {"PF1": 1, "PF2": 1}}
    pools = {"s1": 100, "s2": 200, "s3": 100}
    rows = m.clade_table(species, fams, pools, ["PF1", "PF2"])
    by = {(r["family"], r["clade"]): r for r in rows}
    assert (by[("PF1", 1)]["n_species"], by[("PF1", 1)]["n_present"]) == (2, 1)
    assert (by[("PF1", 2)]["n_species"], by[("PF1", 2)]["n_present"]) == (1, 1)
    assert by[("PF2", 1)]["n_present"] == 0
    assert by[("PF1", 1)]["copies_per_1000_median"] == m.median([20.0, 0.0])


def test_presence_concordance_between_first_and_second_proteome():
    first = {"G1": {"PF1", "PF2"}, "G2": {"PF1"}}
    second = {"G1": {"PF1"}, "G2": {"PF1"}}
    c = m.concordance(first, second, ["PF1", "PF2"])
    assert c["PF1"] == (2, 2)  # same presence call in 2 of 2 genera
    assert c["PF2"] == (
        1,
        2,
    )  # present in the first only for G1 -> differs in G1, same (absent) in G2


def test_second_proteome_substitution_replaces_the_target_genus_data():
    data = {"G1": ("target", 100, {"PF1": 3}), "B1": ("background", 100, {})}
    second = {"G1": (80, {"PF1": 0})}
    out = m.substitute(data, second)
    assert out["G1"] == ("target", 80, {"PF1": 0}) and out["B1"] == data["B1"]
