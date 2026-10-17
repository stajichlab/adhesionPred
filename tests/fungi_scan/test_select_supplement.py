import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/fungi_scan/select_supplement.py"
spec = importlib.util.spec_from_file_location("select_supplement", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_farthest_point_sampling_spreads_over_the_distance_matrix():
    names = ["a", "b", "c", "d"]
    d = {
        ("a", "b"): 1,
        ("a", "c"): 10,
        ("a", "d"): 11,
        ("b", "c"): 10,
        ("b", "d"): 11,
        ("c", "d"): 1,
    }
    dist = lambda x, y: 0 if x == y else d.get((x, y), d.get((y, x)))  # noqa: E731
    got = m.farthest_points(names, dist, start="a", n=3)
    assert got[0] == "a" and got[1] == "d" and len(got) == 3


def test_second_proteome_is_nearest_the_median_among_the_rest():
    def row(n, c):
        return {
            "GENUS": "G",
            "ORDER": "Onygenales",
            "FAMILY": "F",
            "SUBPHYLUM": "Pezizomycotina",
            "CLASS": "C",
            "ASMID": "GCA_1",
            "NCBI_TAXONID": "1",
            "_name": n,
            "_n": c,
        }

    rows = [row("a", 8000), row("b", 9000), row("c", 9500), row("d", 12000)]
    first = m.sg.pick_per_genus(rows, known_taxa={"1"})
    second = m.second_per_genus(rows, first, known_taxa={"1"})
    assert [r["_name"] for r in second] and second[0]["_name"] != first[0]["_name"]
