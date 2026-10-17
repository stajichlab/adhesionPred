import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/fungi_scan/select_proteomes.py"
spec = importlib.util.spec_from_file_location("select_proteomes", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def row(genus, species, strain, order="Onygenales", taxon="1", tag="T", n=9000):
    return {
        "GENUS": genus,
        "SPECIES": f"{genus} {species}",
        "STRAIN": strain,
        "ORDER": order,
        "NCBI_TAXONID": taxon,
        "LOCUSTAG": tag + strain,
        "_n": n,
    }


def test_one_proteome_per_species_up_to_the_limit():
    rows = [
        row("Trichophyton", "rubrum", "A"),
        row("Trichophyton", "rubrum", "B"),
        row("Trichophyton", "tonsurans", "C"),
        row("Trichophyton", "verrucosum", "D"),
    ]
    out = m.choose(rows, [("Trichophyton", None, 2)])
    assert [r["SPECIES"] for r in out] == ["Trichophyton rubrum", "Trichophyton tonsurans"]


def test_fragmentary_proteomes_are_skipped():
    rows = [
        row("Histoplasma", "capsulatum", "A", n=300),
        row("Histoplasma", "capsulatum", "B", n=9000),
    ]
    out = m.choose(rows, [("Histoplasma", None, 1)])
    assert out[0]["STRAIN"] == "B"


def test_named_species_are_chosen_when_given():
    rows = [
        row("Candida", "albicans", "A"),
        row("Candida", "tropicalis", "B"),
        row("Candida", "parapsilosis", "C"),
    ]
    out = m.choose(rows, [("Candida", ["albicans", "parapsilosis"], 2)])
    assert [r["SPECIES"] for r in out] == ["Candida albicans", "Candida parapsilosis"]


def test_a_taxon_that_is_not_in_the_taxonomy_is_dropped():
    rows = [
        row("Aspergillus", "fumigatus", "A", taxon="5"),
        row("Aspergillus", "nidulans", "B", taxon="999999999"),
    ]
    out = m.choose(rows, [("Aspergillus", None, 2)], known_taxa={"5"})
    assert [r["SPECIES"] for r in out] == ["Aspergillus fumigatus"]


def test_selection_is_deterministic():
    rows = [row("Aspergillus", s, s) for s in "dcba"]
    assert m.choose(rows, [("Aspergillus", None, 3)]) == m.choose(
        list(reversed(rows)), [("Aspergillus", None, 3)]
    )
