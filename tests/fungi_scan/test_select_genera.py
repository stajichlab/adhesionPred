import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/fungi_scan/select_genera.py"
spec = importlib.util.spec_from_file_location("select_genera", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def row(
    genus, name, n, order="Onygenales", family="F", asm="GCA_1", taxon="1", sub="Pezizomycotina"
):
    return {
        "GENUS": genus,
        "ORDER": order,
        "FAMILY": family,
        "SUBPHYLUM": sub,
        "CLASS": "Eurotiomycetes",
        "ASMID": asm,
        "NCBI_TAXONID": taxon,
        "_name": name,
        "_n": n,
    }


def test_the_proteome_nearest_the_genus_median_is_chosen_not_the_largest():
    rows = [row("G", "a", 8000), row("G", "b", 9000), row("G", "c", 19000)]
    out = m.pick_per_genus(rows, known_taxa={"1"})
    assert [r["_name"] for r in out] == ["b"]


def test_hybrids_are_excluded_by_name_before_the_median_is_taken():
    rows = [row("G", "a", 8000), row("G", "hyb", 9000), row("G", "c", 9500)]
    out = m.pick_per_genus(rows, known_taxa={"1"}, hybrids={"hyb"})
    assert [r["_name"] for r in out] == ["a"] or [r["_name"] for r in out] == ["c"]
    assert out[0]["_name"] != "hyb"


def test_a_refseq_assembly_is_preferred():
    rows = [row("G", "a", 9000, asm="GCA_1"), row("G", "b", 12000, asm="GCF_2")]
    assert m.pick_per_genus(rows, known_taxa={"1"})[0]["_name"] == "b"


def test_fragmentary_and_empty_genus_rows_are_dropped():
    rows = [row("G", "a", 3999), row("", "b", 9000), row("H", "c", 5000)]
    out = m.pick_per_genus(rows, known_taxa={"1"})
    assert [r["GENUS"] for r in out] == ["H"]


def test_a_merged_taxon_id_is_mapped_and_the_row_is_kept():
    rows = [row("Paracoccidioides", "p", 7000, taxon="502780")]
    assert (
        m.pick_per_genus(rows, known_taxa={"121759"}, merged={"502780": "121759"})[0][
            "NCBI_TAXONID"
        ]
        == "121759"
    )
    assert m.pick_per_genus(rows, known_taxa={"121759"}) == []


def test_groups_target_background_and_excluded():
    assert m.group_of(row("A", "x", 1, order="Onygenales")) == "onygenales"
    assert (
        m.group_of(row("A", "x", 1, order="Eurotiales", family="Aspergillaceae"))
        == "aspergillaceae"
    )
    assert m.group_of(row("A", "x", 1, order="Eurotiales", family="Trichocomaceae")) == "background"
    assert m.group_of(row("A", "x", 1, order="Hypocreales")) == "background"
    assert m.group_of(row("A", "x", 1, order="Saccharomycetales", sub="Saccharomycotina")) == "out"
