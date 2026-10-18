"""Tests for analysis/wall_fusion_truth/map_kiss.py (no network: the UniProt lookup is a stub)."""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/wall_fusion_truth/map_kiss.py"
spec = importlib.util.spec_from_file_location("map_kiss", SCRIPT)
mk = importlib.util.module_from_spec(spec)
sys.modules["map_kiss"] = mk
spec.loader.exec_module(mk)

SEQ = "M" + "ACDEFGHIKL" * 6


def entry(acc, genes, reviewed=True, seq=SEQ, organism_id=330879):
    return {
        "accession": acc,
        "reviewed": reviewed,
        "entry_name": acc + "_X",
        "genes": genes,
        "organism": "x",
        "organism_id": organism_id,
        "length": len(seq),
        "sequence": seq,
    }


# gene_exact term -> hits (the stub ignores the taxon)
FAKE = {
    "AFUA_3G12700": [entry("Q4WRA1", "AFUA_3G12700")],
    "hwp1": [entry("P46593", "HWP1", organism_id=237561)],
}


def lookup(kind, term, taxon):
    assert kind == "gene_exact"
    return FAKE.get(term, [])


def row(gene, jgi, cluster, category="adhesion", afu=""):
    return {
        "gene_or_family_id": gene,
        "protein_id_jgi": jgi,
        "cluster_number": cluster,
        "category": category,
        "ortholog_Afu": afu,
    }


FIXTURE = [
    row("glfB_Afu3g12700", "4407_Aspfu1", "417"),  # locus tag in the row: mapped
    row("hwp1", "55000_Canalb1", "5266"),  # gene name only: candidate, not mapped
    row("unknownX", "", "9999", afu="NA"),  # no JGI ID and no Afu locus: lost
]


def mapped_fixture():
    out = []
    for r in FIXTURE:
        m = mk.map_row(r, "330879", lookup)
        r = dict(r, map_status=m["status"], map_reason=m["reason"], accession=m["accession"])
        out.append(r)
    return out


def test_locus_tag_row_is_mapped_and_gene_name_row_is_only_a_candidate():
    rows = mapped_fixture()
    assert rows[0]["map_status"] == "mapped_locus_tag"
    assert rows[0]["accession"] == "Q4WRA1"
    assert rows[1]["map_status"] == "candidate_gene_name"
    assert rows[2]["map_status"] == "lost"
    assert rows[2]["map_reason"] == "no_jgi_id_and_no_afu_locus"


def test_loss_count():
    c = mk.count_rows(mapped_fixture())
    assert c["rows"] == 3
    assert c["with_jgi_id"] == 2
    assert c["without_jgi_id"] == 1
    assert c["mapped_locus_tag"] == 1
    assert c["candidate_gene_name"] == 1
    assert c["lost"] == 1
    assert c["lost_reasons"] == {"no_jgi_id_and_no_afu_locus": 1}
    # rows are not lost by being candidates, and the three states add up to the rows
    assert (
        c["mapped_locus_tag"] + c["mapped_family_rep"] + c["candidate_gene_name"] + c["lost"] == 3
    )


def test_family_row_without_id_uses_afu_locus():
    r = row("glfB/Afu3g12700", "", "417", afu="glfB/Afu3g12700")
    m = mk.map_row(r, "330879", lookup)
    assert m["status"] == "mapped_family_rep"
    assert m["match_term"] == "AFUA_3G12700"


def test_ambiguous_and_implausible_hits_are_not_mapped():
    def amb(kind, term, taxon):
        return [entry("A1", "g", reviewed=False), entry("A2", "g", reviewed=False)]

    m = mk.map_row(row("hwp1", "1_Canalb1", "1"), "1", amb)
    assert m["status"] == "lost" and m["reason"] == "gene_name_ambiguous"

    def xrun(kind, term, taxon):
        return [entry("A3", "g", seq="M" + "X" * 30 + "ACDEFGHIKL")]

    m = mk.map_row(row("glfB_Afu3g12700", "4407_Aspfu1", "1"), "1", xrun)
    assert m["status"] == "lost" and m["reason"].startswith("implausible_sequence")


def test_reviewed_entry_is_preferred_over_unreviewed():
    hit, why = mk.pick_entry([entry("U1", "g", reviewed=False), entry("R1", "g", reviewed=True)])
    assert hit["accession"] == "R1" and why == ""


def test_organism_code_slice_and_locus_tags():
    assert mk.organism_code("4407_Aspfu1") == "Aspfu1"
    assert mk.organism_code("") == ""
    assert mk.slice_of("pheromone receptor (mating, fusion-adjacent)") == "fusion"
    assert mk.slice_of("cell wall (GPI anchor assembly)") == "wall"
    assert mk.locus_tags("glfB_Afu3g12700") == ["AFUA_3G12700"]
    assert mk.locus_tags("BATDEDRAFT_20500") == ["BATDEDRAFT_20500"]
    assert mk.locus_tags("hwp1") == []


def test_clade_floor_uses_clusters_not_rows():
    rows = []
    for i, cl in enumerate(["1", "1", "2", "3", "4", "5"]):
        rows.append(
            {
                "slice": "wall",
                "cluster_number": cl,
                "clade_coarse": "Pezizomycotina" if i < 5 else "Basidiomycota",
                "clade_fine": "Pezizomycotina",
                "organism_code": "x",
            }
        )
    tabs = mk.clade_tables(rows, lambda r: True)
    pez = tabs[("wall", "coarse", "Pezizomycotina")]
    assert pez["rows"] == 5 and pez["clusters"] == 4
    assert [c for c, _, _ in mk.eligible(tabs, "wall", "coarse", 3)] == ["Pezizomycotina"]
    assert mk.eligible(tabs, "wall", "coarse", 5) == []


def test_lineage_to_clades():
    lin = {
        "phylum": ("Ascomycota", "4890"),
        "subphylum": ("Pezizomycotina", "147538"),
        "class": ("Eurotiomycetes", "147545"),
    }
    assert mk.clades_from_lineage(lin) == ("Pezizomycotina", "Eurotiomycetes")
    bas = {"phylum": ("Basidiomycota", "5204"), "class": ("Agaricomycetes", "155619")}
    assert mk.clades_from_lineage(bas) == ("Basidiomycota", "Agaricomycetes")
