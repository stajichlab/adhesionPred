"""The chunk merging of aga1_domain_modules.py (SYNTHETIC inputs)."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "aga1_domain_modules", ROOT / "analysis/sorting_hat_run/aga1_domain_modules.py"
)
dm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dm)


def chunk(tmp_path, name, **files):
    d = tmp_path / name
    d.mkdir()
    for fname, text in files.items():
        (d / fname).write_text(text)
    return d


def test_chunk_names_are_matched_per_species(tmp_path):
    for n in ("Sp_A_001", "Sp_A_002", "Sp_B_001"):
        (tmp_path / "tmhmm" / n).mkdir(parents=True)
    got = dm.chunks([tmp_path], "tmhmm", "Sp_A")
    assert [p.name for p in got] == ["Sp_A_001", "Sp_A_002"]
    with pytest.raises(dm.DomainError):
        dm.chunks([tmp_path], "tmhmm", "Sp_C")


def test_merge_tmhmm_keeps_one_header_and_refuses_a_different_one(tmp_path):
    a = chunk(tmp_path, "a", **{"tmhmm.tsv": "protein_id\tlen\nP1\t10\n"})
    b = chunk(tmp_path, "b", **{"tmhmm.tsv": "protein_id\tlen\nP2\t20\n"})
    out = tmp_path / "m.tsv"
    assert dm.merge_tmhmm([a, b], out) == 2
    assert out.read_text() == "protein_id\tlen\nP1\t10\nP2\t20\n"
    c = chunk(tmp_path, "c", **{"tmhmm.tsv": "id\tlen\nP3\t5\n"})
    with pytest.raises(dm.DomainError):
        dm.merge_tmhmm([a, c], out)


def test_merge_domtbl_needs_the_trailer_in_every_chunk(tmp_path):
    ok = "# Option settings: hmmsearch --cut_ga\nP1  -  10 F  PF1\n# [ok]\n"
    a = chunk(tmp_path, "a", **{"domtbl.txt": ok})
    b = chunk(tmp_path, "b", **{"domtbl.txt": ok.replace("P1", "P2")})
    out = tmp_path / "m.domtbl"
    dm.merge_domtbl([a, b], out)
    lines = out.read_text().splitlines()
    assert lines[-1] == "# [ok]" and [x.split()[0] for x in lines if not x.startswith("#")] == [
        "P1",
        "P2",
    ]
    bad = chunk(tmp_path, "bad", **{"domtbl.txt": "P3  -  1 F  PF1\n"})
    with pytest.raises(dm.DomainError):
        dm.merge_domtbl([a, bad], out)


aa = importlib.util.module_from_spec(
    importlib.util.spec_from_file_location(
        "aga1_allergen_antigen_modules",
        ROOT / "analysis/sorting_hat_run/aga1_allergen_antigen_modules.py",
    )
)
aa.__spec__.loader.exec_module(aa)


def test_uniprot_protein_map_uses_the_cimg_locus_name(tmp_path):
    table = tmp_path / "u.tsv"
    table.write_text(
        "Entry\tGene Names (ORF)\tGene Names\tRefSeq\n"
        "J3K365\tCIMG_09796\tCIMG_09796\tXP_1;\n"
        "A0A0D8JT44\tCIMG_00001\tCOQ3\tXP_2;\n"  # gene symbol in the names, locus name is still used
        "I7MN50\tCIMG_13598 CIMG_13775\tCIMG_13598 CIMG_13775\tXP_3;\n"  # several: the first
        "X1\t\tFOO\t\n"  # no locus name
    )
    out = tmp_path / "map.tsv"
    assert aa.uniprot_protein_map(table, out) == (3, 4)
    assert out.read_text() == (
        "protein_id\tgene_id\nJ3K365\tCIMG_09796\nA0A0D8JT44\tCIMG_00001\nI7MN50\tCIMG_13598\n"
    )
