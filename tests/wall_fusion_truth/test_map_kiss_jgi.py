"""Tests of the JGI-file step of map_kiss.py (fixture FASTA files only; no JGI data, no network)."""

import csv
import gzip
import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/wall_fusion_truth/map_kiss.py"
spec = importlib.util.spec_from_file_location("map_kiss_jgi", SCRIPT)
mk = importlib.util.module_from_spec(spec)
sys.modules["map_kiss_jgi"] = mk
spec.loader.exec_module(mk)

S1 = "M" + "ACDEFGHIKL" * 6
S2 = "M" + "LKIHGFEDCA" * 6
S3 = "M" + "ACDEFGHIKM" * 6  # similar to S1 (one change per ten)
S4 = "W" * 61


def write_gz(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for head, seq in records:
            fh.write(f">{head}\n{seq[:30]}\n{seq[30:]}*\n")


def make_dir(tmp_path: Path, portals: dict[str, list[tuple[str, str]]]) -> Path:
    d = tmp_path / "jgi"
    rows = []
    for portal, recs in portals.items():
        f = d / portal / "Mycocosm" / "Proteins" / f"{portal}.fasta.gz"
        write_gz(f, recs)
        rows.append(
            [portal, "x", f.name, "use", len(recs), hashlib.sha256(f.read_bytes()).hexdigest()]
        )
    with (d / "FILE_TABLE.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["portal", "organism", "file", "role", "n_proteins", "sha256"])
        w.writerows(rows)
    return d


def row(gene, jgi, status="lost", seq="", genes="", slice_="wall", cluster="1", code=None):
    return {
        "gene_or_family_id": gene,
        "protein_id_jgi": jgi,
        "organism_code": code if code is not None else mk.organism_code(jgi),
        "slice": slice_,
        "cluster_number": cluster,
        "map_status": status,
        "map_reason": ""
        if status != "lost"
        else "jgi_id_has_no_public_sequence_route; gene_name_no_hit",
        "sequence": seq,
        "genes": genes,
    }


def test_portal_rename_and_protein_id():
    assert mk.jgi_portal("Ustma22") == "Ustma2_2"
    assert mk.jgi_portal("CryneH991") == "Cryne_H99_1"
    assert mk.jgi_portal("Ashgo11") == "Ashgo1_1"
    assert mk.jgi_portal("Botci1") == "Bcine1"
    assert mk.jgi_portal("Copci") == "Copci_AmutBmut1"
    assert mk.jgi_portal("Canalb1") == "Canalb1"
    assert mk.jgi_protein_id("12141_Ustma22") == "12141"
    assert mk.jgi_protein_id("") == ""


def test_header_variants_cover_the_portal_naming_styles():
    assert "CAALFM_CR00610WA" in mk.header_variants("CAALFM_CR00610WAT0")
    assert "NCU06242" in mk.header_variants("NCU06242T0")
    assert "BCIN_02G00390" in mk.header_variants("BCIN_02g00390.1 proteinId=6405 transcriptId=6405")
    assert "UMAG_11562" in mk.header_variants("mRNA_UMAG_11562")


def test_fasta_reader_takes_only_wanted_ids_of_the_right_portal(tmp_path):
    d = make_dir(
        tmp_path,
        {
            "Schpo1": [
                ("jgi|Schpo1|3290|bgs1", S1),
                ("jgi|Schpo1|5055|bgs4", S2),
                ("jgi|OtherPortal|7|zzz", S3),
                ("not a jgi header", S3),
            ]
        },
    )
    f = next((d / "Schpo1").rglob("*.fasta.gz"))
    got = mk.read_jgi_fasta(f, "Schpo1", {"3290", "7"})
    assert got == {
        "3290": ("bgs1", S1)
    }  # trailing '*' removed; other portal and bad header skipped


def test_sha256_mismatch_and_duplicate_use_are_refused(tmp_path):
    d = make_dir(tmp_path, {"Schpo1": [("jgi|Schpo1|1|x", S1)]})
    table = mk.read_jgi_table(d / "FILE_TABLE.tsv")
    assert mk.load_jgi_proteins(d, table, {"Schpo1": {"1"}})["Schpo1"]["1"][1] == S1
    table["Schpo1"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="sha256"):
        mk.load_jgi_proteins(d, table, {"Schpo1": {"1"}})
    lines = (d / "FILE_TABLE.tsv").read_text().splitlines()
    (d / "FILE_TABLE.tsv").write_text("\n".join([*lines, lines[1]]) + "\n")
    with pytest.raises(ValueError, match="two files"):
        mk.read_jgi_table(d / "FILE_TABLE.tsv")


def test_gene_name_equal_to_header_name_is_mapped():
    r = row("bgs1", "3290_Schpo1")
    prot = {"Schpo1": {"3290": ("bgs1", S1)}}
    mk.apply_jgi([r], prot, {"Schpo1"})
    assert r["map_status"] == "mapped_jgi_file" and r["jgi_accept"] == "row_verified"
    assert r["jgi_check_how"] == "row_gene_name_equals_header_name"
    assert (
        r["jgi_seq_sha256"] == hashlib.sha256(S1.encode()).hexdigest() and r["jgi_length"] == "61"
    )


def test_candidate_agrees_by_identical_sequence_or_gene_names_and_differs_otherwise():
    same = row("hyr1", "1_Canalb1", "candidate_gene_name", seq=S1, genes="HYR1")
    byname = row("iff4", "2_Canalb1", "candidate_gene_name", seq=S2, genes="IFF4 CAALFM_C1WA")
    diff = row("als9", "3_Canalb1", "candidate_gene_name", seq=S2, genes="ALS9")
    prot = {
        "Canalb1": {
            "1": ("CAALFM_X", S1),
            "2": ("CAALFM_C1WAT0", S1),
            "3": ("CAALFM_Y", S4),
        }
    }
    mk.apply_jgi([same, byname, diff], prot, {"Canalb1"})
    assert (
        same["jgi_check_how"] == "uniprot_sequence_identical"
        and same["map_status"] == "mapped_jgi_file"
    )
    assert byname["jgi_check_how"] == "uniprot_gene_names_carry_header_locus"
    assert byname["map_status"] == "mapped_jgi_file"
    # a different sequence with no locus agreement is a mismatch: the row is not mapped
    assert diff["jgi_check"] == "differs" and diff["map_status"] == "lost"
    assert diff["map_reason"] == "jgi_file_differs_from_uniprot_candidate"


def test_similar_sequence_is_unverifiable_and_candidate_stays_a_candidate():
    r = row("coh2", "9_Copci", "candidate_gene_name", seq=S1, genes="COH2")
    prot = {"Copci_AmutBmut1": {"9": ("fgenesh1_pm.31_#_21", S3)}}
    mk.apply_jgi([r], prot, {"Copci_AmutBmut1"})
    assert r["jgi_vs_uniprot"] == "similar" and r["jgi_check"] == "unverifiable"
    assert r["map_status"] == "candidate_gene_name" and r["jgi_accept"] == ""


def test_locus_tag_of_the_row_must_match_the_header_scheme():
    ok = row("NCU08192", "1_Neucr2", "mapped_locus_tag", seq=S1)
    bad = row("NCU08193", "2_Neucr2", "mapped_locus_tag", seq=S1)
    other = row("CND01240", "3_CryneH991")
    num = row("UM05433", "4_Ustma22")
    prot = {
        "Neucr2": {"1": ("NCU08192T0", S1), "2": ("NCU99999T0", S2)},
        "Cryne_H99_1": {"3": ("CNAG_00996T0", S1)},
        "Ustma2_2": {"4": ("mRNA_UMAG_05433", S1)},
    }
    mk.apply_jgi([ok, bad, other, num], prot, {"Neucr2", "Cryne_H99_1", "Ustma2_2"})
    assert ok["jgi_check"] == "agree" and ok["map_status"] == "mapped_locus_tag"
    # a locus-tag row keeps its status even when the file disagrees; the check records it
    assert bad["jgi_check"] == "differs" and bad["map_status"] == "mapped_locus_tag"
    assert other["jgi_check"] == "unverifiable"
    assert other["jgi_check_how"] == "row_locus_tag_other_naming_scheme"
    assert num["jgi_check"] == "agree" and num["jgi_check_how"] == "row_locus_tag_same_number"


def test_absent_id_and_unprovided_portal_are_not_mapped():
    a = row("x", "5_Canalb1")
    b = row("y", "6_Aspfu1")
    mk.apply_jgi([a, b], {"Canalb1": {}}, {"Canalb1"})
    assert a["jgi_status"] == "absent" and a["map_status"] == "lost"
    assert a["map_reason"].endswith("jgi_id_absent_from_file")
    assert b["jgi_status"] == "portal_not_provided"
    assert b["map_reason"].endswith("jgi_portal_not_provided")


def test_rows_without_name_evidence_need_a_validated_portal():
    agree = [row(f"NCU0000{i}", f"{i}_Neucr2") for i in range(1, 4)]
    named = row("chs1", "10_Neucr2")
    names = {str(i): (f"NCU0000{i}T0", S1) for i in range(1, 4)} | {"10": ("NCU12345T0", S2)}
    mk.apply_jgi([*agree, named], {"Neucr2": names}, {"Neucr2"})
    assert named["jgi_check"] == "unverifiable"
    assert named["jgi_accept"] == "portal_validated" and named["map_status"] == "mapped_jgi_file"

    lone = row("chs2", "11_Hisca1")
    mk.apply_jgi([lone], {"Hisca1": {"11": ("HCAG_1T0", S1)}}, {"Hisca1"})
    assert lone["jgi_accept"] == "" and lone["map_status"] == "lost"
    assert lone["map_reason"].endswith("jgi_file_unverifiable_portal_not_validated")


def test_a_portal_with_too_many_mismatches_is_not_validated():
    rows = [row(f"NCU0000{i}", f"{i}_Neucr2") for i in range(1, 5)]
    names = {
        "1": ("NCU00001T0", S1),
        "2": ("NCU00002T0", S1),
        "3": ("NCU00003T0", S1),
        "4": ("NCU77777T0", S1),
    }
    named = row("chs1", "10_Neucr2")
    names["10"] = ("NCU12345T0", S2)
    stats = mk.apply_jgi([*rows, named], {"Neucr2": names}, {"Neucr2"})
    assert stats["Neucr2"]["agree"] == 3 and stats["Neucr2"]["differs"] == 1
    assert stats["Neucr2"]["portal_validated"] == 0  # 1 of 4 is above the 10 percent limit
    assert named["map_status"] == "lost"


def test_implausible_jgi_sequence_is_not_mapped():
    r = row("bgs1", "1_Schpo1")
    mk.apply_jgi([r], {"Schpo1": {"1": ("bgs1", "M" + "X" * 40)}}, {"Schpo1"})
    assert r["map_status"] == "lost" and "jgi_sequence_implausible" in r["map_reason"]


def test_guard_refuses_a_table_that_carries_a_jgi_sequence(tmp_path):
    ok = tmp_path / "ok.tsv"
    ok.write_text("a\tb\n1\t2\n")
    mk.guard_no_sequences([ok], [S1])
    bad = tmp_path / "bad.tsv"
    bad.write_text(f"a\n{S1}\n")
    with pytest.raises(RuntimeError, match="JGI sequence"):
        mk.guard_no_sequences([ok, bad], [S1])


def test_written_row_table_never_has_a_sequence_column():
    assert "sequence" not in mk.JGI_COLS
    assert all("seq" not in c or c == "jgi_seq_sha256" for c in mk.JGI_COLS)
