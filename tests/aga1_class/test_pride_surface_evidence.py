"""Tests of pride_surface_evidence.py (SYNTHETIC protXML and tables) and of the committed PRIDE record."""

import csv
import hashlib
import importlib.util
import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "pride", ROOT / "analysis/aga1_class/pride_surface_evidence.py"
)
pr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pr)
REC = ROOT / "docs/reports/data/aga1_pride_conidial_surface_2026-10-26"

XML = """<?xml version="1.0"?>
<protein_summary xmlns="http://regis-web.systemsbiology.net/protXML">
<protein_summary_header reference_database="AfumFDB_201102.fasta;0;x" />
<protein_group group_number="0"><protein protein_name="Afu6g07910-T-p1" percent_coverage="5.0">
<annotation protein_description="transcript=T | gene_product=PH domain protein | location=Chr6" />
<peptide peptide_sequence="AAAK" n_instances="4"/><peptide peptide_sequence="CCCK" n_instances="2"/><peptide peptide_sequence="AAAK" n_instances="1"/>
</protein></protein_group>
<protein_group group_number="1"><protein protein_name="Afu1g00010-T-p1" percent_coverage="2.0"><peptide peptide_sequence="GGGK" n_instances="3"/></protein></protein_group>
</protein_summary>"""


def test_locus_tag_conversion():
    assert (
        pr.afu_to_afua("Afu6g07910-T-p1") == "AFUA_6G07910"
        and pr.afu_to_afua("Afu2g00010") == "AFUA_2G00010"
    )
    assert pr.afu_to_afua("jgi|Aspoer1|1") == "" and pr.afu_to_afua("NFIA_011720") == ""


def test_parse_protxml_counts_distinct_peptides_and_spectra(tmp_path):
    f = tmp_path / "1.prot.xml"
    f.write_text(XML)
    species, rows = pr.parse_protxml(f)
    assert species == "A. fumigatus Af293" and len(rows) == 2
    r = rows[0]
    assert (r["protein_id"], r["n_peptides"], r["spectral_count"], r["gene_product"]) == (
        "Afu6g07910-T-p1",
        2,
        7,
        "PH domain protein",
    )


def test_parse_protxml_rejects_unknown_database_and_header_less(tmp_path):
    f = tmp_path / "x.prot.xml"
    f.write_text(XML.replace("AfumFDB", "Other"))
    with pytest.raises(pr.PrideError, match="unknown reference database"):
        pr.parse_protxml(f)
    f.write_text(
        '<protein_summary xmlns="http://regis-web.systemsbiology.net/protXML"></protein_summary>'
    )
    with pytest.raises(pr.PrideError, match="no protein_summary_header"):
        pr.parse_protxml(f)


def test_detections_suffix_and_replicate(tmp_path):
    for n in ("1", "1s", "2"):
        (tmp_path / f"{n}.prot.xml").write_text(XML)
    d = pr.detections(tmp_path)
    assert set(zip(d["run"], d["run_suffix"], d["replicate"], strict=True)) == {
        ("1", "plain", 1),
        ("1s", "s", 1),
        ("2", "plain", 2),
    }
    with pytest.raises(pr.PrideError, match="no"):
        pr.detections(tmp_path / "empty_dir_that_does_not_exist")


def test_kmer_map_picks_the_unique_best_and_leaves_ties_and_poor_matches_empty():
    a = "MKVLAAGIVGLLLAQPAMSTEEDKLRRHHWQ" * 2
    af = {"X1": a, "X2": "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" * 2}
    assert pr.kmer_map(["P"], {"P": a}, af)["P"][0] == "X1"
    assert pr.kmer_map(["P"], {"P": a}, {"X1": a, "X3": a})["P"][0] == ""  # tied best match
    assert (
        pr.kmer_map(["P"], {"P": "MCDEFGHIKLNPQRSTVWYACDEFGHIKLNPQR"}, af)["P"][0] == ""
    )  # nothing shared


def flags(rows):
    cols = [
        "accession",
        "species",
        "sp_yes",
        "gpi_yes",
        "broad_flag",
        "narrow_flag",
        "n_models_call",
    ]
    return pd.DataFrame(rows, columns=cols).astype({"accession": str, "species": str})


def test_join_and_enrichment_on_synthetic_data():
    runs = pd.DataFrame(
        [
            {
                "run": r,
                "run_suffix": "plain",
                "replicate": 1,
                "species": "A. fumigatus Af293",
                "protein_id": "Afu1g00010-T-p1",
                "n_peptides": 3,
                "spectral_count": 5,
                "percent_coverage": "1",
                "gene_product": "gp",
            }
            for r in ("1", "2")
        ]  # fmt: skip
    )
    names = pd.DataFrame({"Entry": ["A1", "A2", "A3"], "Gene Names (ORF)": ["AFUA_1G00010", "AFUA_1G00020", "AFUA_1G00030"],
                          "Protein names": ["n1", "n2", "n3"]})  # fmt: skip
    fl = flags(
        [
            ("A1", "Afum_Af293", 1, 1, 1, 0, 3),
            ("A2", "Afum_Af293", 0, 0, 0, 0, 0),
            ("A3", "Afum_Af293", 1, 0, 0, 0, 1),
        ]
    )
    ev, unmapped, _ = pr.join(runs, names, fl)
    e = ev.set_index("accession")
    assert (
        unmapped == 0
        and e.loc["A1", "runs_detected"] == 2
        and e.loc["A1", "detected_2plus_runs"] == 1
    )
    assert e.loc["A2", "runs_detected"] == 0 and e.loc["A2", "total_spectra"] == 0
    enr = pr.enrichment(ev)
    sp = enr[(enr.protein_set == "signal peptide") & (enr.detection == "detected_any")].iloc[0]
    assert (sp.n_in_set, sp.detected_in_set, sp.detected_outside) == (2, 1, 0)


def test_committed_record_hashes_and_headline_numbers():
    text = (REC / "README.md").read_text()
    found = {p: s for s, p in re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)}
    committed = {p.name for p in REC.iterdir() if p.is_file() and p.name != "README.md"}
    assert committed <= set(found)
    for p in committed:
        assert hashlib.sha256((REC / p).read_bytes()).hexdigest() == found[p], p
    with (REC / "afum_pride_surface_evidence.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 9647
    assert (
        sum(r["detected_any"] == "1" for r in rows) == 380
        and sum(r["detected_2plus_runs"] == "1" for r in rows) == 288
    )
    top = {r["accession"]: r for r in rows}[
        "Q4WPK6"
    ]  # AFUA_4G09600, a tier 1 candidate seen in every run
    assert top["afua"] == "AFUA_4G09600" and top["runs_detected"] == "6"
