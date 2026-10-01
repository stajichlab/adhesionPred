import os
import shutil
from pathlib import Path

import manifest
import paths
import pytest
import truth_table
from conftest import load_script
from gaf_fixture import write_golden_gaf

extract = load_script("01_extract_go_truth")

FIXTURE_SPECIES = [
    {
        "source_id": "Fix_SGD",
        "species": "Fixture yeast",
        "taxon_id": "559292",
        "taxon_filter": "",
        "in_clade": "Saccharomycotina",
        "role": "train",
        "gaf_file": "golden.gaf",
    }
]


def _inputs(tmp_path, fixtures_dir):
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    shutil.copy(fixtures_dir / "mini.obo", input_dir / "go-basic.obo")
    write_golden_gaf(input_dir / "golden.gaf")
    rows = [
        {"file": f, "mode": "strict", "sha256": manifest.sha256_file(input_dir / f)}
        for f in ("go-basic.obo", "golden.gaf")
    ]
    return input_dir, rows


def test_extract_writes_truth_and_counts(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, out)
    truth = truth_table.read_tsv(out / "truth_set.tsv.gz")
    assert tuple(truth[0].keys()) == truth_table.TRUTH_COLUMNS
    assert len(truth) == 20
    counts = truth_table.read_tsv(out / "counts.tsv")[0]
    expected = {
        "primary_db": "SGD",
        "genes_cc": "20",
        "genes_noniea_cc": "19",
        "p_ext": "7",
        "p_ext_wall": "4",
        "p_ext_extonly": "3",
        "n_int": "1",
        "n_sec": "4",
        "ambiguous": "3",
        "pm_candidates": "1",
        "cc_iea_triples": "7",
        "cc_triples": "35",
        "cc_iea_frac": "0.200",
        "all_aspects_iea_frac": "0.209",
        "obsolete_rows": "1",
        "unknown_term_rows": "0",
        "exp_p_ext": "6",
        "exp_n_int": "1",
        "exp_n_sec": "4",
        "exp_ambiguous": "2",
        "nohom_p_ext": "6",
        "nohom_n_int": "1",
        "nohom_n_sec": "4",
        "nohom_ambiguous": "2",
        "direct_p_ext": "5",
        "direct_n_int": "1",
        "direct_n_sec": "4",
        "direct_ambiguous": "2",
        "ambiguous_htp_only": "1",
    }
    assert {k: counts[k] for k in expected} == expected
    assert truth[0]["source_sha256"] == manifest_rows[1]["sha256"]
    assert truth[0]["source_date"] == "2026-05-21T09:00"


def test_extract_output_is_byte_stable(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "a")
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "b")
    a = (tmp_path / "a" / "truth_set.tsv.gz").read_bytes()
    assert a == (tmp_path / "b" / "truth_set.tsv.gz").read_bytes()


def test_extract_stops_on_hash_mismatch(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    with open(input_dir / "golden.gaf", "a") as handle:
        handle.write(
            "SGD\tG99\tNEW\tlocated_in\tGO:0005576\tPMID:1\tIDA\t\tC\tn\t\tprotein"
            "\ttaxon:559292\t20260101\tSGD\t\t\n"
        )
    with pytest.raises(manifest.DownloadError, match="differs from manifest"):
        extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "out")
    assert not (tmp_path / "out" / "truth_set.tsv.gz").exists()


def test_missing_gaf_stops_with_exit_2(tmp_path, fixtures_dir, capsys):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    (input_dir / "golden.gaf").unlink()
    species_path = tmp_path / "species.tsv"
    truth_table.write_tsv(species_path, list(FIXTURE_SPECIES[0]), FIXTURE_SPECIES)
    mpath = tmp_path / "manifest.tsv"
    full = [{c: "" for c in manifest.MANIFEST_COLUMNS} | r for r in manifest_rows]
    manifest.write_manifest(mpath, full)
    argv = ["--species", str(species_path), "--manifest", str(mpath)]
    argv += ["--input-dir", str(input_dir), "--out-dir", str(tmp_path / "out")]
    assert extract.main(argv) == 2
    err = capsys.readouterr().err
    assert "STOP:" in err and "golden.gaf" in err and "file not found" in err


# Spec 3.3 numbers, reproduced from /tmp/glyco_spec/d1_count.py on 2026-10-01.
SPEC_COUNTS = {
    "Scer_SGD": {
        "genes_cc": 6056,
        "genes_noniea_cc": 6012,
        "p_ext": 125,
        "p_ext_wall": 104,
        "p_ext_extonly": 21,
        "n_int": 2645,
        "n_sec": 1548,
        "ambiguous": 48,
        "cc_iea_frac": "0.262",
        "nohom_p_ext": 96,
        "nohom_n_int": 2360,
        "nohom_n_sec": 1457,
        "pm_candidates": 12,
    },
    "Calb_CGD": {
        "genes_cc": 6313,
        "genes_noniea_cc": 6160,
        "p_ext": 259,
        "p_ext_wall": 122,
        "p_ext_extonly": 137,
        "n_int": 1772,
        "n_sec": 961,
        "ambiguous": 100,
        "cc_iea_frac": "0.158",
        "nohom_p_ext": 294,
        "nohom_n_int": 179,
        "nohom_n_sec": 263,
        "pm_candidates": 60,
    },
    "Spom_PomBase": {
        "genes_cc": 5025,
        "genes_noniea_cc": 4998,
        "p_ext": 59,
        "p_ext_wall": 40,
        "p_ext_extonly": 19,
        "n_int": 2962,
        "n_sec": 1214,
        "ambiguous": 14,
        "cc_iea_frac": "0.292",
        "nohom_p_ext": 44,
        "nohom_n_int": 2799,
        "nohom_n_sec": 886,
        "pm_candidates": 10,
    },
    "Spom_SCHPO-mod": {"p_ext": 57, "n_int": 2966, "n_sec": 1212, "ambiguous": 14},
    "Cneo_H99_GOA": {
        "genes_cc": 4320,
        "genes_noniea_cc": 77,
        "p_ext": 11,
        "p_ext_wall": 3,
        "p_ext_extonly": 8,
        "n_int": 23,
        "n_sec": 25,
        "ambiguous": 0,
        "cc_iea_frac": "0.986",
        "nohom_p_ext": 9,
        "nohom_n_int": 17,
        "nohom_n_sec": 15,
    },
    "Cneo_JEC21_GOA": {
        "p_ext": 32,
        "p_ext_wall": 10,
        "p_ext_extonly": 22,
        "n_int": 1813,
        "n_sec": 817,
        "ambiguous": 0,
        "cc_iea_frac": "0.520",
        "nohom_p_ext": 0,
        "nohom_n_int": 3,
        "nohom_n_sec": 1,
    },
    "Cneo_CRYD1": {"p_ext": 32, "n_int": 1813, "n_sec": 817, "ambiguous": 0},
    "Afum_ASPFU": {
        "p_ext": 132,
        "p_ext_wall": 35,
        "p_ext_extonly": 97,
        "n_int": 2018,
        "n_sec": 1062,
        "ambiguous": 1,
        "cc_iea_frac": "0.533",
        "nohom_p_ext": 23,
        "nohom_n_int": 18,
        "nohom_n_sec": 26,
        "pm_candidates": 9,
    },
    "Anid_EMENI": {
        "p_ext": 211,
        "p_ext_wall": 43,
        "p_ext_extonly": 168,
        "n_int": 2054,
        "n_sec": 1066,
        "ambiguous": 38,
        "cc_iea_frac": "0.517",
        "nohom_p_ext": 148,
        "nohom_n_int": 97,
        "nohom_n_sec": 66,
    },
    "Umay_MYCMD": {
        "p_ext": 62,
        "p_ext_wall": 14,
        "p_ext_extonly": 48,
        "n_int": 1712,
        "n_sec": 827,
        "ambiguous": 1,
        "cc_iea_frac": "0.519",
        "nohom_p_ext": 10,
        "nohom_n_int": 6,
        "nohom_n_sec": 21,
    },
}


# M2: genes_cc for the five UniProt-based sources.
GENES_CC_UNIPROT = {
    "Cneo_H99_GOA": 4320,
    "Cneo_JEC21_GOA": 4364,
    "Cneo_CRYD1": 4317,
    "Afum_ASPFU": 5613,
    "Anid_EMENI": 5884,
    "Umay_MYCMD": 4206,
}

# Spec 3.3 direct-evidence truth: label == X and homology_only == "no" (P-ext, N-int, N-sec).
DIRECT_COUNTS = {
    "Scer_SGD": (88, 2360, 1457),
    "Calb_CGD": (211, 179, 263),
    "Spom_PomBase": (42, 2799, 886),
    "Cneo_H99_GOA": (9, 17, 15),
    "Cneo_JEC21_GOA": (0, 3, 1),
    "Afum_ASPFU": (23, 18, 26),
    "Anid_EMENI": (113, 97, 66),
    "Umay_MYCMD": (10, 6, 21),
}


def _real_input_dir() -> Path:
    return Path(os.environ.get("STEP1_GO_DIR") or paths.downloads_dir())


needs_real_files = pytest.mark.skipif(
    not (_real_input_dir() / "go-basic.obo").exists(),
    reason="real GO files absent; run 00_fetch_inputs.py or set STEP1_GO_DIR",
)


@pytest.fixture(scope="module")
def real_counts(tmp_path_factory):
    species = extract.read_species(paths.STEP1_DIR / "species.tsv")
    manifest_rows = manifest.read_manifest(paths.STEP1_DIR / "manifest.tsv")
    out = tmp_path_factory.mktemp("real")
    _, counts = extract.run(species, manifest_rows, _real_input_dir(), out)
    return {c["source_id"]: c for c in counts}


@needs_real_files
def test_reproduces_spec_counts_on_real_files(real_counts):
    for source_id, expected in SPEC_COUNTS.items():
        for key, value in expected.items():
            assert real_counts[source_id][key] == str(value), (source_id, key)
    for source_id, value in GENES_CC_UNIPROT.items():
        assert real_counts[source_id]["genes_cc"] == str(value), source_id


@needs_real_files
def test_direct_evidence_counts_on_real_files(real_counts):
    for source_id, (p_ext, n_int, n_sec) in DIRECT_COUNTS.items():
        got = real_counts[source_id]
        assert (got["direct_p_ext"], got["direct_n_int"], got["direct_n_sec"]) == (
            str(p_ext),
            str(n_int),
            str(n_sec),
        ), source_id
