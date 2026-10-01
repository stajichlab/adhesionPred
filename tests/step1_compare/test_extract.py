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


def _recount(rows):
    """Recompute the label counts of counts.tsv from rows read back from truth_set.tsv.gz."""

    def n(column, value):
        return str(sum(1 for r in rows if r[column] == value))

    return {
        "genes_cc": str(len(rows)),
        "p_ext": n("label", "P-ext"),
        "p_ext_wall": n("subset", "wall"),
        "p_ext_extonly": n("subset", "extracellular-only"),
        "n_int": n("label", "N-int"),
        "n_sec": n("label", "N-sec"),
        "ambiguous": n("label", "ambiguous"),
        "pm_candidates": n("pm_candidate", "yes"),
        "exp_p_ext": n("label_experimental", "P-ext"),
        "exp_n_int": n("label_experimental", "N-int"),
        "exp_n_sec": n("label_experimental", "N-sec"),
        "exp_ambiguous": n("label_experimental", "ambiguous"),
        "nohom_p_ext": n("label_no_homology", "P-ext"),
        "nohom_n_int": n("label_no_homology", "N-int"),
        "nohom_n_sec": n("label_no_homology", "N-sec"),
        "nohom_ambiguous": n("label_no_homology", "ambiguous"),
        "direct_p_ext": str(
            sum(r["label"] == "P-ext" and r["homology_only"] == "no" for r in rows)
        ),
        "direct_n_int": str(
            sum(r["label"] == "N-int" and r["homology_only"] == "no" for r in rows)
        ),
        "direct_n_sec": str(
            sum(r["label"] == "N-sec" and r["homology_only"] == "no" for r in rows)
        ),
        "direct_ambiguous": str(
            sum(r["label"] == "ambiguous" and r["homology_only"] == "no" for r in rows)
        ),
        "ambiguous_htp_only": str(
            sum(
                r["label"] == "ambiguous" and r["internal_evidence_htp_only"] == "yes" for r in rows
            )
        ),
    }


def _two_source_inputs(tmp_path, fixtures_dir):
    """Fixture inputs with a second source (a copy of the golden GAF) after the first."""
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    shutil.copy(input_dir / "golden.gaf", input_dir / "golden2.gaf")
    manifest_rows.append(
        {
            "file": "golden2.gaf",
            "mode": "strict",
            "sha256": manifest.sha256_file(input_dir / "golden2.gaf"),
        }
    )
    species = FIXTURE_SPECIES + [
        dict(FIXTURE_SPECIES[0], source_id="Fix_SGD2", gaf_file="golden2.gaf")
    ]
    return input_dir, manifest_rows, species


def _cli(tmp_path, species, manifest_rows, input_dir, out, extra=()):
    species_path = tmp_path / "species.tsv"
    truth_table.write_tsv(species_path, list(species[0]), species)
    mpath = tmp_path / "manifest.tsv"
    full = [{c: "" for c in manifest.MANIFEST_COLUMNS} | r for r in manifest_rows]
    manifest.write_manifest(mpath, full)
    argv = ["--species", str(species_path), "--manifest", str(mpath)]
    argv += ["--input-dir", str(input_dir), "--out-dir", str(out), *extra]
    return extract.main(argv), argv


def _files(directory):
    return sorted(p.name for p in directory.iterdir()) if directory.exists() else []


def test_truth_file_reads_back_with_expected_content(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "out")
    truth = truth_table.read_tsv(tmp_path / "out" / "truth_set.tsv.gz")
    by_gene = {r["gene_id"]: r for r in truth}
    assert by_gene["G01"]["label"] == "P-ext" and by_gene["G01"]["subset"] == "wall"
    assert by_gene["G02"]["subset"] == "extracellular-only"
    assert by_gene["G04"]["label"] == "ambiguous"
    assert by_gene["G05"]["label"] == "N-int"
    assert by_gene["G07"]["label"] == "N-sec"
    assert by_gene["G10"]["label"] == "ambiguous"
    assert by_gene["G10"]["label_no_homology"] == "P-ext"
    assert by_gene["G10"]["homology_only"] == "yes"
    assert by_gene["G14"]["internal_evidence_htp_only"] == "yes"
    for row in truth:
        assert row["in_clade"] == "Saccharomycotina"
        assert row["role"] == "train"
        assert row["obo_sha256"] == manifest_rows[0]["sha256"]
    counts = truth_table.read_tsv(tmp_path / "out" / "counts.tsv")[0]
    recount = _recount(truth)
    assert {k: counts[k] for k in recount} == recount


def test_output_gzip_header_has_no_name_or_time(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "out")
    data = (tmp_path / "out" / "truth_set.tsv.gz").read_bytes()
    assert data[:2] == b"\x1f\x8b"
    assert data[3] == 0  # FLG: no FNAME, no FEXTRA, no FCOMMENT
    assert data[4:8] == b"\0\0\0\0"  # MTIME


def test_output_does_not_depend_on_file_name_or_directory(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "first")
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "x" / "other_dir_name")
    a = (tmp_path / "first" / "truth_set.tsv.gz").read_bytes()
    assert a == (tmp_path / "x" / "other_dir_name" / "truth_set.tsv.gz").read_bytes()
    rows = [{"a": "1", "b": "2"}]
    truth_table.write_tsv(tmp_path / "name_one.tsv.gz", ["a", "b"], rows)
    truth_table.write_tsv(tmp_path / "different_name.tsv.gz", ["a", "b"], rows)
    assert (tmp_path / "name_one.tsv.gz").read_bytes() == (
        tmp_path / "different_name.tsv.gz"
    ).read_bytes()


def test_log_is_byte_stable_and_complete(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    code, argv = _cli(tmp_path, FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "a")
    assert code == 0
    log_a = (tmp_path / "a" / "extract_log.json").read_bytes()
    assert extract.main(argv) == 0  # same arguments, same output directory
    assert log_a == (tmp_path / "a" / "extract_log.json").read_bytes()
    import json

    log = json.loads(log_a)
    assert log["species_sha256"] == manifest.sha256_file(tmp_path / "species.tsv")
    assert log["manifest_sha256"] == manifest.sha256_file(tmp_path / "manifest.tsv")
    assert log["input_dir"] == str(input_dir)
    assert log["python"] and log["git_commit"]
    assert log["arguments"][:2] == argv[:2]
    source = log["sources"][0]
    assert source["date_generated"] == "2026-05-21T09:00"
    assert source["unknown_term_rows"] == 0
    # run() on its own (fixture, no CLI) also gives identical logs for identical inputs.
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "c")
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "d")
    assert (tmp_path / "c" / "extract_log.json").read_bytes() == (
        tmp_path / "d" / "extract_log.json"
    ).read_bytes()


@pytest.mark.parametrize("victim", ["golden.gaf", "golden2.gaf", "go-basic.obo"])
def test_hash_mismatch_writes_no_output_and_cli_stops(tmp_path, fixtures_dir, victim, capsys):
    input_dir, manifest_rows, species = _two_source_inputs(tmp_path, fixtures_dir)
    with open(input_dir / victim, "a") as handle:
        handle.write("\n! tampered\n" if victim.endswith("obo") else "!tampered\n")
    out = tmp_path / "out"
    with pytest.raises(manifest.DownloadError, match="differs from manifest"):
        extract.run(species, manifest_rows, input_dir, out)
    assert _files(out) == []
    code, _ = _cli(tmp_path, species, manifest_rows, input_dir, out)
    assert code == 2
    assert "STOP:" in capsys.readouterr().err
    assert _files(out) == []


def test_missing_gaf_writes_no_output(tmp_path, fixtures_dir):
    input_dir, manifest_rows, species = _two_source_inputs(tmp_path, fixtures_dir)
    (input_dir / "golden2.gaf").unlink()
    out = tmp_path / "out"
    code, _ = _cli(tmp_path, species, manifest_rows, input_dir, out)
    assert code == 2
    assert _files(out) == []


def test_missing_manifest_row_stops_with_exit_2(tmp_path, fixtures_dir, capsys):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    code, _ = _cli(tmp_path, FIXTURE_SPECIES, manifest_rows[:1], input_dir, out)
    assert code == 2
    assert "golden.gaf: no manifest row" in capsys.readouterr().err
    assert _files(out) == []
    with pytest.raises(manifest.DownloadError, match="go-basic.obo: no manifest row"):
        extract.run(FIXTURE_SPECIES, manifest_rows[1:], input_dir, out)


def test_malformed_gaf_stops_with_exit_2(tmp_path, fixtures_dir, capsys):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    (input_dir / "golden.gaf").write_text("!gaf-version: 2.2\n")
    manifest_rows[1]["sha256"] = manifest.sha256_file(input_dir / "golden.gaf")
    out = tmp_path / "out"
    code, _ = _cli(tmp_path, FIXTURE_SPECIES, manifest_rows, input_dir, out)
    assert code == 2
    assert "STOP:" in capsys.readouterr().err
    assert _files(out) == []


@pytest.mark.parametrize("extra", [["--sources", "Fix_SDG"], ["--sources"]])
def test_unknown_or_empty_source_selection_stops(tmp_path, fixtures_dir, capsys, extra):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    code, _ = _cli(tmp_path, FIXTURE_SPECIES, manifest_rows, input_dir, out, extra)
    assert code == 2
    err = capsys.readouterr().err
    if len(extra) > 1:
        assert "STOP: unknown source Fix_SDG" in err and "Fix_SGD" in err
    else:
        assert "STOP:" in err
    assert _files(out) == []


def test_unknown_source_keeps_previous_outputs(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    assert _cli(tmp_path, FIXTURE_SPECIES, manifest_rows, input_dir, out)[0] == 0
    before = {n: (out / n).read_bytes() for n in extract.OUTPUT_NAMES}
    code, _ = _cli(tmp_path, FIXTURE_SPECIES, manifest_rows, input_dir, out, ["--sources", "no"])
    assert code == 2
    assert {n: (out / n).read_bytes() for n in extract.OUTPUT_NAMES} == before


def test_failed_write_keeps_previous_outputs_and_leaves_no_temp_files(
    tmp_path, fixtures_dir, monkeypatch
):
    input_dir, manifest_rows, species = _two_source_inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    extract.run(species, manifest_rows, input_dir, out)
    before = {n: (out / n).read_bytes() for n in extract.OUTPUT_NAMES}
    real_write = truth_table.write_tsv
    calls = []

    def failing_write(path, columns, rows):
        calls.append(path)
        if len(calls) == 2:  # the truth table is written; counts.tsv fails
            raise OSError("disk full")
        real_write(path, columns, rows)

    monkeypatch.setattr(truth_table, "write_tsv", failing_write)
    with pytest.raises(OSError, match="disk full"):
        extract.run(FIXTURE_SPECIES[:1], manifest_rows, input_dir, out)
    assert {n: (out / n).read_bytes() for n in extract.OUTPUT_NAMES} == before
    assert _files(out) == sorted(extract.OUTPUT_NAMES)
    assert len(calls) == 2


def test_failed_first_write_leaves_empty_directory(tmp_path, fixtures_dir, monkeypatch):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    real_write = truth_table.write_tsv

    def failing_write(path, columns, rows):
        if Path(path).name.endswith("counts.tsv"):
            raise OSError("disk full")
        real_write(path, columns, rows)

    monkeypatch.setattr(truth_table, "write_tsv", failing_write)
    with pytest.raises(OSError):
        extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, out)
    assert _files(out) == []


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
def real_run(tmp_path_factory):
    species = extract.read_species(paths.STEP1_DIR / "species.tsv")
    manifest_rows = manifest.read_manifest(paths.STEP1_DIR / "manifest.tsv")
    out = tmp_path_factory.mktemp("real")
    _, counts = extract.run(species, manifest_rows, _real_input_dir(), out)
    return out, {c["source_id"]: c for c in counts}


@pytest.fixture(scope="module")
def real_counts(real_run):
    return real_run[1]


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


# exp_*: "exp-only" values printed by /tmp/glyco_spec/d1_count.py (2026-10-01).
# direct_ambiguous: d1_count.py "no-homology-codes ambig" (equal to direct_ambiguous here).
# ambiguous_htp_only: S288C 18 and S. pombe 4 are the spec values. The other sources are 0 in
# the first real run (no external source for them).
EXP_AND_AMBIGUOUS_COUNTS = {
    # source: (exp_p_ext, exp_n_int, exp_n_sec, exp_ambiguous, direct_ambiguous, htp_only)
    "Scer_SGD": (93, 2249, 1437, 32, 35, 18),
    "Calb_CGD": (292, 168, 260, 13, 14, 0),
    "Spom_PomBase": (32, 2715, 831, 2, 3, 4),
    "Spom_SCHPO-mod": (31, 2721, 836, 2, 3, 4),
    "Cneo_H99_GOA": (9, 14, 14, 0, 0, 0),
    "Cneo_JEC21_GOA": (0, 2, 1, 0, 0, 0),
    "Cneo_CRYD1": (0, 2, 1, 0, 0, 0),
    "Afum_ASPFU": (23, 18, 26, 1, 1, 0),
    "Anid_EMENI": (148, 91, 66, 3, 3, 0),
    "Umay_MYCMD": (10, 6, 21, 1, 1, 0),
}


@needs_real_files
def test_exp_ambiguous_and_htp_counts_on_real_files(real_counts):
    keys = (
        "exp_p_ext",
        "exp_n_int",
        "exp_n_sec",
        "exp_ambiguous",
        "direct_ambiguous",
        "ambiguous_htp_only",
    )
    for source_id, expected in EXP_AND_AMBIGUOUS_COUNTS.items():
        got = tuple(real_counts[source_id][k] for k in keys)
        assert got == tuple(str(v) for v in expected), source_id


@needs_real_files
def test_real_truth_file_reads_back_and_matches_counts(real_run):
    out, counts = real_run
    truth = truth_table.read_tsv(out / "truth_set.tsv.gz")
    on_disk = {c["source_id"]: c for c in truth_table.read_tsv(out / "counts.tsv")}
    assert on_disk == counts
    by_source = {}
    for row in truth:
        by_source.setdefault(row["source_id"], []).append(row)
    assert set(by_source) == set(counts)
    for source_id, rows in by_source.items():
        assert _recount(rows) == {k: counts[source_id][k] for k in _recount(rows)}, source_id


def test_log_records_whether_all_sources_were_processed(tmp_path, fixtures_dir):
    import json

    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    assert _cli(tmp_path, FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "a")[0] == 0
    log = json.loads((tmp_path / "a" / "extract_log.json").read_text())
    assert log["all_sources"] is True
    two = FIXTURE_SPECIES + [{**FIXTURE_SPECIES[0], "source_id": "Fix_Two"}]
    code, _ = _cli(
        tmp_path, two, manifest_rows, input_dir, tmp_path / "b", ["--sources", "Fix_SGD"]
    )
    assert code == 0
    log = json.loads((tmp_path / "b" / "extract_log.json").read_text())
    assert log["all_sources"] is False
    assert [s["source_id"] for s in log["sources"]] == ["Fix_SGD"]
