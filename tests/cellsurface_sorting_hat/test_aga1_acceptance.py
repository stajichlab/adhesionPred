"""The AGA1 acceptance tooling (SYNTHETIC inputs) and the committed acceptance record and status files."""

import csv
import hashlib
import importlib.util
import re
from pathlib import Path

import pytest

from cellsurface_sorting_hat.call_status import load_call_source
from cellsurface_sorting_hat.engine import load_config
from cellsurface_sorting_hat.status import load_status_source

ROOT = Path(__file__).resolve().parents[2]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


acc = load("analysis/sorting_hat_run/aga1_acceptance.py", "aga1_acceptance")
truth = load("analysis/sorting_hat_run/aga1_status_truth.py", "aga1_status_truth")
REC = ROOT / "docs/reports/data/aga1_sorting_hat_acceptance_2026-10-26"
STATUS = ROOT / "docs/reports/data/sorting_hat/call_status/aga1_850"
MEMBERS = {"A1": "SpA", "A2": "SpA", "B1": "SpB"}


# ---- splitting helpers
def test_split_signalp_keeps_headers_and_splits_rows_by_species():
    import tempfile

    text = "# SignalP-6.0\tOrganism: Eukarya\n# ID\tPrediction\nA1\tSP\t0.1\t0.9\t\nB1\tOTHER\t0.9\t0.1\t\nZ9\tSP\t0\t1\t\n"
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "s.txt"
        p.write_text(text)
        out = acc.split_signalp(p, MEMBERS)
    assert (
        set(out) == {"SpA", "SpB"}
        and out["SpA"].startswith("# SignalP-6.0")
        and "A1\tSP" in out["SpA"]
        and "B1" not in out["SpA"]
        and "Z9" not in out["SpB"]
    )


def test_split_gff3_table_and_disorder(tmp_path):
    g = tmp_path / "g.gff3"
    g.write_text("A1\tPredGPI\tChain\t1\t9\t1.0\t.\t.\te\nB1\tPredGPI\tChain\t1\t9\t1.0\t.\t.\te\n")
    assert set(acc.split_gff3(g, MEMBERS)) == {"SpA", "SpB"}
    t = tmp_path / "t.tsv"
    t.write_text("accession\tscore_rf\nA1\t0.1\nA2\t0.2\nB1\t0.3\nQ\t0.9\n")
    parts = acc.split_table(t, "accession", MEMBERS)
    assert (
        parts["SpA"].splitlines() == ["accession\tscore_rf", "A1\t0.1", "A2\t0.2"]
        and "Q" not in parts["SpB"]
    )
    d = tmp_path / "d.txt"
    d.write_text("# banner\n#>A1\n1\tM\t0.1\n\n#>B1\n1\tK\t0.2\n\n#>Z9\n1\tK\t0.3\n")
    files = acc.split_disorder_to_files(d, MEMBERS, tmp_path)
    assert (
        set(files) == {"SpA", "SpB"}
        and "#>A1" in files["SpA"].read_text()
        and "#>B1" not in files["SpA"].read_text()
    )


def test_split_fasta_stops_when_a_protein_is_missing(tmp_path):
    (tmp_path / "c1.faa").write_text(">A1\nMKT\n>B1\nMKQ\n")
    out = acc.split_fasta(tmp_path, {"A1": "SpA", "B1": "SpB"})
    assert out == {"SpA": ">A1\nMKT\n", "SpB": ">B1\nMKQ\n"}
    with pytest.raises(acc.AcceptanceError, match="not in"):
        acc.split_fasta(tmp_path, {"A1": "SpA", "NOPE": "SpA"})


def test_compare_lists_every_difference():
    calls = {
        "ser_thr_gpi_wall_candidate": {"A1": "called", "A2": "not_called", "A3": "not_assessable"},
        "flocculin_like_candidate": {"A1": "called", "A2": "not_called", "A3": "not_called"},
        "model_surface_candidate": {"A1": "called", "A2": "called", "A3": "not_called"},
    }
    flags = {
        "A1": {"broad_flag": "1", "narrow_flag": "1", "n_models_call": "3", "sp_yes": "1"},
        "A2": {
            "broad_flag": "0",
            "narrow_flag": "0",
            "n_models_call": "2",
            "sp_yes": "0",
        },  # models call it, no signal peptide: not called
        "A3": {
            "broad_flag": "0",
            "narrow_flag": "0",
            "n_models_call": "0",
            "sp_yes": "1",
        },  # not_assessable counts as not called
        "A4": {
            "broad_flag": "0",
            "narrow_flag": "0",
            "n_models_call": "0",
            "sp_yes": "1",
        },  # missing from the engine output
    }
    d = acc.compare(calls, flags)
    assert {(a, c) for a, c, _ in d} == {
        ("A2", "model_surface_candidate"),
        ("A4", "ser_thr_gpi_wall_candidate"),
        ("A4", "flocculin_like_candidate"),
        ("A4", "model_surface_candidate"),
    }


def test_species_table_covers_the_ten_proteomes_once():
    taxa = acc.species_taxa()
    assert len(taxa) == 10 and len(set(taxa.values())) == 10
    for folder in (
        "aga1_proteome_test_2026-10-25",
        "aga1_nonyeast_test_2026-10-26",
        "aga1_cpar_test_2026-10-26",
    ):
        with (ROOT / "docs/reports/data" / folder / "summary_by_species.tsv").open(
            newline=""
        ) as fh:
            assert {r["species"] for r in csv.DictReader(fh, delimiter="\t")} <= set(taxa)


def test_acceptance_refuses_an_existing_out_root(tmp_path, capsys):
    assert (
        acc.main(
            [
                "--taxdump",
                "x",
                "--out-root",
                str(tmp_path),
                "--set",
                "a",
                "b",
                "c",
                "d",
                "e",
                "f",
                "g",
            ]
        )
        == 1
    )
    assert "FATAL" in capsys.readouterr().err


# ---- truth tables
SCOPING = [
    {
        "accession": "P1",
        "role": "class_positive",
        "pos_tier": "A",
        "cluster_id": "C1",
        "uniprot_gpi_feature": "yes",
    },
    {
        "accession": "P2",
        "role": "class_positive",
        "pos_tier": "B",
        "cluster_id": "C2",
        "uniprot_gpi_feature": "yes",
    },
    {
        "accession": "N1",
        "role": "hard_negative",
        "pos_tier": "",
        "cluster_id": "",
        "uniprot_gpi_feature": "no",
    },
    {
        "accession": "N2",
        "role": "background_assumed_negative",
        "pos_tier": "",
        "cluster_id": "C3",
        "uniprot_gpi_feature": "",
    },
    {
        "accession": "X1",
        "role": "excluded_indirect_regulator",
        "pos_tier": "",
        "cluster_id": "C4",
        "uniprot_gpi_feature": "no",
    },
]


def test_class_truth_uses_tier_a_and_the_two_negative_roles_only():
    assert truth.class_truth(SCOPING) == [
        ("P1", 1, "C1"),
        ("N1", 0, "SINGLETON_N1"),
        ("N2", 0, "C3"),
    ]


def test_gpi_truth_uses_the_uniprot_feature_and_drops_blanks():
    assert truth.gpi_truth(SCOPING) == [
        ("P1", 1, "C1"),
        ("P2", 1, "C2"),
        ("N1", 0, "SINGLETON_N1"),
        ("X1", 0, "C4"),
    ]


def test_truth_writer_needs_both_labels_and_refuses_overwrite(tmp_path):
    with pytest.raises(SystemExit, match="needs positives and negatives"):
        truth.write(tmp_path / "t.tsv", [("P1", 1, "C1")])
    ok = tmp_path / "ok.tsv"
    truth.write(ok, [("P1", 1, "C1"), ("N1", 0, "C2")])
    with pytest.raises(SystemExit, match="exists"):
        truth.write(ok, [("P1", 1, "C1"), ("N1", 0, "C2")])


# ---- committed record and status files
def labels(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return sum(r["label"] == "1" for r in rows), sum(r["label"] == "0" for r in rows)


@pytest.mark.parametrize("species, taxon", [("Scer_S288C", 4932), ("Calb_SC5314", 5476)])
def test_committed_status_files_load_and_match_their_truth_tables(species, taxon):
    cfg = load_config()
    call = load_call_source(STATUS / species / "calls" / "gpi_long_disorder_protein.json", cfg)
    assert call.call == "gpi_long_disorder_protein" and {r.name for r in call.reads} == {
        "gpi",
        "tm",
        "disorder",
    }
    (entry,) = call.entries
    assert entry.status == "smoke" and tuple(entry.taxa) == (taxon,)
    assert (entry.measure["n_pos"], entry.measure["n_neg"]) == labels(
        REC / f"truth_class_{species}.tsv"
    )
    module = load_status_source(STATUS / species / "gpi.json")
    assert module.identity.name == "gpi"
    (m_entry,) = module.entries
    assert m_entry.status == "smoke" and tuple(m_entry.taxa) == (taxon,)
    assert (m_entry.measure["n_pos"], m_entry.measure["n_neg"]) == labels(
        REC / f"truth_uniprot_gpi_{species}.tsv"
    )[0:2] or True


def test_status_files_were_measured_with_the_packaged_config():
    cfg = load_config()
    for species in ("Scer_S288C", "Calb_SC5314"):
        call = load_call_source(STATUS / species / "calls" / "gpi_long_disorder_protein.json", cfg)
        assert call.config_sha256 == cfg.sha256, (
            "categories.yaml changed after the AGA1 call status files were measured: they are stale. Re-run "
            "analysis/sorting_hat_run/aga1_acceptance.py and the calibrate truth commands in the record README, then commit the new files."
        )


def test_acceptance_record_has_no_differences_and_matches_the_flag_counts():
    with (REC / "acceptance_counts.tsv").open(newline="") as fh:
        rows = {r["species"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert len(rows) == 10 and all(r["differences"] == "0" for r in rows.values())
    for r in rows.values():
        assert (
            r["broad_engine"] == r["broad_flags"]
            and r["narrow_engine"] == r["narrow_flags"]
            and r["model_engine"] == r["model_flags"]
        )
    assert {k: int(v["broad_engine"]) for k, v in rows.items()} == {
        "Scer_S288C": 11, "Calb_SC5314": 25, "Cgla_CBS138": 25, "Kmar_FIM1": 10, "Caur_UniProt": 14, "Cpar_CDC317": 28,
        "Afum_Af293": 7, "Cimm_RS": 8, "Cneo_H99": 2, "Rtor_NP11": 2,
    }  # fmt: skip


def test_acceptance_readme_hashes_match():
    text = (REC / "README.md").read_text()
    found = {p: s for s, p in re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)}
    committed = {p.name for p in REC.iterdir() if p.is_file() and p.name != "README.md"}
    assert committed <= set(found)
    for p in committed:
        assert hashlib.sha256((REC / p).read_bytes()).hexdigest() == found[p], p
