"""Tests for FASTA discovery, model cards, and prediction CSV output."""

import csv
from pathlib import Path

from sklearn.linear_model import LogisticRegression

from surface_glyco.features import extract_sequence_features
from surface_glyco.io import find_fasta_files, load_sequences_from_dir
from surface_glyco.model import load_model_card, save_model, save_model_card
from surface_glyco.scripts.predict import build_results, write_results


def test_find_fasta_files_returns_each_file_once(tmp_path):
    (tmp_path / "a.fasta").write_text(">x\nMK\n")
    (tmp_path / "b.pep").write_text(">y\nMK\n")
    files = find_fasta_files(tmp_path)
    assert [f.name for f in files] == ["a.fasta", "b.pep"]


def test_load_sequences_from_dir_reads_gzip(tmp_path):
    import gzip

    with gzip.open(tmp_path / "a.fa.gz", "wt") as f:
        f.write(">x\nMKJ*\n")
    seqs = load_sequences_from_dir(tmp_path)
    assert seqs == [{"id": "x", "sequence": "MKL"}]


def test_model_card_round_trip(tmp_path):
    model_path = tmp_path / "m.pkl"
    save_model(LogisticRegression(), model_path)
    assert load_model_card(model_path) is None
    save_model_card(model_path, {"esm_model": "esm2_t6_8M_UR50D", "repr_layer": 6})
    assert load_model_card(model_path)["repr_layer"] == 6


def test_write_results_header_and_quoting(tmp_path):
    out = tmp_path / "r.csv"
    write_results(
        [
            {
                "id": "sp|P1,x",
                "prediction": "surface_glycoprotein",
                "surface_glycoprotein_score": 0.91234,
            }
        ],
        out,
    )
    assert out.read_text().splitlines()[0] == "id,prediction,surface_glycoprotein_score"
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows == [
        {
            "id": "sp|P1,x",
            "prediction": "surface_glycoprotein",
            "surface_glycoprotein_score": "0.9123",
        }
    ]


def test_extract_sequence_features_treats_j_as_l():
    assert extract_sequence_features("JJ") == extract_sequence_features("LL")


def test_dedupe_sequences_keeps_first_id_per_label():
    from surface_glyco.scripts.train import dedupe_sequences

    seqs = [
        {"id": "a", "sequence": "MK", "label": 1},
        {"id": "b", "sequence": "MK", "label": 1},
        {"id": "c", "sequence": "MK", "label": 0},
    ]
    assert [s["id"] for s in dedupe_sequences(seqs)] == ["a", "c"]


def test_build_results_keeps_every_row_with_new_labels():
    rows = build_results(["a", "b"], [1, 0], [[0.1, 0.9], [0.8, 0.2]])
    assert [r["prediction"] for r in rows] == ["surface_glycoprotein", "other"]
    assert rows[1]["surface_glycoprotein_score"] == 0.2
    assert "probability_adhesion" not in rows[0]


def test_default_output_name_uses_new_suffix(tmp_path, monkeypatch):
    from surface_glyco.scripts.predict import default_output_path

    monkeypatch.chdir(tmp_path)
    assert default_output_path(Path("x/Scer.pep.fa")).name == "Scer.pep.surface_glyco.csv"
    (tmp_path / "proteomes").mkdir()
    assert default_output_path(tmp_path / "proteomes").name == "proteomes.surface_glyco.csv"
