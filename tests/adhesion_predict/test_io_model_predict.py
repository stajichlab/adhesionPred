"""Tests for FASTA discovery, model cards, and prediction CSV output."""

import csv

from adhesion_predict.features import extract_sequence_features
from adhesion_predict.io import find_fasta_files, load_sequences_from_dir
from adhesion_predict.model import load_model_card, save_model, save_model_card
from adhesion_predict.scripts.predict import build_results, write_results
from sklearn.linear_model import LogisticRegression


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


def test_write_results_quotes_ids_with_commas(tmp_path):
    out = tmp_path / "r.csv"
    write_results(
        [{"id": "sp|P1,x", "prediction": "Adhesion", "probability_adhesion": 0.91234}], out
    )
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows == [{"id": "sp|P1,x", "prediction": "Adhesion", "probability_adhesion": "0.9123"}]


def test_extract_sequence_features_treats_j_as_l():
    assert extract_sequence_features("JJ") == extract_sequence_features("LL")


def test_dedupe_sequences_keeps_first_id_per_label():
    from adhesion_predict.scripts.train import dedupe_sequences

    seqs = [
        {"id": "a", "sequence": "MK", "label": 1},
        {"id": "b", "sequence": "MK", "label": 1},
        {"id": "c", "sequence": "MK", "label": 0},
    ]
    assert [s["id"] for s in dedupe_sequences(seqs)] == ["a", "c"]


def test_build_results_keeps_non_adhesion_rows():
    rows = build_results(["a", "b"], [1, 0], [[0.1, 0.9], [0.8, 0.2]])
    assert [r["prediction"] for r in rows] == ["Adhesion", "Non-adhesion"]
    assert rows[1]["probability_adhesion"] == 0.2
