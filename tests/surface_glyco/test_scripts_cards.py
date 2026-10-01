"""predict/evaluate/train honour the model card."""

import sys

import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression

from surface_glyco.card import DEFAULT_REPR_LAYER, POOLING_RESIDUE_MEAN, new_card
from surface_glyco.model import save_model, save_model_card
from surface_glyco.scripts import evaluate as evaluate_mod
from surface_glyco.scripts import predict as predict_mod
from surface_glyco.scripts import train as train_mod

ESM = "esm2_t6_8M_UR50D"


def _model(tmp_path, card=None, n_features=4):
    rng = np.random.default_rng(0)
    clf = LogisticRegression().fit(rng.normal(size=(20, n_features)), [0, 1] * 10)
    path = tmp_path / "m.pkl"
    save_model(clf, path)
    if card is not None:
        save_model_card(path, card)
    return path


@pytest.fixture
def fasta(tmp_path):
    p = tmp_path / "in.fa"
    p.write_text(">a\nMKT\n>b\n" + "M" * 1100 + "\n")
    return p


@pytest.fixture
def fake_embed(monkeypatch):
    calls = {}

    def fake(sequences, **kw):
        calls.update(kw)
        return np.ones((len(sequences), 4)), [s["id"] for s in sequences]

    monkeypatch.setattr(predict_mod, "get_esm_embeddings", fake)
    return calls


def test_card_settings_reach_the_embedder(tmp_path, fasta, fake_embed):
    out = tmp_path / "o.csv"
    card = new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1)
    predict_mod.main(fasta, _model(tmp_path, card), out, None, silent=True)
    assert (fake_embed["model_name"], fake_embed["repr_layer"]) == (ESM, 6)
    assert out.read_text().splitlines()[0] == "id,prediction,surface_glycoprotein_score"


@pytest.mark.parametrize(
    "card, cli_model, field",
    [
        (new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1), "esm2_t12_35M_UR50D", "esm_model"),
        (new_card(ESM, 6, "all_tokens_legacy", 1, 1), None, "pooling"),
        (None, None, "card"),
    ],
)
def test_refusals_exit_1_naming_the_field_before_embedding(
    tmp_path, fasta, fake_embed, capsys, card, cli_model, field
):
    with pytest.raises(SystemExit) as e:
        predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", cli_model)
    assert e.value.code == 1
    assert field in capsys.readouterr().err
    assert fake_embed == {}  # nothing was embedded


def test_truncated_sequences_are_reported(tmp_path, fasta, fake_embed, capsys):
    card = new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1)
    predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", None, silent=True)
    assert "1 of 2 sequences" in capsys.readouterr().out


def test_cli_passes_model_and_default_model_name(tmp_path, fasta, monkeypatch):
    got = {}
    monkeypatch.setattr(predict_mod, "main", lambda *a, **k: got.update(args=a, kwargs=k))
    model = _model(tmp_path, new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1))
    monkeypatch.setattr(
        sys, "argv", ["surface_glyco_predict", "--input", str(fasta), "--model", str(model)]
    )
    predict_mod.cli()
    assert got["args"][3] is None  # model_name left to the card


def test_dedupe_drops_sequences_present_in_both_classes():
    seqs = [
        {"id": "a", "sequence": "MK", "label": 1},
        {"id": "b", "sequence": "MK", "label": 0},
        {"id": "c", "sequence": "MKT", "label": 1},
        {"id": "d", "sequence": "MKT", "label": 1},
    ]
    assert [s["id"] for s in train_mod.dedupe_sequences(seqs)] == ["c"]


def test_empty_sequences_are_not_trained_on():
    seqs = [{"id": "a", "sequence": "", "label": 1}, {"id": "b", "sequence": "MK", "label": 0}]
    assert [s["id"] for s in train_mod.dedupe_sequences(seqs)] == ["b"]


def test_removal_counts_are_attributed_correctly_and_sum_to_the_input():
    def rec(i, seq, label):
        return {"id": i, "sequence": seq, "label": label}

    # empties in both classes must not count as conflicts; one empty must not count as a duplicate
    seqs = [
        rec("e1", "", 1),
        rec("e2", "", 0),
        rec("a", "MK", 1),
        rec("b", "MK", 0),
        rec("c", "MKT", 1),
        rec("d", "MKT", 1),
        rec("f", "MKTA", 0),
    ]
    unique, counts = train_mod.dedupe_with_counts(seqs)
    assert [s["id"] for s in unique] == ["c", "f"]
    assert counts == {"n_empty_removed": 2, "n_conflicting_removed": 2, "n_duplicates_removed": 1}
    assert len(unique) + sum(counts.values()) == len(seqs)


def test_end_to_end_with_real_embeddings(tmp_path):
    from surface_glyco.embeddings import get_esm_embeddings

    seqs = [
        {"id": f"s{i}", "sequence": "MKTAYIAKQRQISFVKSHFSRQ"[: 8 + i], "label": i % 2}
        for i in range(6)
    ]
    emb, _ = get_esm_embeddings(seqs, model_name=ESM)
    clf = LogisticRegression(max_iter=200).fit(emb, [s["label"] for s in seqs])
    model = tmp_path / "m.pkl"
    save_model(clf, model)
    save_model_card(model, new_card(ESM, 6, POOLING_RESIDUE_MEAN, 3, 3))
    fa = tmp_path / "in.fa"
    fa.write_text("".join(f">{s['id']}\n{s['sequence']}\n" for s in seqs))
    out = tmp_path / "o.csv"
    predict_mod.main(fa, model, out, None, silent=True)
    rows = out.read_text().splitlines()
    assert len(rows) == 7 and rows[0] == "id,prediction,surface_glycoprotein_score"


def test_trained_model_card_has_environment_and_counts(tmp_path, monkeypatch):
    # prepare_data and embeddings are stubbed (40 records: train_classifier runs 5-fold CV);
    # this pins what train.main writes to the card
    from surface_glyco.model import load_model_card

    counts = {"n_empty_removed": 1, "n_conflicting_removed": 2, "n_duplicates_removed": 3}
    monkeypatch.setattr(
        train_mod,
        "prepare_data",
        lambda p, n: (
            [{"id": f"s{i}", "sequence": "MK" * (i + 2), "label": i % 2} for i in range(40)],
            counts,
        ),
    )
    monkeypatch.setattr(
        train_mod,
        "get_esm_embeddings",
        lambda seqs, **kw: (
            np.random.default_rng(1).normal(size=(len(seqs), 4)),
            [s["id"] for s in seqs],
            list(range(len(seqs))),
        ),
    )
    out = tmp_path / "m.pkl"
    train_mod.main(tmp_path, tmp_path, out, ESM, 0.2)
    card = load_model_card(out)
    assert card["card_version"] == 2 and card["pooling"] == POOLING_RESIDUE_MEAN
    assert (
        card["n_empty_removed"],
        card["n_conflicting_removed"],
        card["n_duplicates_removed"],
    ) == (
        1,
        2,
        3,
    )
    assert set(card["environment"]) == {"numpy", "scikit_learn", "torch", "fair_esm"}


def test_train_passes_repr_layer_explicitly_and_card_matches(tmp_path, monkeypatch):
    from surface_glyco.model import load_model_card

    seen = {}
    monkeypatch.setattr(
        train_mod,
        "prepare_data",
        lambda p, n: (
            [{"id": f"s{i}", "sequence": "MK" * (i + 2), "label": i % 2} for i in range(40)],
            {"n_empty_removed": 0, "n_conflicting_removed": 0, "n_duplicates_removed": 0},
        ),
    )

    def fake(seqs, **kw):
        seen.update(kw)
        return (
            np.random.default_rng(1).normal(size=(len(seqs), 4)),
            [s["id"] for s in seqs],
            list(range(len(seqs))),
        )

    monkeypatch.setattr(train_mod, "get_esm_embeddings", fake)
    out = tmp_path / "m.pkl"
    train_mod.main(tmp_path, tmp_path, out, ESM, 0.2)
    assert seen["repr_layer"] == DEFAULT_REPR_LAYER
    assert load_model_card(out)["repr_layer"] == DEFAULT_REPR_LAYER


def test_predict_checks_the_card_before_unpickling(tmp_path, fasta, fake_embed, capsys):
    bad = tmp_path / "bad.pkl"
    bad.write_bytes(b"not a pickle")
    save_model_card(bad, new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1))
    with pytest.raises(SystemExit) as e:
        predict_mod.main(fasta, bad, tmp_path / "o.csv", "esm2_t12_35M_UR50D")
    assert e.value.code == 1
    assert "esm_model" in capsys.readouterr().err
    assert fake_embed == {}


@pytest.fixture
def eval_dirs(tmp_path):
    pos = tmp_path / "pos"
    neg = tmp_path / "neg"
    pos.mkdir()
    neg.mkdir()
    (pos / "p.fa").write_text(">p1\nMKTAYIAK\n>p2\nMKTAYIAKQR\n>p3\nMKTAYIAKQRQ\n")
    (neg / "n.fa").write_text(">n1\nGGSSGGSA\n>n2\nGGSSGGSAGG\n>n3\nGGSSGGSAGGT\n")
    return pos, neg


@pytest.fixture
def eval_embed(monkeypatch):
    calls = {}

    def fake(sequences, **kw):
        calls.update(kw)
        calls["n"] = len(sequences)
        rng = np.random.default_rng(2)
        emb = rng.normal(size=(len(sequences), 4))
        return emb, [s["id"] for s in sequences], list(range(len(sequences)))

    monkeypatch.setattr(evaluate_mod, "get_esm_embeddings", fake)
    return calls


def test_evaluate_card_settings_reach_the_embedder(tmp_path, eval_dirs, eval_embed):
    card = new_card(ESM, 6, POOLING_RESIDUE_MEAN, 3, 3)
    evaluate_mod.main(*eval_dirs, _model(tmp_path, card), None)
    assert (eval_embed["model_name"], eval_embed["repr_layer"]) == (ESM, 6)
    assert eval_embed["n"] == 6


@pytest.mark.parametrize(
    "card, cli_model, field",
    [
        (new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1), "esm2_t12_35M_UR50D", "esm_model"),
        (None, None, "card"),
    ],
)
def test_evaluate_refusals_exit_1_before_embedding(
    tmp_path, eval_dirs, eval_embed, capsys, card, cli_model, field
):
    with pytest.raises(SystemExit) as e:
        evaluate_mod.main(*eval_dirs, _model(tmp_path, card), cli_model)
    assert e.value.code == 1
    assert field in capsys.readouterr().err
    assert eval_embed == {}


def test_evaluate_cli_leaves_model_name_to_the_card(tmp_path, eval_dirs, monkeypatch):
    got = {}
    monkeypatch.setattr(evaluate_mod, "main", lambda *a, **k: got.update(args=a, kwargs=k))
    model = _model(tmp_path, new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "surface_glyco_evaluate",
            "--positive",
            str(eval_dirs[0]),
            "--negative",
            str(eval_dirs[1]),
            "--model",
            str(model),
        ],
    )
    evaluate_mod.cli()
    assert got["args"][3] is None
