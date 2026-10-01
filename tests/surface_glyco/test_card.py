"""Card contents and settings resolution."""

import pytest

from surface_glyco.card import (
    MAX_RESIDUES,
    POOLING_RESIDUE_MEAN,
    ModelCardError,
    new_card,
    resolve_embedding_settings,
)
from surface_glyco.model import load_model_card

DEFAULT = "esm2_t6_8M_UR50D"


def _card(**over):
    card = new_card(DEFAULT, 6, POOLING_RESIDUE_MEAN, n_positive=10, n_negative=20)
    card.update(over)
    return card


def test_new_card_has_v2_fields_and_claims_nothing_unmeasured():
    card = _card()
    assert card["card_version"] == 2
    assert (card["pooling"], card["max_residues"]) == (POOLING_RESIDUE_MEAN, MAX_RESIDUES)
    assert card["truncation"] == "truncate_at_max_residues"
    assert card["validation"] == []
    assert card["threshold"] == {"value": 0.5, "chosen_on": "default (uncalibrated)"}


def test_card_settings_are_used_when_cli_model_is_omitted():
    s = resolve_embedding_settings(_card(esm_model="esm2_t12_35M_UR50D"), None, DEFAULT)
    assert (s.esm_model, s.repr_layer, s.pooling) == ("esm2_t12_35M_UR50D", 6, POOLING_RESIDUE_MEAN)


def test_matching_cli_model_is_accepted():
    assert resolve_embedding_settings(_card(), DEFAULT, DEFAULT).esm_model == DEFAULT


@pytest.mark.parametrize(
    "card, field",
    [
        (_card(esm_model="esm2_t12_35M_UR50D"), "esm_model"),  # with cli DEFAULT below
        (_card(pooling="all_tokens_legacy"), "pooling"),
        (_card(pooling="max"), "pooling"),
        (_card(max_residues=2000), "max_residues"),
        (_card(repr_layer="six"), "repr_layer"),
    ],
)
def test_bad_cards_are_refused_naming_the_field(card, field):
    with pytest.raises(ModelCardError, match=field):
        resolve_embedding_settings(card, DEFAULT, DEFAULT)


def test_card_without_esm_model_is_refused():
    card = _card()
    del card["esm_model"]
    with pytest.raises(ModelCardError, match="esm_model"):
        resolve_embedding_settings(card, None, DEFAULT)


def test_missing_card_is_refused():
    with pytest.raises(ModelCardError, match="card"):
        resolve_embedding_settings(None, None, DEFAULT)


def test_malformed_card_json_names_the_file(tmp_path):
    model = tmp_path / "m.pkl"
    model.write_bytes(b"x")
    (tmp_path / "m.json").write_text("{not json")
    with pytest.raises(ModelCardError, match="m.json"):
        load_model_card(model)
    (tmp_path / "m.json").unlink()
    assert load_model_card(model) is None
