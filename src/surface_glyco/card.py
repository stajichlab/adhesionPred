"""Model card contents and the rules for turning a card into embedding settings.

This module must not import torch or esm: predict and evaluate validate a model before
loading any weights, and these rules are tested without a GPU stack.
"""

from dataclasses import dataclass

CARD_VERSION = 2
POOLING_RESIDUE_MEAN = "residue_mean"
KNOWN_POOLINGS = (POOLING_RESIDUE_MEAN,)
MAX_RESIDUES = 1022  # ESM-2 context is 1024 tokens including BOS and EOS
DEFAULT_REPR_LAYER = 6


class ModelCardError(Exception):
    """The model card does not allow scoring with the requested settings."""


@dataclass(frozen=True)
class EmbeddingSettings:
    esm_model: str
    repr_layer: int
    pooling: str
    max_residues: int


def new_card(esm_model, repr_layer, pooling, n_positive, n_negative, **extra):
    """Card for a newly trained model. Fields that were not measured are left empty."""
    card = {
        "card_version": CARD_VERSION,
        "esm_model": esm_model,
        "repr_layer": repr_layer,
        "pooling": pooling,
        "max_residues": MAX_RESIDUES,
        "truncation": "truncate_at_max_residues",
        "n_positive": n_positive,
        "n_negative": n_negative,
        "validation": [],
        "threshold": {"value": 0.5, "chosen_on": "default (uncalibrated)"},
    }
    card.update(extra)
    return card


def resolve_embedding_settings(card, cli_model_name, default_model):
    """Return the embedding settings a model must be scored with, or raise ModelCardError."""
    if card is None:
        raise ModelCardError(
            "card: this model has no card, so its embedding settings are unknown. "
            "Train a model with surface_glyco_train, which writes one."
        )
    esm_model = card.get("esm_model")
    if not esm_model:
        raise ModelCardError("esm_model: missing from the model card")
    if cli_model_name is not None and cli_model_name != esm_model:
        raise ModelCardError(
            f"esm_model: the model was trained on {esm_model} embeddings, "
            f"but --model-name is {cli_model_name}"
        )
    pooling = card.get("pooling")
    if pooling not in KNOWN_POOLINGS:
        raise ModelCardError(f"pooling: {pooling!r} is not supported; supported: {KNOWN_POOLINGS}")
    max_residues = card.get("max_residues", MAX_RESIDUES)
    if max_residues != MAX_RESIDUES:
        raise ModelCardError(
            f"max_residues: the model used {max_residues}, this build supports {MAX_RESIDUES}"
        )
    try:
        repr_layer = int(card.get("repr_layer", DEFAULT_REPR_LAYER))
    except (TypeError, ValueError):
        raise ModelCardError(f"repr_layer: {card.get('repr_layer')!r} is not an integer") from None
    return EmbeddingSettings(esm_model, repr_layer, pooling, max_residues)
