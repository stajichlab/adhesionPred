"""Configuration constants for adhesionPredict."""

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent


PACKAGED_MODELS_DIR = Path(__file__).parent / "models"


def get_models_dir():
    """Directory holding trained models.

    The packaged directory is used unless SURFACE_GLYCO_MODELS_DIR is set. A ./models directory
    in the working directory is deliberately ignored so a stray folder cannot replace a model.
    """
    override = os.environ.get("SURFACE_GLYCO_MODELS_DIR")
    return Path(override) if override else PACKAGED_MODELS_DIR


def model_filename(esm_model):
    """File name of the model trained on embeddings from esm_model."""
    return f"surface_glyco_model_{esm_model}.pkl"


DATA_DIR = Path.cwd() / "data"

POSITIVE_DIR = DATA_DIR / "positive"
NEGATIVE_DIR = DATA_DIR / "negative"

DEFAULT_MODEL = "esm2_t6_8M_UR50D"
DEFAULT_TEST_SIZE = 0.2
DEFAULT_BATCH_SIZE = 8

FASTA_EXTENSIONS = [
    ".aa",
    ".faa",
    ".pep",
    ".fa",
    ".fas",
    ".fasta",
    ".aa.gz",
    ".faa.gz",
    ".pep.gz",
    ".fa.gz",
    ".fasta.gz",
]
