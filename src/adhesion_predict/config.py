"""Configuration constants for adhesionPredict."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
# MODELS_DIR = PROJECT_ROOT / "models"
DATA_DIR = Path.cwd() / "data"
MODELS_DIR = Path.cwd() / "models"

POSITIVE_DIR = DATA_DIR / "positive"
NEGATIVE_DIR = DATA_DIR / "negative"

DEFAULT_MODEL = "esm2_t6_8M_UR50D"
DEFAULT_TEST_SIZE = 0.2
DEFAULT_BATCH_SIZE = 8

ESM2_MODELS = {
    "esm2_t6_8M_UR50D": {"params": "8M", "layers": 6},
    "esm2_t6_35M_UR50D": {"params": "35M", "layers": 6},
    "esm2_t6_150M_UR50D": {"params": "150M", "layers": 6},
    "esm2_t12_35M_UR50D": {"params": "35M", "layers": 12},
}

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
