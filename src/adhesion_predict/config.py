"""Configuration constants for adhesionPredict."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent


# Model directory search hierarchy
def get_models_dir():
    """
    Get models directory with search hierarchy:
    1. Current working directory: ./models
    2. System installation: site-packages/adhesion_predict/models
    """
    # First priority: current directory
    current_dir_models = Path.cwd() / "models"
    if current_dir_models.exists():
        return current_dir_models

    # Second priority: system installation directory
    # Look for adhesion_predict package in site-packages
    package_root = Path(__file__).parent
    system_models = package_root / "models"
    if system_models.exists():
        return system_models

    # Fallback: return current directory (will be created if needed)
    return current_dir_models


DATA_DIR = Path.cwd() / "data"
MODELS_DIR = get_models_dir()

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
