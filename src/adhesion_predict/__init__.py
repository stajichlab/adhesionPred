"""adhesionPredict - A bioinformatics package for predicting adhesion proteins."""

try:
    from importlib.metadata import version
except ImportError:
    # Python < 3.8
    from importlib_metadata import version

__version__ = version("adhesion_predict")

from adhesion_predict.io import find_fasta_files, process_fasta_file
from adhesion_predict.model import load_model, predict, save_model, train_classifier

__all__ = [
    "__version__",
    "find_fasta_files",
    "process_fasta_file",
    "load_model",
    "predict",
    "save_model",
    "train_classifier",
]
