"""adhesionPredict - A bioinformatics package for predicting adhesion proteins."""

__version__ = "0.1.0"

from adhesion_predict.io import find_fasta_files, process_fasta_file
from adhesion_predict.model import load_model, predict, save_model, train_classifier

__all__ = [
    "find_fasta_files",
    "process_fasta_file",
    "load_model",
    "predict",
    "save_model",
    "train_classifier",
]
