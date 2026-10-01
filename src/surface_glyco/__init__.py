"""surface_glyco - ESM-2 based scorer for fungal cell-surface glycoproteins."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("surface_glyco")
except PackageNotFoundError:  # a checkout that was not pip-installed under this name
    __version__ = "0+unknown"

from surface_glyco.io import find_fasta_files, process_fasta_file
from surface_glyco.model import load_model, predict, save_model, train_classifier

__all__ = [
    "__version__",
    "find_fasta_files",
    "process_fasta_file",
    "load_model",
    "predict",
    "save_model",
    "train_classifier",
]
