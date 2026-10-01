"""Read prediction CSVs written by surface_glyco_predict."""

import csv
from pathlib import Path

SCORE_COLUMN = "surface_glycoprotein_score"
CALLED_LABEL = "surface_glycoprotein"


def read_all(path):
    """Return (id, label, score) for every row."""
    path = Path(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or SCORE_COLUMN not in reader.fieldnames:
            raise ValueError(
                f"{path}: expected a {SCORE_COLUMN!r} column (header: {reader.fieldnames})"
            )
        return [(row["id"], row["prediction"], float(row[SCORE_COLUMN])) for row in reader]


def read_called(path):
    """Return (id, score) for rows called as surface glycoproteins."""
    return [(i, score) for i, label, score in read_all(path) if label == CALLED_LABEL]
