#!/usr/bin/env python
"""Build the adhesion-protein set and matched background sample, writing
tables/protein_universe.csv (one row per protein, group=adhesion|background)."""

import csv
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from universe import build_protein_universe  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"
INPUT_DIR = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/input")
SUMMARY_CSV = REPO_ROOT / "analysis" / "kingdom_survey" / "tables" / "species_adhesion_summary.csv"
OUT_DIR = Path(__file__).resolve().parent / "tables"

FIELDS = [
    "protein_id",
    "group",
    "locustag",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "probability_adhesion",
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adhesion_rows, background_rows = build_protein_universe(
        RESULTS_DIR, INPUT_DIR, SUMMARY_CSV, seed=42
    )

    with open(OUT_DIR / "protein_universe.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        for row in adhesion_rows:
            writer.writerow({**row, "group": "adhesion"})
        for row in background_rows:
            writer.writerow({**row, "group": "background", "probability_adhesion": ""})

    print(f"Adhesion proteins: {len(adhesion_rows)}")
    print(f"Background proteins: {len(background_rows)}")

    summary = pd.read_csv(SUMMARY_CSV)
    expected = int(summary["adhesion_count"].sum())
    assert len(adhesion_rows) == expected, (
        f"Reconciliation failed: {len(adhesion_rows)} adhesion rows built, "
        f"expected {expected} from species_adhesion_summary.csv"
    )


if __name__ == "__main__":
    main()
