#!/usr/bin/env python
"""Build the master per-species adhesion summary table for the kingdom survey."""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from join import build_species_table, find_missing_results, load_samples_taxonomy  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"
INPUT_DIR = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/input")
SAMPLES_CSV = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv")
OUT_DIR = Path(__file__).resolve().parent / "tables"

MATCHED_FIELDS = [
    "locustag",
    "asmid",
    "species_name",
    "phylum",
    "subphylum",
    "class",
    "order",
    "family",
    "genus",
    "total_proteins",
    "adhesion_count",
    "adhesion_fraction",
    "mean_adhesion_prob",
    "median_adhesion_prob",
]
UNMATCHED_FIELDS = ["stem", "reason"]
MISMATCHED_FIELDS = ["stem", "fai_locustag", "result_locustag"]
MISSING_FIELDS = ["stem", "fai_path"]


def _write_csv(path: Path, rows: list, fieldnames: list) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    taxonomy = load_samples_taxonomy(SAMPLES_CSV)
    matched, unmatched, mismatched = build_species_table(RESULTS_DIR, INPUT_DIR, taxonomy)
    missing = find_missing_results(RESULTS_DIR, INPUT_DIR)

    _write_csv(OUT_DIR / "species_adhesion_summary.csv", matched, MATCHED_FIELDS)
    _write_csv(OUT_DIR / "unmatched_species.csv", unmatched, UNMATCHED_FIELDS)
    _write_csv(OUT_DIR / "mismatched_locustag.csv", mismatched, MISMATCHED_FIELDS)
    _write_csv(OUT_DIR / "missing_results.csv", missing, MISSING_FIELDS)

    n_results = sum(1 for _ in RESULTS_DIR.glob("*.adhesion_predict.csv"))
    print(f"Result files found: {n_results}")
    print(f"Matched: {len(matched)}  Unmatched: {len(unmatched)}  Mismatched: {len(mismatched)}")
    print(f"Missing results (fai present, no result file): {len(missing)}")

    assert (
        len(matched) + len(unmatched) + len(mismatched) == n_results
    ), "Reconciliation failed: matched+unmatched+mismatched != result file count"


if __name__ == "__main__":
    main()
