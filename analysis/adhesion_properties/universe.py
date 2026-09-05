"""Build the adhesion-protein set and a matched per-species background
sample of non-adhesion proteins, for the protein-properties analysis.

Reuses kingdom_survey's join.py locustag-matching primitives (the .fai-
anchored join key, verified during that project's design review) rather
than re-deriving them.
"""

import csv
import random
import sys
from pathlib import Path

_KINGDOM_SURVEY_DIR = Path(__file__).resolve().parent.parent / "kingdom_survey"
sys.path.insert(0, str(_KINGDOM_SURVEY_DIR))
from join import first_id_from_fai, locustag_from_id  # noqa: E402

RESULT_SUFFIX = ".adhesion_predict.csv"
FAI_SUFFIX = ".proteins.fa.fai"
_TAXONOMY_FIELDS = ["phylum", "class", "order", "family", "genus"]


def load_taxonomy_by_locustag(summary_csv_path: Path) -> dict[str, dict]:
    """Load kingdom_survey's species_adhesion_summary.csv keyed by locustag,
    for attaching taxonomy to individual protein rows."""
    lookup = {}
    with open(summary_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lookup[row["locustag"]] = row
    return lookup


def read_result_ids(result_csv_path: Path) -> list[tuple[str, float]]:
    """Return [(protein_id, probability_adhesion), ...] for one result CSV."""
    ids = []
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ids.append((row["id"], float(row["probability_adhesion"])))
    return ids


def read_fai_ids(fai_path: Path) -> list[str]:
    """Return every protein id in a .fai index (first column)."""
    ids = []
    with open(fai_path) as fh:
        for line in fh:
            ids.append(line.split("\t", 1)[0])
    return ids


def sample_background_ids(
    all_ids: list[str], adhesion_ids: set[str], n_wanted: int, rng: random.Random
) -> list[str]:
    """Sample up to n_wanted protein ids from all_ids that are NOT in
    adhesion_ids. Returns all available non-adhesion ids if fewer than
    n_wanted exist."""
    candidates = [pid for pid in all_ids if pid not in adhesion_ids]
    if len(candidates) <= n_wanted:
        return candidates
    return rng.sample(candidates, n_wanted)


def build_protein_universe(
    results_dir: Path, input_dir: Path, summary_csv_path: Path, seed: int = 42
) -> tuple[list[dict], list[dict]]:
    """Build (adhesion_rows, background_rows) for every species with a
    result file whose locustag is in the taxonomy lookup. Species not
    found in the lookup (unmatched/mismatched per kingdom_survey) are
    skipped — this analysis only covers species kingdom_survey already
    successfully joined.
    """
    rng = random.Random(seed)
    taxonomy = load_taxonomy_by_locustag(summary_csv_path)

    adhesion_rows: list[dict] = []
    background_rows: list[dict] = []

    for result_path in sorted(results_dir.glob(f"*{RESULT_SUFFIX}")):
        stem = result_path.name[: -len(RESULT_SUFFIX)]
        fai_path = input_dir / f"{stem}{FAI_SUFFIX}"
        if not fai_path.exists():
            continue
        locustag = locustag_from_id(first_id_from_fai(fai_path))
        if locustag not in taxonomy:
            continue
        tax = taxonomy[locustag]
        tax_fields = {field: tax[field] for field in _TAXONOMY_FIELDS}

        result_ids = read_result_ids(result_path)
        adhesion_id_set = {pid for pid, _ in result_ids}
        for pid, prob in result_ids:
            adhesion_rows.append(
                {
                    "locustag": locustag,
                    **tax_fields,
                    "protein_id": pid,
                    "probability_adhesion": prob,
                }
            )

        all_ids = read_fai_ids(fai_path)
        background_ids = sample_background_ids(all_ids, adhesion_id_set, len(result_ids), rng)
        for pid in background_ids:
            background_rows.append(
                {
                    "locustag": locustag,
                    **tax_fields,
                    "protein_id": pid,
                    "probability_adhesion": None,
                }
            )

    return adhesion_rows, background_rows
