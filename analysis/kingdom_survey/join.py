"""Join logic for building the kingdom-wide adhesion species table.

Join key: the LOCUSTAG prefix embedded in protein ids (e.g.
"F07B100A_000481-T1" -> "F07B100A"). Anchored on each species' .fai index
(verified 100%-reliable during design review), with the result CSV's own
first-row LOCUSTAG used as a cross-check that catches annotation mismatches.
"""

import csv
import statistics
from pathlib import Path
from typing import Optional

RESULT_SUFFIX = ".adhesion_predict.csv"
FAI_SUFFIX = ".proteins.fa.fai"


def locustag_from_id(seq_id: str) -> str:
    """Extract the LOCUSTAG prefix from a protein id like 'F07B100A_000481-T1'."""
    return seq_id.split("_", 1)[0]


def first_id_from_fai(fai_path: Path) -> str:
    """Return the id in the first column of the first line of a .fai index."""
    with open(fai_path) as fh:
        first_line = fh.readline()
    return first_line.split("\t", 1)[0]


def count_fai_lines(fai_path: Path) -> int:
    """Count total records (= total proteins) in a .fai index."""
    with open(fai_path) as fh:
        return sum(1 for _ in fh)


def first_id_from_result_csv(result_csv_path: Path) -> Optional[str]:
    """Return the id in the first data row of a result CSV, or None if header-only."""
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            return row["id"]
    return None


def count_result_rows(result_csv_path: Path) -> int:
    """Count adhesion-called proteins (data rows) in a result CSV."""
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        return sum(1 for _ in reader)


def probability_stats(result_csv_path: Path) -> tuple[float, float]:
    """Return (mean, median) of probability_adhesion over adhesion-called proteins."""
    probs = []
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            probs.append(float(row["probability_adhesion"]))
    if not probs:
        return (float("nan"), float("nan"))
    return (statistics.mean(probs), statistics.median(probs))


def load_samples_taxonomy(samples_csv_path: Path) -> dict[str, dict]:
    """Load samples.csv into a dict keyed by LOCUSTAG.

    Blank taxonomy fields are stored as None (not "") so downstream
    per-rank grouping can exclude a species from a rank rather than
    pooling it into a fake blank group. SUBCLASS is intentionally not
    loaded (>50% blank in samples.csv; dropped from analysis entirely
    per design review).
    """
    lookup: dict[str, dict] = {}
    with open(samples_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            locustag = row["LOCUSTAG"].strip()
            if not locustag:
                continue
            lookup[locustag] = {
                "asmid": row["ASMID"],
                "phylum": row["PHYLUM"].strip() or None,
                "subphylum": row["SUBPHYLUM"].strip() or None,
                "class": row["CLASS"].strip() or None,
                "order": row["ORDER"].strip() or None,
                "family": row["FAMILY"].strip() or None,
                "genus": row["GENUS"].strip() or None,
                "species_name": row["SPECIES"].strip() or None,
            }
    return lookup


def build_species_table(
    results_dir: Path, input_dir: Path, taxonomy: dict[str, dict]
) -> tuple[list[dict], list[dict], list[dict]]:
    """Build the master species table plus QC rows.

    Returns (matched, unmatched, mismatched), where each is a list of
    dict rows. matched + unmatched + mismatched always sums to the
    number of result files found in results_dir.
    """
    matched: list[dict] = []
    unmatched: list[dict] = []
    mismatched: list[dict] = []

    for result_path in sorted(results_dir.glob(f"*{RESULT_SUFFIX}")):
        stem = result_path.name[: -len(RESULT_SUFFIX)]
        fai_path = input_dir / f"{stem}{FAI_SUFFIX}"

        if not fai_path.exists():
            unmatched.append({"stem": stem, "reason": "missing_fai"})
            continue

        fai_locustag = locustag_from_id(first_id_from_fai(fai_path))
        if fai_locustag not in taxonomy:
            unmatched.append(
                {"stem": stem, "reason": f"locustag_{fai_locustag}_not_in_samples_csv"}
            )
            continue

        result_first_id = first_id_from_result_csv(result_path)
        result_locustag = locustag_from_id(result_first_id) if result_first_id else None
        if result_locustag is not None and result_locustag != fai_locustag:
            mismatched.append(
                {
                    "stem": stem,
                    "fai_locustag": fai_locustag,
                    "result_locustag": result_locustag,
                }
            )
            continue

        tax = taxonomy[fai_locustag]
        total_proteins = count_fai_lines(fai_path)
        adhesion_count = count_result_rows(result_path)
        mean_prob, median_prob = probability_stats(result_path)

        matched.append(
            {
                "locustag": fai_locustag,
                "asmid": tax["asmid"],
                "species_name": tax["species_name"],
                "phylum": tax["phylum"],
                "subphylum": tax["subphylum"],
                "class": tax["class"],
                "order": tax["order"],
                "family": tax["family"],
                "genus": tax["genus"],
                "total_proteins": total_proteins,
                "adhesion_count": adhesion_count,
                "adhesion_fraction": adhesion_count / total_proteins,
                "mean_adhesion_prob": mean_prob,
                "median_adhesion_prob": median_prob,
            }
        )

    return matched, unmatched, mismatched


def find_missing_results(results_dir: Path, input_dir: Path) -> list[dict]:
    """Species with a .fai in the input dir but no corresponding result file."""
    result_stems = {p.name[: -len(RESULT_SUFFIX)] for p in results_dir.glob(f"*{RESULT_SUFFIX}")}
    missing = []
    for fai_path in sorted(input_dir.glob(f"*{FAI_SUFFIX}")):
        stem = fai_path.name[: -len(FAI_SUFFIX)]
        if stem not in result_stems:
            missing.append({"stem": stem, "fai_path": str(fai_path)})
    return missing
