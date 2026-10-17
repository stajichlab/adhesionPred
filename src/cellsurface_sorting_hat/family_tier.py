"""Evidence tier of a Pfam family, derived from a committed measurement file.

The tier (``untested``, ``smoke``, ``estimated``) describes the family model's own evidence. The call
status (``status.py``) describes the measurement of a call. Neither implies the other.

``data/sorting_hat/family_measure.tsv`` has one row per family with the counts of a measurement.
Nobody edits the ``specificity_tier`` column of ``family_table.tsv`` by hand: this module derives it
and CI fails when the table disagrees (``tests/cellsurface_sorting_hat/test_family_tier.py``).

    python -m cellsurface_sorting_hat.family_tier --check      # exit 1 when the table disagrees
    python -m cellsurface_sorting_hat.family_tier --write      # rewrite the tier column only
    python -m cellsurface_sorting_hat.family_tier --convert FAMILY_TESTS.tsv --out family_measure.tsv

The derivation calls ``status.status_from_measure``; it does not copy the thresholds.
"""

import argparse
import csv
import sys
from pathlib import Path

from cellsurface_sorting_hat.status import (
    ESTIMATED,
    SMOKE,
    UNVALIDATED,
    status_from_measure,
    weakest,
)

UNTESTED = "untested"
TIERS = (UNTESTED, SMOKE, ESTIMATED)
# descriptive_shared_sets: one shared set of positives of all classes and one set of hard negatives,
# not made for this family. It gives counts and a coverage, not a sensitivity: capped at smoke.
# per_family_reviewed: the family's own members and near-neighbour negatives, reviewed by a person.
MEASURE_BASES = ("descriptive_shared_sets", "per_family_reviewed")
DATA = Path(__file__).resolve().parents[2] / "data" / "sorting_hat"
TABLE_PATH = DATA / "family_table.tsv"
MEASURE_PATH = DATA / "family_measure.tsv"
MEASURE_COLUMNS = [
    "pfam_acc",
    "name",
    "measure_basis",
    "denominator",
    "n_pos",  # the family's own independent members; empty when the set is shared
    "pos_hits",  # shared positives that hit the family (coverage, never sensitivity)
    "pos_denominator",
    "coverage",
    "n_neg",
    "n_clusters_pos",
    "n_clusters_neg",
    "neg_hits",
    "neg_hits_n2",
    "neg_clusters_hit",
    "negative_basis",
    "sensitivity",
    "sens_lo",
    "sens_hi",
    "specificity",
    "spec_lo",
    "spec_hi",
    "n_pos_in_reference",
    "n_pos_unknown",
    "source",
    "date",
    "pfam_release",
    "hmmer_version",
    "notes",
]
TIER_COLUMN = "specificity_tier"


def _number(row, key, kind=int):
    text = (row.get(key) or "").strip()
    return kind(text) if text else None


def _rate(row, value_key, lo_key, hi_key):
    value, lo, hi = (_number(row, k, float) for k in (value_key, lo_key, hi_key))
    if value is None or lo is None or hi is None:
        return None
    return {"value": value, "lo": lo, "hi": hi}


def derive_tier(row):
    """The tier of one family from its measurement row (a dict of strings), or ``untested`` for None.

    Positives are the family's own members (``n_pos``); a shared set has none, so its hits
    (``pos_hits``) stand in. No positive or no negative gives ``untested``. A row on shared sets is
    capped at ``smoke``. A row whose positives are mostly ``in_reference`` or ``unknown`` is capped at
    ``smoke`` (the same cap as other statuses). A missing cluster count does not meet the floor
    (``status_from_measure``). Negatives that are not ``curated`` cannot reach ``estimated``."""
    if row is None:
        return UNTESTED
    basis = (row.get("measure_basis") or "").strip()
    if basis not in MEASURE_BASES:
        raise ValueError(
            f"{row.get('pfam_acc', '?')}: measure_basis must be one of {MEASURE_BASES} (got {basis!r})"
        )
    n_pos = _number(row, "n_pos")
    if n_pos is None:
        n_pos = _number(row, "pos_hits") or 0
    n_neg = _number(row, "n_neg") or 0
    if n_pos < 1 or n_neg < 1:
        return UNTESTED
    measure = {"n_pos": n_pos, "n_neg": n_neg}
    for key in ("n_clusters_pos", "n_clusters_neg"):
        value = _number(row, key)
        if value is not None:
            measure[key] = value
    negatives = (row.get("negative_basis") or "").strip()
    if negatives:
        measure["negative_basis"] = {"all": {"basis": negatives, "n": n_neg}}
    # no interval was measured on a shared set: the widest interval says so
    coverage = _number(row, "coverage", float) or 0.0
    measure["sensitivity"] = _rate(row, "sensitivity", "sens_lo", "sens_hi") or {
        "value": coverage,
        "lo": 0.0,
        "hi": 1.0,
    }
    specificity = _rate(row, "specificity", "spec_lo", "spec_hi")
    if specificity:
        measure["specificity"] = specificity
    status = status_from_measure(measure)
    if status == UNVALIDATED:
        return UNTESTED
    if basis != "per_family_reviewed":
        status = weakest([status, SMOKE])
    leaky = (_number(row, "n_pos_in_reference") or 0) + (_number(row, "n_pos_unknown") or 0)
    if 2 * leaky > n_pos:
        status = weakest([status, SMOKE])
    return status


def read_measure(path):
    """``{pfam_acc: row}`` of the measurement file. A missing column or a duplicate is refused."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in MEASURE_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing column(s) {missing}")
        rows = {}
        for r in reader:
            if r["pfam_acc"] in rows:
                raise ValueError(f"{path}:{reader.line_num}: duplicate {r['pfam_acc']}")
            rows[r["pfam_acc"]] = r
    return rows


def _read_table(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "pfam_acc" not in (reader.fieldnames or []):
            raise ValueError(f"{path}: no pfam_acc column")
        return list(reader.fieldnames), list(reader)


def derive_tiers(table_path=TABLE_PATH, measure_path=MEASURE_PATH):
    """``{pfam_acc: tier}`` for every family of the table."""
    measures = read_measure(measure_path)
    _, rows = _read_table(table_path)
    return {r["pfam_acc"]: derive_tier(measures.get(r["pfam_acc"])) for r in rows}


def check(table_path=TABLE_PATH, measure_path=MEASURE_PATH):
    """A list of problems; empty when the table agrees with the derived tiers."""
    fields, rows = _read_table(table_path)
    measures = read_measure(measure_path)
    problems = []
    if TIER_COLUMN not in fields:
        return [f"{table_path}: no {TIER_COLUMN} column"]
    known = {r["pfam_acc"] for r in rows}
    for acc in sorted(set(measures) - known):
        problems.append(f"{measure_path}: {acc} is not in the family table")
    for r in rows:
        derived = derive_tier(measures.get(r["pfam_acc"]))
        if r[TIER_COLUMN] != derived:
            problems.append(
                f"{r['pfam_acc']} {r['name']}: {TIER_COLUMN} is {r[TIER_COLUMN]!r} but the "
                f"measurement derives {derived!r}"
            )
    return problems


def write_tiers(table_path=TABLE_PATH, measure_path=MEASURE_PATH):
    """Rewrite the tier column of the table. Every other column keeps its value."""
    fields, rows = _read_table(table_path)
    tiers = derive_tiers(table_path, measure_path)
    if TIER_COLUMN not in fields:
        fields.append(TIER_COLUMN)
    for r in rows:
        r[TIER_COLUMN] = tiers[r["pfam_acc"]]
    with open(table_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ---- conversion of the descriptive family test (2026-10-09) ---------------------------------------

FAMILY_TESTS_REPORT = "docs/reports/2026-10-09-family-tests.md"
FAMILY_TESTS_DATE = "2026-10-09"
# Set sizes of that test (see the report): 137 N1 hard negatives in 76 clusters. The positives are
# not clustered.
N1_ROWS, N1_CLUSTERS = 137, 76


def convert_family_tests(path):
    """Measurement rows (dicts) from ``family_tests.tsv`` of the 2026-10-09 descriptive test.

    Every row has ``measure_basis=descriptive_shared_sets``. ``pos_E1_hit`` (``9/54``) becomes
    ``pos_hits`` and ``pos_denominator`` and a ``coverage``: the 54 E1 adhesins are of all classes,
    so a hit count is the share of a mixed panel, not the sensitivity of the family. The family's
    own members are not known, so ``n_pos`` and ``sensitivity`` stay empty."""
    out = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            hits, denominator = (int(x) for x in r["pos_E1_hit"].split("/"))
            n1, n1_clusters, n2 = r["neg_N1_hit"], r["neg_N1_clusters_hit"], r["neg_N2_hit"]
            row = dict.fromkeys(MEASURE_COLUMNS, "")
            row.update(
                {
                    "pfam_acc": r["pfam_acc"],
                    "name": r["name"],
                    "measure_basis": "descriptive_shared_sets",
                    "denominator": "e1_adhesins_all_classes",
                    "pos_hits": hits,
                    "pos_denominator": denominator,
                    "coverage": f"{hits / denominator:.3f}",
                    "n_neg": N1_ROWS,
                    "n_clusters_neg": N1_CLUSTERS,
                    "neg_hits": n1,
                    "neg_hits_n2": n2,
                    "neg_clusters_hit": n1_clusters,
                    "negative_basis": "curated",
                    "source": FAMILY_TESTS_REPORT,
                    "date": FAMILY_TESTS_DATE,
                    "pfam_release": "38.2",
                    "hmmer_version": "not recorded",
                    "notes": (
                        f"E1 hits {hits} of {denominator} adhesins of all classes (coverage, not "
                        f"sensitivity); N1 hits {n1} in {n1_clusters} of {N1_CLUSTERS} clusters; "
                        f"N2 hits {n2}; lists are agent-curated and not owner-reviewed; E1 "
                        "positives are not clustered"
                    ),
                }
            )
            out.append(row)
    return out


def write_measure(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, MEASURE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat.family_tier", description=__doc__)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="exit 1 when the table disagrees")
    mode.add_argument("--write", action="store_true", help="rewrite the tier column of the table")
    mode.add_argument("--convert", metavar="FAMILY_TESTS_TSV", help="write the measurement file")
    ap.add_argument("--table", default=str(TABLE_PATH))
    ap.add_argument("--measure", default=str(MEASURE_PATH))
    ap.add_argument("--out", help="with --convert: the measurement file to write")
    args = ap.parse_args(argv)
    try:
        if args.convert:
            write_measure(args.out or args.measure, convert_family_tests(args.convert))
            return 0
        if args.write:
            write_tiers(args.table, args.measure)
            return 0
        problems = check(args.table, args.measure)
    except (ValueError, OSError) as err:
        print(f"family_tier: error: {err}", file=sys.stderr)
        return 2
    for problem in problems:
        print(problem, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
