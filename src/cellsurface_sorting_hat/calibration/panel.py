"""Report-only check of calls against a panel of proteins with expected values."""

import csv
import gzip
from collections import Counter

EXPECTED = ("called", "not_called", "known_miss", "excluded")
LEAKY = ("in_reference", "tuned", "partial")
CALL_COLUMNS = ("protein", "call", "variant", "value")
PANEL_COLUMNS = ("protein", "call", "variant", "expected")


def _need(fieldnames, columns, path):
    missing = [c for c in columns if c not in (fieldnames or [])]
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}")


def read_calls(path):
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        _need(reader.fieldnames, CALL_COLUMNS, path)
        return {(r["protein"], r["call"], r["variant"]): r["value"] for r in reader}


def panel_check(calls_long, panel_tsv):
    """Compare ``calls_long`` with ``panel_tsv`` (columns protein, call, variant, expected, source; optional column tuning:
    ``in_reference``, ``tuned`` or ``partial`` for a protein that the module, its reference set or its
    cutoffs have seen).

    ``known_miss`` and ``excluded`` rows are listed and not counted as agreement or disagreement.
    A tuning value that is not empty and not one of ``LEAKY`` is refused: it must not be read as
    "not seen". Returns ``(rows, summary)``. Nothing here is an acceptance gate and no status changes.
    """
    calls = read_calls(calls_long)
    rows, summary = [], Counter()
    with open(panel_tsv, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        _need(reader.fieldnames, PANEL_COLUMNS, panel_tsv)
        for r in reader:
            where = f"{panel_tsv}:{reader.line_num}"
            if r["expected"] not in EXPECTED:
                raise ValueError(f"{where}: {r['protein']}: expected must be one of {EXPECTED}")
            tuning = (r.get("tuning") or "").strip()
            if tuning and tuning not in LEAKY:
                raise ValueError(
                    f"{where}: tuning must be empty or one of {LEAKY} (got {tuning!r})"
                )
            observed = calls.get((r["protein"], r["call"], r["variant"]), "missing")
            if r["expected"] in ("known_miss", "excluded"):
                verdict = r["expected"]
            elif tuning in LEAKY:
                verdict = "excluded_leakage"  # seen in tuning or part of the reference set
            elif observed == "missing":
                verdict = "not_in_run"  # the protein is not in the proteome that was run
            else:
                verdict = "agree" if observed == r["expected"] else "disagree"
            summary[verdict] += 1
            rows.append({**r, "observed": observed, "verdict": verdict})
    return rows, dict(summary)
