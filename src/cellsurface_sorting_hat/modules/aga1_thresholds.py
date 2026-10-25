"""The signed thresholds table of the AGA1 model call (``data/sorting_hat/aga1_model_thresholds.tsv``) and the choice of a row by lineage.

One row per (clade, model). A protein's taxon picks, for each model, the row of the DEEPEST clade that is the taxon itself or one of its ancestors
(``taxonomy.Lineage.ancestors``). Taxon 1 is allowed in THIS table (the "all other taxa" fallback; the status schema's minimum taxon of 2 does not apply
here). A row with ``calibrated=no`` says that its numbers were not measured for the clade. Adding a calibrated clade is adding signed rows.
"""

import csv
import datetime
import math
from dataclasses import dataclass

HEADER = [
    "clade_taxon",
    "clade_name",
    "model",
    "threshold",
    "basis",
    "calibrated",
    "signed_by",
    "signed_date",
    "evidence",
]
MODELS = ("logreg", "rf", "gboost")


class ThresholdError(ValueError):
    """The thresholds table or a choice from it is not valid."""


@dataclass(frozen=True)
class Row:
    clade_taxon: int
    clade_name: str
    model: str
    threshold: float
    calibrated: str
    basis: str
    signed_by: str
    signed_date: str
    evidence: str


def load_thresholds(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header != HEADER:
            raise ThresholdError(f"{path}: header must be {HEADER}, got {header}")
        rows, seen = [], set()
        for n, f in enumerate(reader, 2):
            where = f"{path}:{n}"
            if len(f) != len(HEADER):
                raise ThresholdError(f"{where}: expected {len(HEADER)} fields, got {len(f)}")
            d = dict(zip(HEADER, f, strict=True))
            try:
                taxon = int(d["clade_taxon"])
            except ValueError:
                raise ThresholdError(f"{where}: clade_taxon is not an integer") from None
            if taxon < 1:
                raise ThresholdError(f"{where}: clade_taxon must be at least 1")
            if d["model"] not in MODELS:
                raise ThresholdError(f"{where}: model {d['model']!r} is not one of {MODELS}")
            try:
                thr = float(d["threshold"])
            except ValueError:
                raise ThresholdError(f"{where}: threshold is not a number") from None
            if not math.isfinite(thr) or not 0 < thr < 1:
                raise ThresholdError(f"{where}: threshold must be in (0, 1), got {d['threshold']}")
            if d["calibrated"] not in ("yes", "no"):
                raise ThresholdError(f"{where}: calibrated must be yes or no")
            if not d["signed_by"].strip():
                raise ThresholdError(f"{where}: signed_by is empty (every row is signed)")
            try:
                datetime.date.fromisoformat(d["signed_date"])
            except ValueError:
                raise ThresholdError(f"{where}: signed_date must be YYYY-MM-DD") from None
            if not d["basis"].strip() or not d["evidence"].strip():
                raise ThresholdError(f"{where}: basis and evidence must not be empty")
            if (taxon, d["model"]) in seen:
                raise ThresholdError(
                    f"{where}: duplicate row for clade {taxon} and model {d['model']}"
                )
            seen.add((taxon, d["model"]))
            rows.append(
                Row(
                    taxon,
                    d["clade_name"],
                    d["model"],
                    thr,
                    d["calibrated"],
                    d["basis"],
                    d["signed_by"],
                    d["signed_date"],
                    d["evidence"],
                )
            )
    if not rows:
        raise ThresholdError(f"{path}: no rows")
    for clade in sorted({r.clade_taxon for r in rows}):
        missing = [m for m in MODELS if (clade, m) not in seen]
        if missing:
            raise ThresholdError(f"{path}: clade {clade} lacks rows for {missing}")
    return rows


def select_thresholds(rows, ancestors):
    """The row of each model for a protein whose taxon has the given ``ancestors`` (taxon first, root last).

    Returns ``{model: Row}``. Stops if a model has no row for any ancestor (the table has no fallback for taxon 1) or if the rows
    of one protein come from different clades (the three models of a clade are signed together).
    """
    depth = {t: i for i, t in enumerate(ancestors)}
    chosen = {}
    for m in MODELS:
        cands = [r for r in rows if r.model == m and r.clade_taxon in depth]
        if not cands:
            raise ThresholdError(
                f"no threshold row for model {m!r} and any of the taxa {ancestors[:3]}..."
            )
        chosen[m] = min(cands, key=lambda r: depth[r.clade_taxon])
    if len({r.clade_taxon for r in chosen.values()}) != 1:
        raise ThresholdError(
            "the thresholds of the three models come from different clades; sign them together"
        )
    return chosen
