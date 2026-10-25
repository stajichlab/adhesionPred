"""Classifier scores of ``analysis/aga1_class/feature_classifier.py apply`` -> module ``aga1_model`` (at least 2 of 3 models at their threshold).

The engine and this wrapper never import scikit-learn: the wrapper reads the score tables (one per model, columns ``accession`` and ``score_<model>``) and
the signed thresholds table, and writes the scores, the thresholds used, the clade of the row, whether that row is calibrated, the number of models that
call the protein and the ``call``. The identity of the module (``artefact_digest``) covers the threshold rows used and the sha256 of every score table, so a
changed threshold or score table makes status files stale.
"""

import csv
import hashlib
import json
import math

from cellsurface_sorting_hat.modules.aga1_thresholds import (
    MODELS,
    select_thresholds,
)
from cellsurface_sorting_hat.modules.base import invalid_row, sha256_file

MODULE = "aga1_model"
MIN_MODELS = 2
COLUMNS = [
    "call", "score_logreg", "score_rf", "score_gboost", "thr_logreg", "thr_rf", "thr_gboost", "n_models_call", "threshold_clade", "calibrated",
]  # fmt: skip


class ScoreFormatError(ValueError):
    """A score table does not have the expected shape."""


def parse_scores(path, model):
    """Return ``{protein ID: score}`` from a TSV with columns ``accession`` and ``score_<model>``."""
    col = f"score_{model}"
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if (
            reader.fieldnames is None
            or "accession" not in reader.fieldnames
            or col not in reader.fieldnames
        ):
            raise ScoreFormatError(
                f"{path}: needs the columns 'accession' and {col!r}, has {reader.fieldnames}"
            )
        for n, r in enumerate(reader, 2):
            pid = (r["accession"] or "").strip()
            if not pid:
                raise ScoreFormatError(f"{path}:{n}: empty accession")
            if pid in out:
                raise ScoreFormatError(f"{path}:{n}: duplicate accession {pid!r}")
            try:
                s = float(r[col])
            except (TypeError, ValueError):
                raise ScoreFormatError(f"{path}:{n}: {col} is not a number: {r[col]!r}") from None
            if not math.isfinite(s) or not 0.0 <= s <= 1.0:
                raise ScoreFormatError(f"{path}:{n}: {col} {s} is outside [0, 1]")
            out[pid] = s
    if not out:
        raise ScoreFormatError(f"{path}: no rows")
    return out


def model_rows(proteins, taxa, scores, threshold_rows, lineage):
    """``scores`` is ``{model: {id: score}}``; ``taxa`` is ``{id: taxon}``. Returns (rows, rows used for the identity)."""
    rows, used = [], set()
    cache = {}
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        have = [scores[m].get(p.id) for m in MODELS]
        if any(v is None for v in have):
            rows.append({"id": p.id, "state": "error"})
            continue
        taxon = taxa[p.id]
        if taxon not in cache:
            cache[taxon] = select_thresholds(threshold_rows, lineage.ancestors(taxon))
        chosen = cache[taxon]
        for r in chosen.values():
            used.add((r.clade_taxon, r.model, r.threshold, r.calibrated))
        called = sum(1 for m, v in zip(MODELS, have, strict=True) if v >= chosen[m].threshold)
        clade = next(iter(chosen.values()))
        row = {"id": p.id, "state": "ok", "call": "called" if called >= MIN_MODELS else "not_called", "n_models_call": called,
               "threshold_clade": clade.clade_taxon, "calibrated": clade.calibrated}  # fmt: skip
        for m, v in zip(MODELS, have, strict=True):
            row[f"score_{m}"] = f"{v:.5f}"
            row[f"thr_{m}"] = f"{chosen[m].threshold:g}"
        rows.append(row)
    return rows, sorted(used)


def identity_digest(used, score_paths):
    """Hash of the threshold rows used and of the score tables (so a change of either changes the module identity)."""
    payload = {
        "thresholds": [list(map(str, u)) for u in used],
        "scores": {m: sha256_file(p) for m, p in sorted(score_paths.items())},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
