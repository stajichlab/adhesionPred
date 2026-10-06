"""Provenance record for a proteome FASTA that a run uses."""

import json
import math
from datetime import date
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.fasta import OK, FastaError, read_fasta
from cellsurface_sorting_hat.modules.base import sha256_file


class ProteomeError(ValueError):
    """The proteome does not look like the one that was expected."""


def _text(field, value):
    if not isinstance(value, str) or not value.strip():
        raise ProteomeError(f"{field} must be a non-empty string, got {value!r}")


def _bound(field, value):
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ProteomeError(f"{field} must be a number, got {value!r}")
    if not math.isfinite(value) or value < 0:
        raise ProteomeError(f"{field} must be finite and not negative, got {value!r}")


def write_provenance(path, name, source, fasta, retrieved, expected_min, expected_max):
    """Check the FASTA and write ``path`` (JSON). ``source`` is a URL or a file path; ``retrieved``
    is the date as ``YYYY-MM-DD`` text. The protein count must lie in
    ``[expected_min, expected_max]``. All checks run before ``path`` is written."""
    _text("name", name)
    _text("source", source)
    _text("retrieved", retrieved)
    try:
        date.fromisoformat(retrieved)
    except ValueError:
        raise ProteomeError(f"retrieved must be YYYY-MM-DD, got {retrieved!r}") from None
    _bound("expected_min", expected_min)
    _bound("expected_max", expected_max)
    if expected_min > expected_max:
        raise ProteomeError(f"expected_min {expected_min} is above expected_max {expected_max}")
    try:
        proteins = read_fasta(fasta)
    except (FastaError, FileNotFoundError) as err:
        raise ProteomeError(f"{fasta}: {err}") from err
    n = len(proteins)
    if not expected_min <= n <= expected_max:
        raise ProteomeError(f"{name}: {n} proteins, expected {expected_min} to {expected_max}")
    record = {
        "name": name,
        "source": source,
        "retrieved": retrieved,
        "fasta": str(Path(fasta).resolve()),
        "sha256": sha256_file(fasta),
        "n_proteins": n,
        "n_invalid": sum(p.state != OK for p in proteins),
        "n_unique_sequences": len({p.sha256 for p in proteins}),
        "first_ids": [p.id for p in proteins[:3]],
    }
    write_atomic(path, (json.dumps(record, indent=2, sort_keys=True) + "\n").encode())
    return record
