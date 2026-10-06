"""Repeat detector tables -> modules ``repeat02`` and ``repeat14`` (call = repeat per 03_repeat_surface_candidates.py)."""

import csv
import math
import subprocess
import sys
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

DEFAULT_MIN_COVERAGE = 0.25
DEFAULT_MIN_COPIES = 2.5
DETECTOR_MIN_LEN = 80  # the --min-len default of both detector scripts
SCRIPTS = {"repeat02": "02_repeat_profile.py", "repeat14": "14_repeat_detect_general.py"}


def _number(text, column):
    """A finite number that is not negative; an empty value is 0 (the detectors write it for no repeat).
    A short row gives ``None`` here and is refused."""
    if text is None:
        raise ValueError(f"column {column!r} is missing (short row)")
    value = float(text or 0)
    if not math.isfinite(value):
        raise ValueError(f"column {column!r} value {text!r} is not finite")
    if value < 0:
        raise ValueError(f"column {column!r} value {text!r} is negative")
    return value


def parse_repeat_table(path):
    """Return ``{protein ID: {"period": int, "copies": float, "coverage": float}}``."""
    out = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"protein", "rep_period", "rep_n_copies", "rep_coverage"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
        for n, r in enumerate(reader, 2):
            if r["protein"] in out:
                raise ValueError(f"{path}:{n}: duplicate protein {r['protein']!r}")
            try:
                period = _number(r["rep_period"], "rep_period")
                copies = _number(r["rep_n_copies"], "rep_n_copies")
                coverage = _number(r["rep_coverage"], "rep_coverage")
                out[r["protein"]] = {"period": int(period), "copies": copies, "coverage": coverage}
            except (ValueError, OverflowError) as exc:
                raise ValueError(
                    f"{path}:{n}: bad numeric field for {r['protein']!r}: {exc}"
                ) from exc
    return out


def repeat_rows(
    proteins,
    parsed,
    min_coverage=DEFAULT_MIN_COVERAGE,
    min_copies=DEFAULT_MIN_COPIES,
    min_len=DETECTOR_MIN_LEN,
):
    """The detectors skip proteins shorter than ``min_len`` (default 80). Such a protein cannot hold
    a repeat array that the detectors can see, so it is ``not_called`` with ``period`` 0."""
    stray = sorted(set(parsed) - {p.id for p in proteins})
    if stray:
        raise ValueError(
            f"{len(stray)} repeat-table ID(s) are not in the FASTA, for example {stray[0]!r}"
        )
    rows = []
    for p in proteins:
        # the detectors strip trailing "*" before the length check
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed and len(p.sequence.rstrip("*")) < min_len:
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "not_called",
                    "period": 0,
                    "copies": "0.00",
                    "coverage": "0.000",
                }
            )
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            r = parsed[p.id]
            called = r["period"] > 0 and r["coverage"] >= min_coverage and r["copies"] >= min_copies
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if called else "not_called",
                    "period": r["period"],
                    "copies": f"{r['copies']:.2f}",
                    "coverage": f"{r['coverage']:.3f}",
                }
            )
    return rows


COLUMNS = ["call", "period", "copies", "coverage"]


def run_detector(repo_root, module, fasta, out_tsv, python=sys.executable, extra=()):
    """Run an existing detector script from ``analysis/cocci_repeats`` on one FASTA.

    ``--min-len`` is passed from ``DETECTOR_MIN_LEN`` so the wrapper and the detector agree.
    """
    if module not in SCRIPTS:
        raise ValueError(f"unknown repeat module {module!r}; expected one of {sorted(SCRIPTS)}")
    script = Path(repo_root) / "analysis" / "cocci_repeats" / SCRIPTS[module]
    cmd = [
        python,
        str(script),
        str(fasta),
        "--out",
        str(out_tsv),
        "--min-len",
        str(DETECTOR_MIN_LEN),
        *extra,
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"{module} failed on {fasta} (exit {exc.returncode}). Both detector scripts fail with"
            f" IndexError when no protein is at least {DETECTOR_MIN_LEN} aa (empty output)."
        ) from exc
    if not Path(out_tsv).exists():
        raise RuntimeError(f"{module} ran on {fasta} but wrote no output file {out_tsv}")
    return cmd
