"""Repeat detector tables -> modules ``repeat02`` and ``repeat14`` (call = repeat per 03_repeat_surface_candidates.py)."""

import csv
import subprocess
import sys
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

DEFAULT_MIN_COVERAGE = 0.25
DEFAULT_MIN_COPIES = 2.5
DETECTOR_MIN_LEN = 80  # the --min-len default of both detector scripts
SCRIPTS = {"repeat02": "02_repeat_profile.py", "repeat14": "14_repeat_detect_general.py"}


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
            out[r["protein"]] = {
                "period": int(float(r["rep_period"] or 0)),
                "copies": float(r["rep_n_copies"] or 0),
                "coverage": float(r["rep_coverage"] or 0),
            }
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
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed and len(p.sequence) < min_len:
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
    """Run an existing detector script from ``analysis/cocci_repeats`` on one FASTA."""
    script = Path(repo_root) / "analysis" / "cocci_repeats" / SCRIPTS[module]
    cmd = [python, str(script), str(fasta), "--out", str(out_tsv), *extra]
    subprocess.run(cmd, check=True)
    return cmd
