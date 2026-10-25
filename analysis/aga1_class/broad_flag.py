"""Candidate broad flag "Ser/Thr-rich GPI cell-wall glycoprotein" from the AGA1-class record table (stdlib only).

The flag is `yes` only if all three parts are yes: SignalP 6 signal peptide, PredGPI GPI-anchor call, and an AIUPred disordered
run (consecutive residues scoring at or above the record's threshold, 0.5) of at least --min-disordered-run residues.
Three-valued: any part `no` gives `no`; otherwise any part `NA` (PredGPI never scores a protein of 40 aa or fewer) gives `NA`;
otherwise `yes`. The Ser/Thr stalk rule is NOT part of this candidate; it is still an owner decision.

--min-disordered-run is the one tunable; the owner set 200 on 2026-10-24 (100 is the pre-declared sensitivity value). The
value is written into every row (`broad_min_run`), so a table cannot be read without its setting. The record must hold the
expected disorder tool (--disorder-tool, default aiupred) or the script stops.

Usage: python3 -I broad_flag.py --record RECORD.tsv --out FLAGS.tsv [--min-disordered-run 200] [--summary-out SUMMARY.tsv]
Refuses an existing output file. Fails loudly (FATAL:, exit 1) on a missing column, a mixed or wrong disorder tool, or a value
that is not yes/no/NA (parts) or a non-negative integer (run length).
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

DEFAULT_MIN_RUN = 200
PARTS = (
    "signal_peptide_pred",
    "gpi_pred",
    "disorder_longest_run",
    "disorder_tool",
    "role_group",
    "accession",
)
OUT_COLS = [
    "accession", "role_group", "signal_peptide_pred", "gpi_pred", "disorder_longest_run", "broad_min_run",
    "disordered_run_ok", "broad_flag",
]  # fmt: skip


class FlagError(Exception):
    pass


def _yes_no_na(row, col):
    v = row[col]
    if v not in ("yes", "no", "NA"):
        raise FlagError(f"{row['accession']}: {col} is {v!r}, expected yes, no or NA")
    return v


def _run_ok(row, min_run):
    v = row["disorder_longest_run"]
    if v == "NA":
        return "NA"
    if not v.isdigit():
        raise FlagError(
            f"{row['accession']}: disorder_longest_run is {v!r}, expected a non-negative integer or NA"
        )
    return "yes" if int(v) >= min_run else "no"


def combine(parts):
    """Three-valued AND: no beats NA, NA beats yes."""
    if "no" in parts:
        return "no"
    return "NA" if "NA" in parts else "yes"


def flag_rows(rows, min_run, tool="aiupred"):
    if min_run < 1:
        raise FlagError(f"--min-disordered-run must be at least 1, not {min_run}")
    if not rows:
        raise FlagError("the record table has no rows")
    missing = [c for c in PARTS if c not in rows[0]]
    if missing:
        raise FlagError(
            f"the record table lacks columns: {missing} (is this a table made with a disorder file?)"
        )
    tools = {r["disorder_tool"] for r in rows}
    if tools != {tool}:
        raise FlagError(f"disorder_tool is {sorted(tools)}, expected only {tool!r}")
    out = []
    for r in rows:
        run = _run_ok(r, min_run)
        sp, gpi = _yes_no_na(r, "signal_peptide_pred"), _yes_no_na(r, "gpi_pred")
        out.append({
            "accession": r["accession"], "role_group": r["role_group"], "signal_peptide_pred": sp, "gpi_pred": gpi,
            "disorder_longest_run": r["disorder_longest_run"], "broad_min_run": min_run, "disordered_run_ok": run,
            "broad_flag": combine([sp, gpi, run]),
        })  # fmt: skip
    return out


def summarise(flagged):
    n, c = Counter(), {}
    for r in flagged:
        n[r["role_group"]] += 1
        c.setdefault(r["role_group"], Counter())[r["broad_flag"]] += 1
    return [
        {"role_group": g, "n": n[g], "yes": c[g]["yes"], "no": c[g]["no"], "NA": c[g]["NA"]}
        for g in sorted(n)
    ]


def _write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--record", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-disordered-run", type=int, default=DEFAULT_MIN_RUN)
    ap.add_argument("--disorder-tool", default="aiupred")
    ap.add_argument("--summary-out")
    a = ap.parse_args(argv)
    for p in (a.out, a.summary_out):
        if p and Path(p).exists():
            print(f"FATAL: {p} exists; nothing is overwritten", file=sys.stderr)
            return 1
    try:
        with open(a.record, newline="") as fh:
            flagged = flag_rows(
                list(csv.DictReader(fh, delimiter="\t")), a.min_disordered_run, a.disorder_tool
            )
    except (FlagError, KeyError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    _write(a.out, flagged, OUT_COLS)
    summary = summarise(flagged)
    if a.summary_out:
        _write(a.summary_out, summary, ["role_group", "n", "yes", "no", "NA"])
    print(f"wrote {a.out}: {len(flagged)} proteins, minimum disordered run {a.min_disordered_run}")
    for s in summary:
        print("{role_group:45s} n={n:4d} yes={yes:3d} no={no:4d} NA={NA}".format(**s))
    return 0


if __name__ == "__main__":
    sys.exit(main())
