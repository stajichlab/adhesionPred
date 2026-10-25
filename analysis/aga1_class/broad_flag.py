"""Candidate broad flag "Ser/Thr-rich GPI cell-wall glycoprotein" from the AGA1-class record table (stdlib only).

The flag is `yes` only if all three parts are yes: SignalP 6 signal peptide, PredGPI GPI-anchor call, and a disordered run of at
least --min-disordered-run residues. With --stalk-table (owner decision 2026-10-24) the run must also hold a Ser/Thr stalk window
(`stalk_run_len` of stalk_features.py: longest AIUPred run >= 0.5 holding a 30-residue window that is at least 50% Ser/Thr);
with --no-stalk-rule it is the longest AIUPred run (`disorder_longest_run`). One of the two must be named.
Three-valued: any part `no` gives `no`; otherwise any part `NA` (PredGPI never scores a protein of 40 aa or fewer) gives `NA`;
otherwise `yes`.

Optional, off by default and exploratory: --tail-rescue-min-hydro X (needs --stalk-table). A protein with a signal peptide and a
PredGPI `no` whose C-terminal tail hydrophobicity (`ctail_hydro`) is at least X gets its GPI part set to yes and `gpi_rescued`
yes. The rescue is a rule found by reading the PredGPI misses of the same 850 proteins, so it is not independent evidence.

--min-disordered-run is the one tunable; the owner set 200 on 2026-10-24 (100 is the pre-declared sensitivity value). The
value is written into every row (`broad_min_run`), so a table cannot be read without its setting. The record must hold the
expected disorder tool (--disorder-tool, default aiupred) or the script stops.

Usage: python3 -I broad_flag.py --record RECORD.tsv --out FLAGS.tsv (--stalk-table STALK.tsv | --no-stalk-rule)
         [--min-disordered-run 200] [--tail-rescue-min-hydro X] [--summary-out SUMMARY.tsv]
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
    "accession", "role_group", "signal_peptide_pred", "gpi_pred", "disorder_longest_run", "stalk_rule", "stalk_run_len",
    "ctail_hydro", "tail_rescue_min_hydro", "gpi_rescued", "broad_min_run", "disordered_run_ok", "broad_flag",
]  # fmt: skip


class FlagError(Exception):
    pass


def _yes_no_na(row, col):
    v = row[col]
    if v not in ("yes", "no", "NA"):
        raise FlagError(f"{row['accession']}: {col} is {v!r}, expected yes, no or NA")
    return v


def _run_ok(acc, value, label, min_run):
    if value == "NA":
        return "NA"
    if not value.isdigit():
        raise FlagError(f"{acc}: {label} is {value!r}, expected a non-negative integer or NA")
    return "yes" if int(value) >= min_run else "no"


def combine(parts):
    """Three-valued AND: no beats NA, NA beats yes."""
    if "no" in parts:
        return "no"
    return "NA" if "NA" in parts else "yes"


def _stalk_rule_text(stalk):
    if stalk is None:
        return "none"
    keys = {(r["stalk_window"], r["stalk_min_frac"], r["stalk_threshold"]) for r in stalk.values()}
    if len(keys) != 1:
        raise FlagError("the stalk table mixes several window / fraction / threshold settings")
    w, f, t = keys.pop()
    return f"window={w};min_st_frac={f};disorder_threshold={t}"


def flag_rows(rows, min_run, tool="aiupred", stalk=None, rescue=None):
    """`stalk` is None (no stalk rule) or {accession: row of stalk_features.py}; `rescue` is a minimum ctail_hydro or None."""
    if min_run < 1:
        raise FlagError(f"--min-disordered-run must be at least 1, not {min_run}")
    if rescue is not None and stalk is None:
        raise FlagError("--tail-rescue-min-hydro needs a stalk table (it holds ctail_hydro)")
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
    if stalk is not None and set(stalk) != {r["accession"] for r in rows}:
        raise FlagError("the stalk table and the record table do not hold the same accessions")
    rule = _stalk_rule_text(stalk)
    out = []
    for r in rows:
        acc = r["accession"]
        s = stalk[acc] if stalk is not None else None
        run = _run_ok(acc, s["stalk_run_len"] if s else r["disorder_longest_run"],
                      "stalk_run_len" if s else "disorder_longest_run", min_run)  # fmt: skip
        sp, gpi = _yes_no_na(r, "signal_peptide_pred"), _yes_no_na(r, "gpi_pred")
        rescued = "no"
        if (
            rescue is not None
            and sp == "yes"
            and gpi == "no"
            and s["ctail_hydro"] != "NA"
            and float(s["ctail_hydro"]) >= rescue
        ):
            gpi, rescued = "yes", "yes"
        out.append({
            "accession": acc, "role_group": r["role_group"], "signal_peptide_pred": sp, "gpi_pred": r["gpi_pred"],
            "disorder_longest_run": r["disorder_longest_run"], "stalk_rule": rule,
            "stalk_run_len": s["stalk_run_len"] if s else "NA", "ctail_hydro": s["ctail_hydro"] if s else "NA",
            "tail_rescue_min_hydro": "NA" if rescue is None else rescue, "gpi_rescued": rescued, "broad_min_run": min_run,
            "disordered_run_ok": run, "broad_flag": combine([sp, gpi, run]),
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
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument(
        "--stalk-table",
        help="stalk_features.py output; the disordered run must hold a Ser/Thr stalk window",
    )
    g.add_argument(
        "--no-stalk-rule",
        action="store_true",
        help="use the longest AIUPred run without the stalk rule",
    )
    ap.add_argument("--tail-rescue-min-hydro", type=float, default=None)
    ap.add_argument("--summary-out")
    a = ap.parse_args(argv)
    for p in (a.out, a.summary_out):
        if p and Path(p).exists():
            print(f"FATAL: {p} exists; nothing is overwritten", file=sys.stderr)
            return 1
    try:
        stalk = None
        if a.stalk_table:
            with open(a.stalk_table, newline="") as fh:
                stalk = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
        with open(a.record, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        flagged = flag_rows(
            rows, a.min_disordered_run, a.disorder_tool, stalk, a.tail_rescue_min_hydro
        )
    except (FlagError, KeyError, OSError) as e:
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
