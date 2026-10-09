#!/usr/bin/env python3
"""Exploratory: sensitivity and specificity of the repeat call on the truth tables, by subset and cutoff.

NOT a tuning step and nothing here changes the module. The cutoffs of the repeat call are module
parameters (coverage and copies, modules/repeats.py). The truth tables are small and several
positives tuned the detectors, so a scan on them can only describe the trade-off. Applying a new
cutoff needs a larger independent truth set and the owner's decision.

Usage: threshold_scan.py --work WORKDIR_ROOT --truth DIR [--prefix truth_v2] --out FILE
Reads <work>/<proteome>/raw/repeats/repeat02.tsv and repeat14.tsv (detector tables) and
<work>/<proteome>/out_calls/calls.long.tsv.gz (for the R0 call).
"""

import argparse
import csv
import gzip
import itertools
import sys

SPECIES = {"Scer_S288C": "S. cerevisiae S288C", "Calb_SC5314": "C. albicans SC5314"}


def read(path):
    return {r["protein"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def sp_calls(path):
    out = {}
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["call"] == "signal_peptide_protein" and r["variant"] == "R0":
                out[r["protein"]] = r["value"] == "called"
    return out


def detector_call(row, cov, copies):
    if not row:
        return False
    return (
        int(row["rep_period"]) > 0
        and float(row["rep_coverage"]) >= cov
        and float(row["rep_n_copies"]) >= copies
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--work", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--prefix", default="truth", help="truth or truth_v2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for prot in SPECIES:
        truth = list(
            csv.DictReader(open(f"{a.truth}/{a.prefix}.{prot}.annotated.tsv"), delimiter="\t")
        )
        r02 = read(f"{a.work}/{prot}/raw/repeats/repeat02.tsv")
        r14 = read(f"{a.work}/{prot}/raw/repeats/repeat14.tsv")
        sp = sp_calls(f"{a.work}/{prot}/out_calls/calls.long.tsv.gz")
        subsets = {
            "all": lambda t: True,
            "secreted (R0 called)": lambda t, sp=sp: sp.get(t["id"], False),
            "untuned": lambda t: t["tuned_or_homolog"] == "no",
        }
        for (cov, copies), (sname, keep) in itertools.product(
            itertools.product([0.05, 0.10, 0.15, 0.20, 0.25], [2.5, 3.0, 4.0]), subsets.items()
        ):
            tp = fp = fn = tn = 0
            for t in truth:
                if not keep(t):
                    continue
                called = detector_call(r02.get(t["id"]), cov, copies) or detector_call(
                    r14.get(t["id"]), cov, copies
                )
                y = t["label"] == "1"
                tp += y and called
                fn += y and not called
                fp += (not y) and called
                tn += (not y) and not called
            rows.append(
                {
                    "proteome": prot,
                    "subset": sname,
                    "min_coverage": cov,
                    "min_copies": copies,
                    "positives": tp + fn,
                    "negatives": fp + tn,
                    "tp": tp,
                    "fp": fp,
                    "sensitivity": f"{tp / (tp + fn):.3f}" if tp + fn else "NA",
                    "specificity": f"{tn / (fp + tn):.3f}" if fp + tn else "NA",
                }
            )
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
