#!/usr/bin/env python3
"""Exploratory: add a sequence-level periodicity rule to the repeat call and measure it on the truth tables.

Current call: either detector gives a period, with region coverage >= 0.25 and >= 2.5 copies.
Candidate extra rule: detector 14 `rep_z_seq` >= Z (strong periodicity of the whole sequence, even when no
repeat region could be extracted). NOTHING IS CHANGED in the module. The truth tables are small and
partly tuned; an extra rule needs a held-out truth set and the owner's decision.

Usage: zseq_rule_scan.py --work WORKDIR_ROOT --truth DIR --prefix truth_v2 --out FILE
"""

import argparse
import csv
import sys


def read(path):
    return {r["protein"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def current(row, cov=0.25, copies=2.5):
    return (
        bool(row)
        and int(row["rep_period"]) > 0
        and float(row["rep_coverage"]) >= cov
        and float(row["rep_n_copies"]) >= copies
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--work", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--prefix", default="truth_v2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = []
    for prot in ("Scer_S288C", "Calb_SC5314"):
        truth = list(
            csv.DictReader(open(f"{a.truth}/{a.prefix}.{prot}.annotated.tsv"), delimiter="\t")
        )
        r02 = read(f"{a.work}/{prot}/raw/repeats/repeat02.tsv")
        r14 = read(f"{a.work}/{prot}/raw/repeats/repeat14.tsv")
        for z in (None, 6.0, 5.5, 5.0, 4.5, 4.0, 3.5):
            tp = fp = fn = tn = 0
            gained, falsepos = [], []
            for t in truth:
                base = current(r02.get(t["id"])) or current(r14.get(t["id"]))
                row14 = r14.get(t["id"])  # proteins under the detectors' minimum length have no row
                extra = z is not None and bool(row14) and float(row14["rep_z_seq"]) >= z
                called = base or extra
                y = t["label"] == "1"
                tp += y and called
                fn += y and not called
                fp += (not y) and called
                tn += (not y) and not called
                if extra and not base:
                    (gained if y else falsepos).append(t["gene"] or t["accession"])
            out.append(
                {
                    "proteome": prot,
                    "z_seq_min": "none" if z is None else z,
                    "positives": tp + fn,
                    "sensitivity": f"{tp / (tp + fn):.3f}",
                    "negatives": fp + tn,
                    "false_positives": fp,
                    "specificity": f"{tn / (fp + tn):.3f}",
                    "gained_positives": ";".join(gained),
                    "added_false_positives": ";".join(falsepos),
                }
            )
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
