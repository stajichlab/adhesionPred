#!/usr/bin/env python3
"""Exploratory: what happens to the repeat call if detector 14 used a lower region cut?

Loads analysis/cocci_repeats/14_repeat_detect_general.py as a library, changes REGION_THRESHOLD in memory
only, runs detect() on the truth proteins, and applies the call rule of the sorting hat (period > 0,
coverage >= 0.25, copies >= 2.5; detector 14 only, so this is a lower bound on the combined call).
Nothing on disk changes. A different cut needs a held-out truth set and the owner's decision.

Usage: region_cut_scan.py --detector FILE --fasta PROTEOME=PATH ... --truth DIR --prefix truth_v2 --out FILE
"""

import argparse
import csv
import importlib.util
import sys


def load(path):
    spec = importlib.util.spec_from_file_location("detector14", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).rstrip("*")
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf).rstrip("*")
    return seqs


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--detector", required=True)
    ap.add_argument("--fasta", action="append", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--prefix", default="truth_v2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    det = load(a.detector)
    rows = []
    for spec in a.fasta:
        prot, path = spec.split("=", 1)
        seqs = read_fasta(path)
        truth = list(
            csv.DictReader(open(f"{a.truth}/{a.prefix}.{prot}.annotated.tsv"), delimiter="\t")
        )
        for cut in (0.5, 0.4, 0.3, 0.2):
            det.REGION_THRESHOLD = cut
            tp = fp = fn = tn = 0
            gained, falsepos = [], []
            for t in truth:
                seq = seqs[t["id"]]
                r = det.detect(seq) if len(seq) >= 80 else {}
                called = (
                    bool(r)
                    and r.get("period", 0) > 0
                    and r.get("coverage", 0) >= 0.25
                    and r.get("n_copies", 0) >= 2.5
                )
                y = t["label"] == "1"
                tp += y and called
                fn += y and not called
                fp += (not y) and called
                tn += (not y) and not called
                if called:
                    (gained if y else falsepos).append(t["gene"] or t["accession"])
            rows.append(
                {
                    "proteome": prot,
                    "region_threshold": cut,
                    "positives": tp + fn,
                    "sensitivity_detector14_only": f"{tp / (tp + fn):.3f}",
                    "negatives": fp + tn,
                    "false_positives": fp,
                    "specificity": f"{tn / (fp + tn):.3f}",
                    "called_positives": ";".join(gained),
                    "called_negatives": ";".join(falsepos),
                }
            )
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
