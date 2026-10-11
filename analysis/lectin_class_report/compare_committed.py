#!/usr/bin/env python3.12
"""Compare the hit table of lectin_report.py with the committed whole-proteome table.

Usage: /usr/bin/python3.12 compare_committed.py NEW_HITS_TSV COMMITTED_HITS_TSV
The committed table is analysis/lectin_class_spec/review2/whole_proteome_hits.tsv.
Only proteomes that appear in NEW_HITS_TSV or in the list of names (--names FILE, one name
per line) are compared. Use --names so that proteomes without any hit are not lost.
Models in the committed table: Fungal_lectin, FB_lectin, LDL, Gal-bind_lectin, CCL2-like.
"""

import argparse
import csv

MODEL_NAMES = {
    "PF07938": "Fungal_lectin",
    "PF07367": "FB_lectin",
    "PF21691": "LDL",
    "PF00337": "Gal-bind_lectin",
    "PF21595": "CCL2-like",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("new")
    ap.add_argument("committed")
    ap.add_argument("--names", help="file with the proteome names of the run")
    args = ap.parse_args()
    by_name = {v: k for k, v in MODEL_NAMES.items()}
    names = set()
    if args.names:
        names = {ln.strip().split("\t")[0] for ln in open(args.names) if ln.strip()}
    new = {}
    with open(args.new) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            names.add(r["proteome"]) if not args.names else None
            if r["model_acc"] in MODEL_NAMES:
                new[(r["model_acc"], r["proteome"], r["protein"])] = float(r["seq_score"])
    com = {}
    with open(args.committed) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["proteome"] in names and r["model"] in by_name:
                com[(by_name[r["model"]], r["proteome"], r["protein"])] = float(r["seq_score"])
    ok = True
    for acc in MODEL_NAMES:
        n = {k for k in new if k[0] == acc}
        c = {k for k in com if k[0] == acc}
        diff = max((abs(new[k] - com[k]) for k in n & c), default=0.0)
        same = n == c
        ok = ok and same and diff == 0.0
        print(
            "%s %s: new %d, committed %d, same set %s, max score difference %.1f, "
            "proteomes with hit new %d committed %d"
            % (
                acc,
                MODEL_NAMES[acc],
                len(n),
                len(c),
                same,
                diff,
                len({k[1] for k in n}),
                len({k[1] for k in c}),
            )
        )
        if not same:
            print("  only new:", sorted(n - c)[:10])
            print("  only committed:", sorted(c - n)[:10])
    print("proteomes compared: %d" % len(names))
    print("MATCH" if ok else "DIFFERENCE")


if __name__ == "__main__":
    main()
