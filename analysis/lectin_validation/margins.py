#!/usr/bin/env python3.12
"""Score minus GA for every positive hit of the 15 M1/M2 models (how close to the cutoff each hit is).

Usage: margins.py TRUTH_TSV SCAN_HITS_TSV MODELS_HMM OUT_TSV
"""

import csv
import sys

truth_path, hits_path, hmm_path, out_path = sys.argv[1:5]
ga = {}
acc = None
for line in open(hmm_path):
    if line.startswith("ACC"):
        acc = line.split()[1].split(".")[0]
    if line.startswith("GA"):
        ga[acc] = float(line.split()[1])
truth = {r["accession"]: r for r in csv.DictReader(open(truth_path), delimiter="\t")}
rows = []
for r in csv.DictReader(open(hits_path), delimiter="\t"):
    t = truth[r["accession"]]
    if t["label"] != "positive":
        continue
    s = float(r["score"])
    rows.append(
        {
            "accession": r["accession"],
            "label_class": t["label_class"],
            "model": r["model"],
            "score": s,
            "ga": ga[r["model"]],
            "score_minus_ga": round(s - ga[r["model"]], 1),
        }
    )
rows.sort(key=lambda x: x["score_minus_ga"])
with open(out_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(len(rows), "positive hits; smallest margin", rows[0]["score_minus_ga"])
