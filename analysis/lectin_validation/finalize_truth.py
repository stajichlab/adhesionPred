#!/usr/bin/env python3.12
"""Add cluster ids and leakage to truth_all.tsv.

Leakage rule (COMMON-RULES rule 6, applied to the Pfam seed sequences of the 20 models):
  in_reference  the sequence has 30% identity or more (target coverage 0.5 or more, E <= 1e-3,
                mmseqs2 easy-search -s 7.5) to a seed sequence of at least one of the 20 models
  none          no such hit. This does NOT exclude remote homology that mmseqs2 at 30% misses.
Columns added: cluster30 (representative accession), seed_hits (MODEL:pident;...), leakage.
Usage: finalize_truth.py TRUTH_IN CLUSTERS_TSV SEED_HITS_M8 TRUTH_OUT CLUSTERS_OUT
"""

import csv
import sys
from collections import defaultdict

truth_in, clu, hits, truth_out, clu_out = sys.argv[1:6]
rep = {}
for line in open(clu):
    r, m = line.rstrip("\n").split("\t")
    rep[m] = r
seed = defaultdict(dict)
for line in open(hits):
    q, t, pid, aln, qc, tc, ev, bits = line.rstrip("\n").split("\t")
    model = t.split("|")[0]
    p = float(pid)
    if p < 1.5:  # mmseqs reports a fraction
        p *= 100
    if p >= 30.0 and float(tc) >= 0.5 and float(ev) <= 1e-3:
        if p > seed[q].get(model, 0):
            seed[q][model] = p
rows = list(csv.DictReader(open(truth_in), delimiter="\t"))
for r in rows:
    r["cluster30"] = rep[r["accession"]]
    sh = seed.get(r["accession"], {})
    r["seed_hits"] = ";".join(f"{m}:{p:.1f}" for m, p in sorted(sh.items()))
    r["leakage"] = "in_reference" if sh else "none"
cols = list(rows[0].keys())
with open(truth_out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
with open(clu_out, "w") as fh:
    for r in rows:
        fh.write(r["accession"] + "\t" + r["cluster30"] + "\n")
print(len(rows), "rows;", sum(r["leakage"] == "in_reference" for r in rows), "in_reference")
