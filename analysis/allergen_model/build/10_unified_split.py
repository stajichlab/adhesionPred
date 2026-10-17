#!/usr/bin/env python3.12
"""Cluster-level split over one shared clustering of all 916 sequences (200 positives, 200 length-matched weak
negatives, 516 family-matched negatives). Positive clusters are split 70/15/15 by protein count (seed 20261010);
every other cluster is assigned by seeded draw. No cluster is split across train, validation and test.
Input clustering: unified/clu_cluster.tsv (representative<TAB>member). Output: unified/unified_splits.tsv.
"""

import collections
import csv
import random
import sys

SEED = 20261010
rep = {}
for ln in open("unified/clu_cluster.tsv"):
    r, m = ln.rstrip("\n").split("\t")
    rep[m] = r
rows = []
for r in csv.DictReader(open("splits/splits.tsv"), delimiter="\t"):
    if r["class"] == "positive" or r["class"] == "weak_negative":
        rows.append((r["id"], r["class"], r["accession"]))
for r in csv.DictReader(open("family_matched/family_negatives.tsv"), delimiter="\t"):
    rows.append((f"F_{r['accession']}", "family_negative", r["accession"]))
missing = [i for i, _, _ in rows if i not in rep]
if missing:
    sys.exit(f"FATAL: {len(missing)} ids not in the unified clustering, e.g. {missing[:3]}")
pos_clusters = collections.Counter(rep[i] for i, c, _ in rows if c == "positive")
rng = random.Random(SEED)
order = sorted(pos_clusters)
rng.shuffle(order)
total = sum(pos_clusters.values())
split_of = {}
acc = 0
for c in order:
    f = acc / total
    split_of[c] = "train" if f < 0.70 else ("validation" if f < 0.85 else "test")
    acc += pos_clusters[c]
for i, _c, _ in rows:
    k = rep[i]
    if k not in split_of:
        x = rng.random()
        split_of[k] = "train" if x < 0.70 else ("validation" if x < 0.85 else "test")
with open("unified/unified_splits.tsv", "w") as fh:
    fh.write("id\tclass\taccession\tsplit\tcluster\n")
    for i, c, a in rows:
        fh.write(f"{i}\t{c}\t{a}\t{split_of[rep[i]]}\t{rep[i]}\n")
cnt = collections.Counter((s, c) for s, c in [(split_of[rep[i]], c) for i, c, _ in rows])
print("clusters", len({rep[i] for i, _, _ in rows}), dict(cnt))
