"""Review 3: compare an independent PF07938 rescan with the committed review 2 hit table.

Usage: python compare_rescan.py WORK_DIR COMMITTED_HITS_TSV
"""

import csv
import os
import sys

d, hits = sys.argv[1:3]
sample = [line.split("\t")[0] for line in open(os.path.join(d, "sample.txt"))]
com = {}
for r in csv.DictReader(open(hits), delimiter="\t"):
    if r["model"] == "Fungal_lectin" and r["proteome"] in sample:
        com[(r["proteome"], r["protein"])] = float(r["seq_score"])
new = {}
for s in sample:
    for line in open(os.path.join(d, "tbl", s + ".tbl")):
        if line.startswith("#"):
            continue
        p = line.split()
        new[(s, p[0])] = float(p[5])
gen_new = {k[0] for k in new}
gen_com = {k[0] for k in com}
print("sample proteomes", len(sample))
print("proteins: rescan", len(new), "committed", len(com), "same set", set(new) == set(com))
print("max score difference", max((abs(new[k] - com[k]) for k in set(new) & set(com)), default=0))
print("proteomes with hit: rescan", len(gen_new), "committed", len(gen_com))
print("only rescan", sorted(set(new) - set(com)))
print("only committed", sorted(set(com) - set(new)))
