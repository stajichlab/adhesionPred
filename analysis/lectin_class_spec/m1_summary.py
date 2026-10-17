#!/usr/bin/env python3.12
"""Summarise hmmsearch domtbl files made by m1_scan.sh.

Usage: m1_summary.py SCANDIR OUT_HITS_TSV
Reads SCANDIR/models.txt (accession:set) and every SCANDIR/*.domtbl.
Writes one row per (FASTA, protein, model) and prints counts of distinct proteins per FASTA and set.
"""

import collections
import glob
import os
import sys

scan, out = sys.argv[1:3]
sets = {}
for line in open(os.path.join(scan, "models.txt")):
    acc, s = line.strip().split(":")
    sets[acc] = s
hits = collections.defaultdict(set)  # (fasta, set) -> proteins
rows = []
for f in sorted(glob.glob(os.path.join(scan, "*.domtbl"))):
    fasta = os.path.basename(f)[: -len(".domtbl")]
    seen = set()
    for line in open(f):
        if line.startswith("#"):
            continue
        p = line.split()
        target, query_name, query_acc = p[0], p[3], p[4].split(".")[0]
        seq_score, dom_score = float(p[7]), float(p[13])
        key = (target, query_acc)
        s = sets.get(query_acc, "?")
        hits[(fasta, s)].add(target)
        if key in seen:
            continue
        seen.add(key)
        rows.append((fasta, target, query_acc, query_name, s, seq_score, dom_score))
with open(out, "w") as fh:
    fh.write("fasta\tprotein\tmodel_acc\tmodel_name\tset\tseq_score\tdomain_score\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
per_model = collections.Counter((r[0], r[2]) for r in rows)
names = {r[2]: r[3] for r in rows}
with open(out.replace(".tsv", "_per_model.tsv"), "w") as fh:
    fh.write("model_acc\tmodel_name\tset\t" + "\t".join(sorted({r[0] for r in rows})) + "\n")
    for acc in sets:
        fh.write(
            f"{acc}\t{names.get(acc, '')}\t{sets[acc]}\t"
            + "\t".join(str(per_model.get((fa, acc), 0)) for fa in sorted({r[0] for r in rows}))
            + "\n"
        )
fastas = sorted({k[0] for k in hits})
print("fasta\t" + "\t".join(["M1", "M2", "M3", "own"]))
for fa in fastas:
    print(fa + "\t" + "\t".join(str(len(hits.get((fa, s), ()))) for s in ["M1", "M2", "M3", "own"]))
