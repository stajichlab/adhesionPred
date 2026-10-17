"""Pool S proteins with a Pfam 38.2 PF07938 GA hit, and those without a Pfam 37.x PF07938 label in the pool table."""

import csv
import gzip
import os
import sys

pools, sel, hits_tsv = sys.argv[1:4]
grp = {r["name"]: r["group"] for r in csv.DictReader(open(sel), delimiter="\t")}
hits = {}
for r in csv.DictReader(open(hits_tsv), delimiter="\t"):
    if r["model"] == "Fungal_lectin":
        hits.setdefault(r["proteome"], {})[r["protein"]] = float(r["seq_score"])
inpool, nopfam = [], []
for name, h in hits.items():
    with gzip.open(os.path.join(pools, name + ".pools.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["protein"] in h:
                inpool.append((name, r["protein"], h[r["protein"]], r["pool"]))
                if "PF07938" not in (r["pfam"] or "").split(";"):
                    nopfam.append((grp[name], name, r["protein"], h[r["protein"]], r["length"]))
print("Pfam 38.2 GA PF07938 hits that are in pool S:", len(inpool))
print("of these, without a 37.x PF07938 label:", len(nopfam))
for x in sorted(nopfam, key=lambda x: -x[3]):
    print(*x)
