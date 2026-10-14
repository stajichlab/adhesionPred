"""Review 3: whole-proteome PF07938 scores of the pool S PF07938 proteins (option (a) population).

Usage: python pools_pf07938_scores.py POOLS_DIR SELECTION_TSV WHOLE_PROTEOME_HITS_TSV
Prints each pool S protein with a PF07938 annotation, its pool, and its Pfam 38.2 GA score
from the review 2 whole-proteome table (blank if not a GA hit there).
"""

import csv
import gzip
import os
import sys

pools_dir, sel_tsv, hits_tsv = sys.argv[1:4]
sel = list(csv.DictReader(open(sel_tsv), delimiter="\t"))
grp = {r["name"]: r["group"] for r in sel}
score = {}
for r in csv.DictReader(open(hits_tsv), delimiter="\t"):
    if r["model"] == "Fungal_lectin":
        score[(r["proteome"], r["protein"])] = float(r["seq_score"])
rows = []
for name in grp:
    f = os.path.join(pools_dir, name + ".pools.tsv.gz")
    with gzip.open(f, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if not rows:
            print("columns", rd.fieldnames)
        for r in rd:
            if "PF07938" in (r.get("pfam") or "").split(";"):
                pid = r.get("protein") or r.get("protein_id") or list(r.values())[0]
                rows.append((name, grp[name], pid, r.get("pool"), score.get((name, pid))))
for x in sorted(rows, key=lambda x: -(x[4] or 0)):
    print(*x, sep="\t")
sc = [x[4] for x in rows if x[4] is not None]
print("pool S PF07938 proteins", len(rows), "with GA hit in whole-proteome table", len(sc))
for t in (50, 100, 200, 300):
    print(
        f">={t}",
        sum(1 for s in sc if s >= t),
        "proteins",
        len({x[0] for x in rows if x[4] is not None and x[4] >= t}),
        "genera",
    )
