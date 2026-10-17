"""Review 3: recompute whole-proteome genus counts from committed hit tables.

Usage: python recount_u.py SELECTION_TSV HITS_TSV [HITS_TSV ...]
HITS_TSV columns: proteome, protein, model, model_acc, seq_evalue, seq_score.
"""

import collections
import csv
import sys

sel = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
n = len(sel)
grp = {r["name"]: r["group"] for r in sel}
print("proteomes", n, "unique genera", len({r["genus"] for r in sel}))
print("groups", dict(collections.Counter(grp.values())))
print("total proteins", sum(int(r["proteins"]) for r in sel))
for path in sys.argv[2:]:
    h = list(csv.DictReader(open(path), delimiter="\t"))
    for m in sorted({x["model"] for x in h}):
        rows = [x for x in h if x["model"] == m]
        for t in [0, 50, 100, 200, 300]:
            rr = [x for x in rows if float(x["seq_score"]) >= t]
            g = {x["proteome"] for x in rr}
            prot = len({(x["proteome"], x["protein"]) for x in rr})
            notsel = sum(1 for p in g if p not in grp)
            by = dict(collections.Counter(grp.get(p, "?") for p in g))
            print(
                m,
                f">={t}",
                "proteins",
                prot,
                "genera",
                len(g),
                f"{100 * len(g) / n:.1f}%",
                by,
                "not_in_selection",
                notsel,
            )
