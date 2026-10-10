import collections
import csv
import os
import sys

sel_tsv, scan_dir = sys.argv[1:3]
sel = list(csv.DictReader(open(sel_tsv), delimiter="\t"))
genera = collections.defaultdict(lambda: collections.defaultdict(set))
prots = collections.Counter()
missing = 0
for r in sel:
    b = os.path.basename(r["fasta"])
    if b.endswith(".proteins.fa"):
        b = b[: -len(".proteins.fa")]
    f = os.path.join(scan_dir, b + ".tbl")
    if not os.path.exists(f):
        missing += 1
        continue
    seen = set()
    for line in open(f):
        if line.startswith("#"):
            continue
        p = line.split()
        key = (p[0], p[2])
        if key in seen:
            continue
        seen.add(key)
        genera[p[2]][r["group"]].add(r["name"])
        prots[p[2]] += 1
print("missing tbl", missing, "of", len(sel))
for m in sorted(genera):
    g = genera[m]
    tot = sum(len(v) for v in g.values())
    print(
        m,
        "proteins",
        prots[m],
        "genera",
        tot,
        f"({100*tot/len(sel):.1f}%)",
        {k: len(v) for k, v in sorted(g.items())},
    )
