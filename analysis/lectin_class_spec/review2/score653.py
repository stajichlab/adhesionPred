import glob
import os
import sys

scan_dir = sys.argv[1]
model = sys.argv[2]
best = {}  # (genome, protein) -> score
for f in glob.glob(os.path.join(scan_dir, "*.tbl")):
    g = os.path.basename(f)[:-4]
    for line in open(f):
        if line.startswith("#"):
            continue
        p = line.split()
        if p[2] != model:
            continue
        best[(g, p[0])] = max(best.get((g, p[0]), 0), float(p[5]))
bins = [25.3, 50, 100, 200, 300]
for t in bins:
    gs = {g for (g, _), s in best.items() if s >= t}
    n = sum(1 for s in best.values() if s >= t)
    print(f"{model} seq score >= {t}: proteins {n}, genera {len(gs)} ({100*len(gs)/653:.1f}%)")
