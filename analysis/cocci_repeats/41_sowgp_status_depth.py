#!/usr/bin/env python3.12
"""Join SOWgp pangenome status (script 13 table) to read depth at CIMG_04613.

Tests whether strains with no SOWgp gene model, or only a fragment, have low read depth at the
RS locus (as a deletion would) or normal depth (as an assembly/annotation gap would).

Inputs : sowgp_tree_units.tsv, CIMG_04613.coverage.tsv (whole-gene depth, RS reference)
Outputs: sowgp_status_depth.tsv (all joined strains), sowgp_missing_model_depth.tsv (no-model strains)
"""

import csv
import re

import numpy as np
from scipy import stats

d = {r["strain"]: r for r in csv.DictReader(open("CIMG_04613.coverage.tsv"), delimiter="\t")}
tips = list(csv.DictReader(open("sowgp_tree_units.tsv"), delimiter="\t"))


def key(t):
    return re.sub(r"^Coccidioides_(immitis|posadasii)_", "", t)


status = {key(t["tip"]): t["status"] for t in tips if key(t["tip"]) in d}
unmatched_tips = [t["tip"] for t in tips if key(t["tip"]) not in d]
for k in d:
    status.setdefault(k, "not in tree")
print(
    f"tree tips {len(tips)}; depth strains {len(d)}; tips without depth {len(unmatched_tips)}: {unmatched_tips}"
)


def num(k, c):
    return float(d[k][c])


with open("sowgp_status_depth.tsv", "w") as f:
    f.write("strain\tsowgp_status\tmean_depth\tgenome_mean\tratio\n")
    for k in sorted(d):
        f.write(
            f"{k}\t{status[k]}\t{num(k, 'mean_depth')}\t{num(k, 'genome_mean')}\t{num(k, 'ratio')}\n"
        )
with open("sowgp_missing_model_depth.tsv", "w") as f:
    f.write("strain\tmean_depth\tgenome_mean\tratio\n")
    for k in sorted(
        (k for k in d if status[k] == "no SOWgp gene model"), key=lambda k: num(k, "ratio")
    ):
        f.write(f"{k}\t{num(k, 'mean_depth')}\t{num(k, 'genome_mean')}\t{num(k, 'ratio')}\n")

groups = ["full-length", "fragment only", "no SOWgp gene model", "not in tree"]
print("group, n, median ratio, ratio<0.5, ratio<0.5 and genome>=10x, median genome depth")
for s in groups:
    ks = [k for k in d if status[k] == s]
    r = np.array([num(k, "ratio") for k in ks])
    g = np.array([num(k, "genome_mean") for k in ks])
    print(
        f"{s:22s} {len(ks):4d} {np.median(r):.2f} {(r < 0.5).sum():3d} {((r < 0.5) & (g >= 10)).sum():3d} {np.median(g):.1f}"
    )
ok = [k for k in d if num(k, "genome_mean") >= 10]
for a, b in (
    ("no SOWgp gene model", "full-length"),
    ("fragment only", "full-length"),
    ("no SOWgp gene model", "fragment only"),
):
    x = [num(k, "ratio") for k in ok if status[k] == a]
    y = [num(k, "ratio") for k in ok if status[k] == b]
    print(
        f"{a} vs {b} (genome>=10x): medians {np.median(x):.2f} / {np.median(y):.2f}, "
        f"Mann-Whitney p={stats.mannwhitneyu(x, y).pvalue:.2g}"
    )
print(
    "ratio<0.5 and genome>=10x:",
    [(k, status[k], round(num(k, "ratio"), 2)) for k in ok if num(k, "ratio") < 0.5],
)
