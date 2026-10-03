#!/usr/bin/env python3.12
"""Whole-gene read depth at CIMG_04613 vs SOWgp unit count (report 2026-10-01, section 4.2).

Inputs : sowgp_tree_units.tsv, CIMG_04613.coverage.tsv
Output : stdout (Spearman rho for the depth ratio, and for the fraction of reads kept at Q20)
Strains: full-length SOWgp copy and genome-wide depth >= 10x.
"""

import csv
import re

import numpy as np
from scipy import stats


def key(t):
    return re.sub(r"^Coccidioides_(immitis|posadasii)_", "", t)


d = {r["strain"]: r for r in csv.DictReader(open("CIMG_04613.coverage.tsv"), delimiter="\t")}
rows = []
for t in csv.DictReader(open("sowgp_tree_units.tsv"), delimiter="\t"):
    k = key(t["tip"])
    if k in d and t["status"] == "full-length" and t["n_units"]:
        x = d[k]
        m = float(x["mean_depth"])
        if float(x["genome_mean"]) >= 10 and m > 0:
            rows.append(
                {
                    "u": int(float(t["n_units"])),
                    "sp": t["clade_species"],
                    "ratio": float(x["ratio"]),
                    "q20": float(x["mean_depth_Q20"]) / m,
                }
            )
print(f"n={len(rows)}")
for label, sub in (
    ("both species", rows),
    ("immitis", [r for r in rows if r["sp"] == "immitis"]),
    ("posadasii", [r for r in rows if r["sp"] == "posadasii"]),
):
    for y in ("ratio", "q20"):
        rho, p = stats.spearmanr([r["u"] for r in sub], [r[y] for r in sub])
        print(f"{label:14s} n={len(sub):3d} {y:6s} rho={rho:+.3f} p={p:.3g}")
print("median ratio by unit count:")
for sp in ("immitis", "posadasii"):
    for u in range(2, 7):
        x = [r["ratio"] for r in rows if r["sp"] == sp and r["u"] == u]
        if x:
            print(f"  {sp:10s} u={u} n={len(x):3d} median={np.median(x):.2f}")
