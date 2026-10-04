#!/usr/bin/env python3.12
"""Read depth over the SOWgp repeat array vs repeat unit count.

Inputs
  mosdepth per-strain output from pipeline/11c_mosdepth_regions.sh (Genotyping project):
    <regions-dir>/<strain>.CIMG_04613_array.{all,Q20}.regions.bed.gz
    BED lines named nflank / array / cflank (CDS segments of CIMG_04613, RS reference)
  sowgp_tree_units.tsv      unit count per strain (script 13)
  CIMG_04613.coverage.tsv   whole-gene depth and genome-wide depth

Prediction (written before the data were seen): RS has 4 units. If reads from a strain with u
units map inside the array, array depth / flank depth is about u/4.

Outputs: sowgp_array_depth.tsv, sowgp_array_depth.png, and statistics on stdout.
"""

import argparse
import csv
import glob
import gzip
import os
import re
from collections import defaultdict

import matplotlib
import numpy as np
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RS_UNITS = 4
D = "/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS"
ap = argparse.ArgumentParser()
ap.add_argument("--regions-dir", default=f"{D}/coverage/mosdepth_regions")
ap.add_argument("--tag", default="CIMG_04613_array")
ap.add_argument("--units", default="sowgp_tree_units.tsv")
ap.add_argument("--gene", default="CIMG_04613.coverage.tsv")
ap.add_argument("--out", default="sowgp_array_depth")
a = ap.parse_args()


def weighted(path):
    """Length-weighted mean depth per region name."""
    num, den = defaultdict(float), defaultdict(float)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, name, m = line.rstrip("\n").split("\t")
            num[name] += float(m) * (int(e) - int(s))
            den[name] += int(e) - int(s)
    return {k: num[k] / den[k] for k in num}


strains = sorted(
    os.path.basename(p)[: -len(f".{a.tag}.all.regions.bed.gz")]
    for p in glob.glob(f"{a.regions_dir}/*.{a.tag}.all.regions.bed.gz")
)
gene = {r["strain"]: r for r in csv.DictReader(open(a.gene), delimiter="\t")}
tips = {}
for r in csv.DictReader(open(a.units), delimiter="\t"):
    tips[re.sub(r"^Coccidioides_(immitis|posadasii)_", "", r["tip"])] = r

rows = []
for s in strains:
    al = weighted(f"{a.regions_dir}/{s}.{a.tag}.all.regions.bed.gz")
    q = weighted(f"{a.regions_dir}/{s}.{a.tag}.Q20.regions.bed.gz")
    fl_all = (al["nflank"] * 249 + al["cflank"] * 165) / 414
    fl_q = (q["nflank"] * 249 + q["cflank"] * 165) / 414
    t = tips.get(s)
    rows.append(
        {
            "strain": s,
            "status": t["status"] if t else "not in tree",
            "clade": t["clade_species"] if t else "",
            "n_units": t["n_units"] if t and t["n_units"] else "",
            "genome_mean": float(gene[s]["genome_mean"]),
            "array": al["array"],
            "flank": fl_all,
            "array_q20": q["array"],
            "flank_q20": fl_q,
            "array_flank": al["array"] / fl_all if fl_all > 0 else float("nan"),
            "array_flank_q20": q["array"] / fl_q if fl_q > 0 else float("nan"),
            "q20_frac_array": q["array"] / al["array"] if al["array"] > 0 else float("nan"),
            "q20_frac_flank": fl_q / fl_all if fl_all > 0 else float("nan"),
        }
    )

cols = list(rows[0])
with open(a.out + ".tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write(
            "\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n"
        )
print(f"{len(rows)} strains with array depth")

# ---- statistics, full-length strains with >= 10x genome depth and flank depth > 0
fl = [
    r
    for r in rows
    if r["status"] == "full-length"
    and r["n_units"] != ""
    and r["genome_mean"] >= 10
    and r["flank"] > 0
]
for r in fl:
    r["u"] = int(float(r["n_units"]))
print(f"full-length, genome >= 10x: n={len(fl)}")


def spear(label, sub, y):
    if len(sub) < 8:
        print(f"{label:34s} n={len(sub)}")
        return
    rho, p = stats.spearmanr([r["u"] for r in sub], [r[y] for r in sub])
    print(f"{label:34s} n={len(sub):3d} {y:16s} rho={rho:+.3f} p={p:.3g}")


for y in ("array_flank", "array_flank_q20", "q20_frac_array", "q20_frac_flank"):
    spear("both species", fl, y)
    for sp in ("immitis", "posadasii"):
        spear(sp, [r for r in fl if r["clade"] == sp], y)

print("\nmedian array/flank (unfiltered) by unit count; prediction u/4")
for sp in ("immitis", "posadasii"):
    for u in range(2, 7):
        x = [r for r in fl if r["clade"] == sp and r["u"] == u]
        if x:
            v = np.array([r["array_flank"] for r in x])
            print(
                f"{sp:10s} u={u} n={len(x):3d} median={np.median(v):.2f} "
                f"IQR={np.percentile(v, 25):.2f}-{np.percentile(v, 75):.2f} predicted={u / RS_UNITS:.2f}"
            )

# linear slope of array/flank on u, both species, with 95% CI
for sp in ("immitis", "posadasii"):
    x = [r for r in fl if r["clade"] == sp]
    if len(x) > 8:
        res = stats.linregress([r["u"] for r in x], [r["array_flank"] for r in x])
        ci = 1.96 * res.stderr
        print(
            f"{sp}: slope per unit = {res.slope:.3f} (95% CI {res.slope - ci:.3f} to {res.slope + ci:.3f}); "
            f"predicted slope = {1 / RS_UNITS:.3f}"
        )

# ---- figure
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for axi, sp in zip(ax, ("immitis", "posadasii"), strict=True):
    us = sorted({r["u"] for r in fl if r["clade"] == sp})
    data = [[r["array_flank"] for r in fl if r["clade"] == sp and r["u"] == u] for u in us]
    axi.boxplot(data, positions=us, widths=0.5, showfliers=True)
    axi.plot([2, 6], [0.5, 1.5], "r--", lw=1, label="prediction u/4")
    axi.set_title(f"{sp} (n={sum(len(d) for d in data)})")
    axi.set_xlabel("SOWgp repeat units (full-length strains)")
    axi.set_xlim(1.5, 6.5)
ax[0].set_ylabel("array depth / flank depth (CDS, RS ref)")
ax[0].legend()
fig.tight_layout()
os.makedirs(os.path.join(os.path.dirname(a.out) or ".", "plots"), exist_ok=True)
fig.savefig(os.path.join(os.path.dirname(a.out) or ".", "plots", os.path.basename(a.out) + ".png"), dpi=150)
