#!/usr/bin/python3.12
"""Depth-implied unit count for the full-length strains, calibrated by the read simulation.

Calibration: OLS of array/flank on true units for the simulated RS edits with u = 2-5
(sim_array_depth.tsv from 66), all read lengths and coverages. Then for each real strain:
    units_depth = (array_flank - intercept) / slope
This is indicative only: the simulation has uniform coverage and no GC or PCR bias, and the
calibration uses RS (immitis) sequence; for posadasii the real Silveira 2022 allele (u = 4)
gave array/flank 0.955 against 1.015 for RS u = 4 (see sim_array_summary.log).

Also: the array/flank slope on gene-model units, restricted to strains whose array is
contiguous and N-free in the annotated assembly (asm_locus_summary.tsv, 63).

Inputs : sim_array_depth.tsv, sowgp_array_depth.tsv, asm_locus_summary.tsv
Output : stdout
Run    : /usr/bin/python3.12 67_unit_depth_calibration.py | tee unit_depth_calibration.log
"""

import numpy as np
import pandas as pd
from scipy import stats

S = pd.read_csv("sim_array_depth.tsv", sep="\t")
cal = S[(S.source == "RS edit") & (S.true_units <= 5)]
c = stats.linregress(cal.true_units, cal.array_flank)
print(
    f"calibration (sim RS edits, u 2-5, n={len(cal)}): array/flank = {c.intercept:.3f} + {c.slope:.3f} x units"
)
resid_sd = np.std(cal.array_flank - (c.intercept + c.slope * cal.true_units), ddof=2)
print(f"residual SD in simulation {resid_sd:.3f} (= {resid_sd / c.slope:.2f} units)")

a = pd.read_csv("sowgp_array_depth.tsv", sep="\t")
a = a[
    (a.status == "full-length") & a.n_units.notna() & (a.genome_mean >= 10) & (a.flank > 0)
].copy()
a["u"] = a.n_units.astype(int)
a["units_depth"] = (a.array_flank - c.intercept) / c.slope
a["diff"] = a.units_depth - a.u
print("\nmedian depth-implied units by species and gene-model units")
print(
    a.groupby(["clade", "u"])
    .agg(
        n=("strain", "size"),
        median_units_depth=("units_depth", "median"),
        q25=("units_depth", lambda x: x.quantile(0.25)),
        q75=("units_depth", lambda x: x.quantile(0.75)),
    )
    .round(2)
    .to_string()
)
for sp, g in a.groupby("clade"):
    print(
        f"{sp}: depth-implied units exceed gene-model units by >= 1 in {(g['diff'] >= 1).sum()} of {len(g)} "
        f"strains ({(g['diff'] >= 1).mean():.0%}); fall short by >= 1 in {(g['diff'] <= -1).sum()}"
    )

rl = pd.read_csv("cram_readlen.tsv", sep="\t")[["strain", "modal_len"]]
a = a.merge(rl, on="strain", how="left")
a["rl_class"] = pd.cut(
    a.modal_len, [0, 101, 126, 151, 251, 301], labels=["<=101", "126", "150", "251", "301"]
)
print(
    "\nreal data: median array/flank by read-length class (full-length strains); the simulation shows"
)
print(
    "little read-length dependence (sim_array_summary.log), so a trend here is not a mapping effect"
)
print(
    a.groupby(["clade", "rl_class"], observed=True)
    .agg(n=("strain", "size"), af=("array_flank", "median"), u=("u", "mean"))
    .round(3)
    .to_string()
)
for sp, g in a.groupby("clade"):
    h = stats.kruskal(
        *[x.array_flank for _, x in g.groupby("rl_class", observed=True) if len(x) > 2]
    )
    print(f"{sp}: Kruskal-Wallis array/flank across read-length classes p={h.pvalue:.3g}")

L = pd.read_csv("asm_locus_summary.tsv", sep="\t")[["strain", "ann_class", "ann_units_asm"]]
b = a.merge(L, on="strain", how="left")
print(
    "\narray/flank slope on gene-model units, by assembly class of the array (annotated scaffolds)"
)
for sp in ("immitis", "posadasii"):
    for lab, sub in (
        ("all", b[b.clade == sp]),
        ("contiguous, N-free", b[(b.clade == sp) & (b.ann_class == "contiguous")]),
        ("N gap or break", b[(b.clade == sp) & (b.ann_class != "contiguous")]),
    ):
        if len(sub) > 5 and sub.u.nunique() > 1:
            r = stats.linregress(sub.u, sub.array_flank)
            tq = stats.t.ppf(0.975, len(sub) - 2)
            print(
                f"{sp:10s} {lab:20s} n={len(sub):3d} slope {r.slope:+.3f} (95% CI {r.slope - tq * r.stderr:+.3f} "
                f"to {r.slope + tq * r.stderr:+.3f}); median array/flank {sub.array_flank.median():.3f}"
            )
        else:
            print(f"{sp:10s} {lab:20s} n={len(sub):3d} (too few for a slope)")
