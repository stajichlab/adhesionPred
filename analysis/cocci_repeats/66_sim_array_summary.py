#!/usr/bin/python3.12
"""Summarise the read simulation (65): what array/flank slope does the pipeline give when the
true unit count is known?

array/flank is computed exactly as in 40_sowgp_array_depth.py: length-weighted mean depth over
the 'array' BED lines, divided by the length-weighted mean over 'nflank' + 'cflank'.

Inputs : sim_array_depth.regions.tsv.gz, sim_alleles.tsv (both from 65_sim_array_depth.sh)
Outputs: sim_array_depth.tsv (one row per simulated run), sim_array_depth.png, stdout
Run    : /usr/bin/python3.12 66_sim_array_summary.py | tee sim_array_summary.log
"""

import matplotlib
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = pd.read_csv("sim_array_depth.regions.tsv.gz", sep="\t")
A = pd.read_csv("sim_alleles.tsv", sep="\t")
R["w"] = R.end - R.start
R["grp"] = R.region.where(R.region == "array", "flank")
g = (
    R.assign(dw=R.depth * R.w)
    .groupby(["run", "allele", "read_len", "cov", "rep", "mapq", "grp"])[["dw", "w"]]
    .sum()
    .reset_index()
)
g["mean"] = g.dw / g.w
W = g.pivot_table(
    index=["run", "allele", "read_len", "cov", "rep"], columns=["mapq", "grp"], values="mean"
)
W.columns = [f"{m}_{r}" for m, r in W.columns]
W = W.reset_index().merge(A, on="allele")
W["array_flank"] = W.all_array / W.all_flank
W["array_flank_q20"] = W.q20_array / W.q20_flank
W["q20_frac_array"] = W.q20_array / W.all_array
W["predicted"] = W.true_units / 4
W.to_csv("sim_array_depth.tsv", sep="\t", index=False, float_format="%.4f")
print(f"runs: {len(W)}; flank depth range {W.all_flank.min():.1f}-{W.all_flank.max():.1f}")


def slope(df, y="array_flank"):
    r = stats.linregress(df.true_units, df[y])
    tq = stats.t.ppf(0.975, len(df) - 2)
    return r.slope, r.slope - tq * r.stderr, r.slope + tq * r.stderr, r.intercept


rs = W[W.source == "RS edit"]
print("\n== RS edits: median array/flank by true unit count (all read lengths and coverages)")
print(
    rs.groupby("true_units")[["array_flank", "array_flank_q20", "q20_frac_array"]]
    .median()
    .round(3)
    .to_string()
)
for y in ("array_flank", "array_flank_q20"):
    s, lo, hi, b0 = slope(rs, y)
    print(
        f"\nRS edits, all runs (n={len(rs)}): {y} slope per unit {s:.3f} (95% CI {lo:.3f} to {hi:.3f}); "
        f"intercept {b0:.3f}; predicted slope 0.250"
    )
    r5 = rs[rs.true_units <= 5]
    s, lo, hi, b0 = slope(r5, y)
    print(
        f"RS edits, u = 2-5 only (the range in the real data; n={len(r5)}): {y} slope {s:.3f} "
        f"(95% CI {lo:.3f} to {hi:.3f}); intercept {b0:.3f}"
    )
    rho, p = stats.spearmanr(r5.true_units, r5[y])
    print(f"   Spearman rho {rho:+.3f} (p={p:.2g})")
print("\n== RS edits: slope by read length and coverage (array_flank; array_flank_q20)")
rows = []
for (L, c), df in rs.groupby(["read_len", "cov"]):
    s, lo, hi, _ = slope(df)
    sq, loq, hiq, _ = slope(df, "array_flank_q20")
    rows.append(
        (
            L,
            c,
            len(df),
            round(s, 3),
            round(lo, 3),
            round(hi, 3),
            round(sq, 3),
            round(loq, 3),
            round(hiq, 3),
        )
    )
print(
    pd.DataFrame(
        rows, columns=["read_len", "cov", "n", "slope", "lo", "hi", "slope_q20", "lo_q20", "hi_q20"]
    ).to_string(index=False)
)
print("\n== RS edits: median array/flank by units and read length")
print(
    rs.pivot_table(index="true_units", columns="read_len", values="array_flank", aggfunc="median")
    .round(3)
    .to_string()
)
print("\n== real long-read alleles: median array/flank (all reads) and Q20, against u/4")
re_ = W[W.source != "RS edit"]
print(
    re_.groupby(["allele", "true_units"])
    .agg(
        af=("array_flank", "median"),
        af_q20=("array_flank_q20", "median"),
        q20_frac=("q20_frac_array", "median"),
        n=("run", "size"),
    )
    .assign(predicted=lambda x: x.index.get_level_values(1) / 4)
    .round(3)
    .to_string()
)
print("\n== real alleles by read length (median array/flank)")
print(
    re_.pivot_table(index="allele", columns="read_len", values="array_flank", aggfunc="median")
    .round(3)
    .to_string()
)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for axi, y, t in zip(
    ax, ("array_flank", "array_flank_q20"), ("all reads", "MAPQ >= 20"), strict=True
):
    for L, mk in zip((100, 150, 250, 300), ("o", "s", "^", "D"), strict=True):
        df = rs[rs.read_len == L]
        axi.scatter(
            df.true_units + (L - 200) / 1000, df[y], s=10, marker=mk, alpha=0.6, label=f"{L} nt"
        )
    axi.plot([2, 6], [0.5, 1.5], "r--", lw=1, label="u/4")
    axi.set_xlabel("true SOWgp units (simulated RS edits)")
    axi.set_title(f"simulated reads mapped to RS, {t}")
ax[0].set_ylabel("array depth / flank depth")
ax[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig("plots/sim_array_depth.png", dpi=150)
print("\nwrote sim_array_depth.tsv, plots/sim_array_depth.png")
