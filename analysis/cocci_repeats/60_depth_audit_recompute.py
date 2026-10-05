#!/usr/bin/python3.12
"""Audit of REPORT_2026-10-01_sowgp_depth_vs_repeats.md: recompute every number from source.

Sections (each prints a block on stdout; tables are written as TSV):
  R0  rebuild CIMG_04613.coverage.tsv and the array/flank table from the raw mosdepth files,
      and compare with the committed tables
  R1  tip matching and status counts (report section 2 and table 4.1)
  R2  table 4.1, the 39 no-model strains, and the three Mann-Whitney tests (script 41),
      with effect sizes (difference in medians, Hodges-Lehmann shift) and bootstrap 95% CIs
  R3  the 7 "deletion candidates": sensitivity to the ratio and genome-depth thresholds
  R4  section 4.2, whole-gene ratio vs units (script 42), with bootstrap CIs for rho
  R5  section 4.3, array/flank vs units (script 40): medians, 12 Spearman tests, slopes with
      95% CI, and t-tests of slope = 0.25 and slope = 0
  R6  "Species explains this" (pooled Q20 rho = -0.17): partial Spearman given species, and a
      within-species permutation test
  R7  "base length 281 against 277 aa": observed base lengths (length - 47 x units)
  R8  every hypothesis test the scripts ran, in families; raw p, BH q and Holm p

Inputs : CIMG_04613.coverage.tsv, sowgp_array_depth.tsv, sowgp_tree_units.tsv, and the raw
         mosdepth outputs under Genotyping/all_C_immitis_ref_RS/coverage/
Outputs: audit_tests.tsv (R8), audit_deletion_sensitivity.tsv (R3), stdout
Run    : /usr/bin/python3.12 60_depth_audit_recompute.py | tee audit_recompute.log
"""

import gzip

import numpy as np
import pandas as pd
from scipy import stats

RNG = np.random.default_rng(20261004)
NBOOT = 4000
D = (
    "/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/"
    "Genotyping/all_C_immitis_ref_RS/coverage"
)
TESTS = []  # rows for R8


def hdr(s):
    print(f"\n{'=' * 8} {s} {'=' * 8}")


def add_test(family, tid, desc, script, reported, n, effect, lo, hi, p, note=""):
    TESTS.append(
        {
            "family": family,
            "test_id": tid,
            "description": desc,
            "script": script,
            "in_report": reported,
            "n": n,
            "effect": effect,
            "ci_low": lo,
            "ci_high": hi,
            "p_raw": p,
            "note": note,
        }
    )


def boot_ci(fn, *arrays, n=NBOOT):
    """Percentile bootstrap 95% CI of fn over independent resamples of each array."""
    vals = []
    for _ in range(n):
        res = [a[RNG.integers(0, len(a), len(a))] for a in arrays]
        vals.append(fn(*res))
    vals = np.array(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


def boot_ci_paired(fn, x, y, n=NBOOT):
    vals = []
    for _ in range(n):
        i = RNG.integers(0, len(x), len(x))
        vals.append(fn(x[i], y[i]))
    vals = np.array(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


def hl_shift(x, y):
    """Hodges-Lehmann estimate of the location shift x - y."""
    return np.median(np.subtract.outer(x, y))


def rho(x, y):
    return stats.spearmanr(x, y)[0]


def key(t):
    return t.str.replace(r"^Coccidioides_(immitis|posadasii)_", "", regex=True)


# ------------------------------------------------------------------ R0
hdr("R0 rebuild tables from raw mosdepth files")
cov = pd.read_csv("CIMG_04613.coverage.tsv", sep="\t")
rows = []
for s in cov.strain:

    def region(path):
        with gzip.open(path, "rt") as fh:
            return float(fh.readline().split("\t")[-1])

    g = None
    with open(f"{D}/mosdepth/{s}.5000bp.mosdepth.summary.txt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] == "total":
                g = float(f[3])
    m = region(f"{D}/mosdepth_gene/{s}.all.regions.bed.gz")
    mq = region(f"{D}/mosdepth_gene/{s}.Q20.regions.bed.gz")
    rows.append((s, m, mq, g))
raw = pd.DataFrame(rows, columns=["strain", "m", "mq", "g"]).merge(cov, on="strain")
print(
    f"strains: {len(raw)}; max |diff| mean_depth {np.max(np.abs(raw.m - raw.mean_depth)):.3f}, "
    f"Q20 {np.max(np.abs(raw.mq - raw.mean_depth_Q20)):.3f}, genome {np.max(np.abs(raw.g - raw.genome_mean)):.3f}, "
    f"ratio {np.max(np.abs(raw.m / raw.g - raw.ratio)):.4f}"
)

arr = pd.read_csv("sowgp_array_depth.tsv", sep="\t")
diffs = []
for s in arr.strain:
    out = {}
    for tag in ("all", "Q20"):
        num, den = {}, {}
        with gzip.open(
            f"{D}/mosdepth_regions/{s}.CIMG_04613_array.{tag}.regions.bed.gz", "rt"
        ) as fh:
            for line in fh:
                c, st, en, nm, v = line.rstrip("\n").split("\t")
                w = int(en) - int(st)
                num[nm] = num.get(nm, 0) + float(v) * w
                den[nm] = den.get(nm, 0) + w
        mean = {k: num[k] / den[k] for k in num}
        fl = (mean["nflank"] * den["nflank"] + mean["cflank"] * den["cflank"]) / (
            den["nflank"] + den["cflank"]
        )
        out[tag] = (mean["array"], fl, den)
    diffs.append(
        (
            s,
            out["all"][0] / out["all"][1] if out["all"][1] else np.nan,
            out["Q20"][0] / out["Q20"][1] if out["Q20"][1] else np.nan,
        )
    )
rd = pd.DataFrame(diffs, columns=["strain", "af", "afq"]).merge(arr, on="strain")
print(
    f"array table: {len(rd)} strains; max |diff| array/flank {np.nanmax(np.abs(rd.af - rd.array_flank)):.4f}, "
    f"Q20 {np.nanmax(np.abs(rd.afq - rd.array_flank_q20)):.4f}; BED lengths {out['all'][2]}"
)
print(
    "Note: script 40 hard-codes flank lengths 249 and 165; the BED lengths above are what mosdepth used."
)

# ------------------------------------------------------------------ R1
hdr("R1 tip matching and status counts")
tu = pd.read_csv("sowgp_tree_units.tsv", sep="\t")
tu["strain"] = key(tu.tip)
print(f"tree tips {len(tu)}; duplicated keys {tu.strain.duplicated().sum()}")
m = tu[tu.strain.isin(cov.strain)]
print(f"tips matched to a depth strain: {len(m)} (report: 484)")
print("unmatched tips:", list(tu.loc[~tu.strain.isin(cov.strain), "tip"]))
print(f"depth strains with no tree tip: {(~cov.strain.isin(tu.strain)).sum()} (report: 75)")
print("status among tree tips:", tu.status.value_counts().to_dict())
print("status among matched tips:", m.status.value_counts().to_dict())
st = cov.merge(tu[["strain", "status", "clade_species", "n_units"]], on="strain", how="left")
st["status"] = st.status.fillna("not in tree")
print(f"table total: {st.status.value_counts().to_dict()} sum={len(st)}")
fl = st[st.status == "full-length"]
print("full-length with depth, by species and units:")
print(fl.groupby(["clade_species", "n_units"]).size().to_string())
print(
    f"full-length with genome depth >= 10x: {(fl.genome_mean >= 10).sum()} of {len(fl)}; min genome depth {fl.genome_mean.min():.2f}"
)

# ------------------------------------------------------------------ R2
hdr("R2 table 4.1 and status vs ratio (script 41)")
groups = ["full-length", "fragment only", "no SOWgp gene model", "not in tree"]
for g in groups:
    x = st[st.status == g]
    print(
        f"{g:22s} n={len(x):3d} median ratio={x.ratio.median():.3f} ratio<0.5={int((x.ratio < 0.5).sum())} "
        f"ratio<0.5&genome>=10={int(((x.ratio < 0.5) & (x.genome_mean >= 10)).sum())} "
        f"median genome={x.genome_mean.median():.1f}"
    )
nm = st[st.status == "no SOWgp gene model"].sort_values("ratio")
print(
    f"no-model ratio range {nm.ratio.min():.3f}-{nm.ratio.max():.3f}; lowest {nm.strain.iloc[0]} {nm.ratio.iloc[0]:.3f}"
)
ok = st[st.genome_mean >= 10]
print("Mann-Whitney on genome depth >= 10x (as script 41):")
for i, (a, b) in enumerate(
    (
        ("no SOWgp gene model", "full-length"),
        ("fragment only", "full-length"),
        ("no SOWgp gene model", "fragment only"),
    ),
    1,
):
    x = ok.loc[ok.status == a, "ratio"].to_numpy()
    y = ok.loc[ok.status == b, "ratio"].to_numpy()
    p = stats.mannwhitneyu(x, y).pvalue
    dmed = np.median(x) - np.median(y)
    lo, hi = boot_ci(lambda u, v: np.median(u) - np.median(v), x, y)
    hl = hl_shift(x, y)
    hlo, hhi = boot_ci(hl_shift, x, y, n=1000)
    print(
        f"  {a} (n={len(x)}) vs {b} (n={len(y)}): medians {np.median(x):.3f}/{np.median(y):.3f}; "
        f"diff {dmed:+.3f} (95% CI {lo:+.3f} to {hi:+.3f}); HL shift {hl:+.3f} ({hlo:+.3f} to {hhi:+.3f}); p={p:.3g}"
    )
    add_test(
        "F1 status vs whole-gene ratio",
        f"F1.{i}",
        f"MW ratio {a} vs {b}",
        "41",
        "yes" if i < 3 else "no (run by 41, not reported)",
        len(x) + len(y),
        hl,
        hlo,
        hhi,
        p,
        "effect = Hodges-Lehmann shift in ratio",
    )
print(
    "genome depth by group (all strains): ",
    {g: round(st.loc[st.status == g, "genome_mean"].median(), 1) for g in groups},
)
x = ok.loc[ok.status == "no SOWgp gene model", "genome_mean"]
y = ok.loc[ok.status == "full-length", "genome_mean"]
print(f"genome depth no-model vs full-length (>=10x): MW p={stats.mannwhitneyu(x, y).pvalue:.3g}")
r, p = stats.spearmanr(ok.genome_mean, ok.ratio)
print(f"ratio vs genome depth, all strains >=10x: Spearman rho={r:+.3f} p={p:.3g} n={len(ok)}")

# ------------------------------------------------------------------ R3
hdr("R3 deletion candidates: threshold sensitivity")
srt = st.sort_values("ratio")
print("lowest 20 ratios:")
print(srt.head(20)[["strain", "status", "ratio", "genome_mean"]].to_string(index=False))
sens = []
for rt in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
    for gt in (3, 5, 10, 15, 20, 30):
        sel = st[(st.ratio < rt) & (st.genome_mean >= gt)]
        sens.append(
            {
                "ratio_lt": rt,
                "genome_ge": gt,
                "n": len(sel),
                "strains": ",".join(sorted(sel.strain)),
            }
        )
S = pd.DataFrame(sens)
S.to_csv("audit_deletion_sensitivity.tsv", sep="\t", index=False)
print(S.pivot(index="ratio_lt", columns="genome_ge", values="n").to_string())
base = set(S[(S.ratio_lt == 0.5) & (S.genome_ge == 10)].strains.iloc[0].split(","))
print("report set (0.5, 10x):", sorted(base))
gaps = np.diff(np.sort(st.ratio.to_numpy())[:30])
print("gaps between consecutive sorted ratios (first 30 strains):", np.round(gaps, 3).tolist())

# ------------------------------------------------------------------ R4
hdr("R4 section 4.2: whole-gene ratio vs unit count (script 42)")
f = fl[(fl.genome_mean >= 10) & (fl.mean_depth > 0)].copy()
f["q20"] = f.mean_depth_Q20 / f.mean_depth
for _i, (lab, sub) in enumerate(
    (
        ("both species", f),
        ("immitis", f[f.clade_species == "immitis"]),
        ("posadasii", f[f.clade_species == "posadasii"]),
    ),
    1,
):
    for y, fam, rep in (
        ("ratio", "F2 depth vs unit count", "yes"),
        ("q20", "F4 Q20 retention vs unit count", "no (run by 42, not reported)"),
    ):
        u, v = sub.n_units.to_numpy(), sub[y].to_numpy()
        r, p = stats.spearmanr(u, v)
        lo, hi = boot_ci_paired(rho, u, v)
        print(
            f"{lab:13s} n={len(sub):3d} {y:6s} rho={r:+.3f} (95% CI {lo:+.3f} to {hi:+.3f}) p={p:.3g}"
        )
        add_test(
            fam,
            f"{'F2' if y == 'ratio' else 'F4'}.gene.{lab}",
            f"Spearman whole-gene {y} vs units, {lab}",
            "42",
            rep,
            len(sub),
            r,
            lo,
            hi,
            p,
            "effect = Spearman rho",
        )

# ------------------------------------------------------------------ R5
hdr("R5 section 4.3: array/flank vs unit count (script 40)")
a = arr[
    (arr.status == "full-length") & arr.n_units.notna() & (arr.genome_mean >= 10) & (arr.flank > 0)
].copy()
a["u"] = a.n_units.astype(int)
print(f"n={len(a)}")
for sp in ("immitis", "posadasii"):
    for u in range(2, 7):
        v = a[(a.clade == sp) & (a.u == u)].array_flank
        if len(v):
            print(
                f"{sp:10s} u={u} n={len(v):3d} median={v.median():.2f} IQR={v.quantile(0.25):.2f}-{v.quantile(0.75):.2f} "
                f"(numpy linear percentile) predicted={u / 4:.2f}"
            )
reported = {
    ("array_flank", "immitis"),
    ("array_flank", "posadasii"),
    ("array_flank_q20", "immitis"),
    ("array_flank_q20", "posadasii"),
    ("q20_frac_array", "both species"),
    ("q20_frac_array", "immitis"),
    ("q20_frac_array", "posadasii"),
}
for y in ("array_flank", "array_flank_q20", "q20_frac_array", "q20_frac_flank"):
    fam = (
        "F2 depth vs unit count"
        if y.startswith("array_flank")
        else "F4 Q20 retention vs unit count"
    )
    for lab, sub in (
        ("both species", a),
        ("immitis", a[a.clade == "immitis"]),
        ("posadasii", a[a.clade == "posadasii"]),
    ):
        u, v = sub.u.to_numpy(), sub[y].to_numpy()
        r, p = stats.spearmanr(u, v)
        lo, hi = boot_ci_paired(rho, u, v)
        print(
            f"{lab:13s} n={len(sub):3d} {y:16s} rho={r:+.3f} (95% CI {lo:+.3f} to {hi:+.3f}) p={p:.3g}"
        )
        add_test(
            fam,
            f"{fam[:2]}.{y}.{lab}",
            f"Spearman {y} vs units, {lab}",
            "40",
            "yes" if (y, lab) in reported else "no (run by 40, not reported)",
            len(sub),
            r,
            lo,
            hi,
            p,
            "effect = Spearman rho",
        )
for sp in ("immitis", "posadasii"):
    sub = a[a.clade == sp]
    res = stats.linregress(sub.u, sub.array_flank)
    df = len(sub) - 2
    tq = stats.t.ppf(0.975, df)
    p25 = 2 * stats.t.sf(abs((res.slope - 0.25) / res.stderr), df)
    print(
        f"{sp}: slope {res.slope:.3f} (95% CI z: {res.slope - 1.96 * res.stderr:.3f} to {res.slope + 1.96 * res.stderr:.3f}; "
        f"t: {res.slope - tq * res.stderr:.3f} to {res.slope + tq * res.stderr:.3f}); "
        f"p(slope=0.25)={p25:.3g}; p(slope=0)={res.pvalue:.3g}; n={len(sub)}"
    )
    add_test(
        "F3 u/4 prediction",
        f"F3.slope.{sp}",
        f"OLS slope of array/flank on units = 0.25, {sp}",
        "40",
        "yes (as a CI)",
        len(sub),
        res.slope,
        res.slope - tq * res.stderr,
        res.slope + tq * res.stderr,
        p25,
        "effect = slope per unit; H0 slope = 0.25",
    )
    # robustness: drop u=2 and u=5 tails; Theil-Sen slope
    ts = stats.theilslopes(sub.array_flank, sub.u)
    print(f"   Theil-Sen slope {ts.slope:.3f} (95% CI {ts.low_slope:.3f} to {ts.high_slope:.3f})")

# ------------------------------------------------------------------ R6
hdr("R6 'Species explains this' (pooled q20_frac_array vs units)")
for y in ("q20_frac_array", "q20_frac_flank"):
    xi = a.loc[a.clade == "immitis", y]
    xp = a.loc[a.clade == "posadasii", y]
    print(
        f"{y}: median immitis {xi.median():.3f} (n={len(xi)}), posadasii {xp.median():.3f} (n={len(xp)}); "
        f"MW p={stats.mannwhitneyu(xi, xp).pvalue:.3g}"
    )
print("units by species (mean):", a.groupby("clade").u.mean().round(3).to_dict())
# partial Spearman: ranks residualised on species
ru = stats.rankdata(a.u)
ry = stats.rankdata(a.q20_frac_array)
sp = (a.clade == "posadasii").astype(float).to_numpy()
X = np.column_stack([np.ones(len(a)), sp])


def resid(v):
    return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]


pr, _ = stats.pearsonr(resid(ru), resid(ry))
dfp = len(a) - 3
tstat = pr * np.sqrt(dfp / (1 - pr**2))
pp = 2 * stats.t.sf(abs(tstat), dfp)
pb = []
for _ in range(NBOOT):
    i = RNG.integers(0, len(a), len(a))
    Xi = X[i]
    rr = [
        v - Xi @ np.linalg.lstsq(Xi, v, rcond=None)[0]
        for v in (stats.rankdata(a.u.to_numpy()[i]), stats.rankdata(a.q20_frac_array.to_numpy()[i]))
    ]
    pb.append(stats.pearsonr(*rr)[0])
plo, phi = np.percentile(pb, [2.5, 97.5])
print(
    f"partial Spearman (q20_frac_array, units | species) = {pr:+.3f} (95% CI {plo:+.3f} to {phi:+.3f}), "
    f"p={pp:.3g}, n={len(a)}"
)
# within-species permutation: permute units within species, statistic = pooled rho
obs = rho(a.u.to_numpy(), a.q20_frac_array.to_numpy())
null = []
uu = a.u.to_numpy().copy()
idx = [np.where(a.clade.to_numpy() == s)[0] for s in ("immitis", "posadasii")]
for _ in range(10000):
    perm = uu.copy()
    for ix in idx:
        perm[ix] = RNG.permutation(uu[ix])
    null.append(rho(perm, a.q20_frac_array.to_numpy()))
null = np.array(null)
pperm = (np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (len(null) + 1)
print(
    f"pooled rho {obs:+.3f}; within-species permutation null mean {null.mean():+.3f}, "
    f"95% range {np.percentile(null, 2.5):+.3f} to {np.percentile(null, 97.5):+.3f}; two-sided p={pperm:.3g}"
)
add_test(
    "F4 Q20 retention vs unit count",
    "F4.partial.q20_frac_array",
    "partial Spearman q20_frac_array vs units | species",
    "60 (new)",
    "no (new test of the 'species explains this' claim)",
    len(a),
    pr,
    plo,
    phi,
    pp,
    "effect = partial rho",
)

# ------------------------------------------------------------------ R7
hdr("R7 base length (length - 47 x units), full-length tips")
t = tu[tu.status == "full-length"].copy()
t["base"] = t.length - 47 * t.n_units
for sp_, g in t.groupby("clade_species"):
    vc = g.base.value_counts().sort_index()
    print(
        f"{sp_}: n={len(g)} modal base {g.base.mode()[0]:.0f}; range {g.base.min():.0f}-{g.base.max():.0f}; {vc.to_dict()}"
    )
    print(
        f"   modal length by units: {g.groupby('n_units').length.agg(lambda s: s.mode()[
                0
            ]).to_dict()}"
    )

# ------------------------------------------------------------------ R9
hdr("R9 section 3: gene structure from the RS GFF and genome")
REF = D.replace("/coverage", "/genome/FungiDB-68_CimmitisRS_Genome.fasta")
seq, keep = [], False
with open(REF) as fh:
    for line in fh:
        if line.startswith(">"):
            if keep:
                break
            keep = line[1:].split()[0] == "GG704914"
        elif keep:
            seq.append(line.strip())
g = "".join(seq).upper()
cds_seg = [
    (970129, 970235),
    (970316, 970423),
    (970477, 970599),
    (970656, 970796),
    (970853, 970993),
    (971050, 971190),
    (971247, 971460),
]
cds = "".join(g[s - 1 : e] for s, e in cds_seg)
codon = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
}
prot = "".join(codon[cds[i : i + 3]] for i in range(0, len(cds) - 2, 3))
anch = [i + 1 for i in range(len(prot)) if prot.startswith("PTDCYGDC", i)]
print(f"CDS {len(cds)} nt, protein {len(prot.rstrip('*'))} aa, anchors PTDCYGDC at aa {anch}")
print(
    f"exon lengths (CDS parts): {[e - s + 1 for s, e in cds_seg]}; "
    f"introns between CDS parts: {[cds_seg[i + 1][0] - cds_seg[i][1] - 1 for i in range(6)]}"
)
# CDS position (nt, 1-based) of each exon end -> aa position of the intron
cum = np.cumsum([e - s + 1 for s, e in cds_seg])
for k in (2, 3, 4, 5):
    aa = cum[k] / 3
    unit_of = [i for i, s in enumerate(anch) if s <= aa + 1]
    pos_in_unit = aa - anch[unit_of[-1]] + 1 if unit_of else None
    print(
        f"intron after CDS part {k + 1}: after aa {aa:.2f}; within anchored unit {len(unit_of)} at unit aa {pos_in_unit}"
    )
arr_bed = [(970510, 970599), (970655, 970796), (970852, 970993), (971049, 971190), (971246, 971295)]
print(f"array BED total {sum(e - s for s, e in arr_bed)} nt; flank 249 + 165 = 414 nt")

# ------------------------------------------------------------------ R8
hdr("R8 multiple testing")
T = pd.DataFrame(TESTS)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def holm(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum(1, np.maximum.accumulate(p[o] * (n - np.arange(n))))
    return adj


T["q_bh_family"] = np.nan
T["p_holm_family"] = np.nan
for _, g in T.groupby("family"):
    T.loc[g.index, "q_bh_family"] = bh(g.p_raw)
    T.loc[g.index, "p_holm_family"] = holm(g.p_raw)
orig = T[~T.script.str.startswith("60")]
T["q_bh_all_original"] = np.nan
T.loc[orig.index, "q_bh_all_original"] = bh(orig.p_raw)
primary = [
    "F1.1",
    "F1.2",
    "F1.3",
    "F3.slope.immitis",
    "F3.slope.posadasii",
    "F2.array_flank.immitis",
    "F2.array_flank.posadasii",
]
T["primary"] = T.test_id.isin(primary)
pr_ = T[T.primary]
T.loc[pr_.index, "q_bh_primary"] = bh(pr_.p_raw)
T.loc[pr_.index, "p_holm_primary"] = holm(pr_.p_raw)
T.to_csv("audit_tests.tsv", sep="\t", index=False, float_format="%.4g")
print(
    f"tests run by scripts 40-42: {len(orig)}; shown in the report: {(orig.in_report.str.startswith('yes')).sum()}"
)
pd.set_option("display.width", 250)
print(
    T[
        [
            "test_id",
            "in_report",
            "n",
            "effect",
            "ci_low",
            "ci_high",
            "p_raw",
            "q_bh_family",
            "p_holm_family",
            "q_bh_all_original",
            "q_bh_primary",
            "p_holm_primary",
        ]
    ].to_string(index=False, float_format=lambda v: f"{v:.3g}")
)
