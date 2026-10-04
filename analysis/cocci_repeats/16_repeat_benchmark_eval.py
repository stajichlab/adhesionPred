#!/usr/bin/python3.12
"""Score the old and new repeat detectors on the benchmark from 15_repeat_benchmark.py.

Three arms, so the improvement can be attributed:

  old            02_repeat_profile.py                         exact match, raw rate,
                                                              fractional copy count
  new-exact      14_repeat_detect_general.py --mode exact     exact match, background
                                                              correction, PSSM copy count
  new            14_repeat_detect_general.py (similarity)     BLOSUM62 match, background
                                                              correction, PSSM copy count

"new-exact" isolates what the scoring change buys from what the copy-counting change buys.

A sequence counts as DETECTED under the thresholds 03_repeat_surface_candidates.py uses:
period > 0, coverage >= 0.25, copies >= 2.5. The same rule is applied to all three arms.

Reported:
  - recall against measured mean unit-to-unit identity (the divergence floor)
  - period accuracy: exact (+-1 aa), a multiple of the truth, a divisor, or wrong
  - copy-number bias: reported minus true
  - false-positive rate on the negatives, split by kind
  - the real SOWgp alleles, where truth is the anchored unit count

Outputs: `repeat_benchmark_calls.tsv` (one row per sequence per arm),
`repeat_benchmark_metrics.tsv` (the numbers in the report), `repeat_benchmark.png`.

Usage: ./16_repeat_benchmark_eval.py
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).parent


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, str(HERE / fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


OLD = _load("rp_old", "02_repeat_profile.py")
NEW = _load("rp_new", "14_repeat_detect_general.py")

ARMS = ["old", "new-exact", "new"]
INK, MUTED, GRID, SURFACE = "#222220", "#6b6a64", "#e6e5df", "#fcfcfb"
ARM_COLOR = {"old": "#2a78d6", "new-exact": "#eda100", "new": "#eb6834"}
MIN_COVERAGE, MIN_COPIES = 0.25, 2.5


def call(arm, seq):
    if arm == "old":
        r = OLD.best_repeat(seq)
        r = {**r, "n_units": 0, "unit_period": 0, "region_entropy": 0.0}
    else:
        r = NEW.detect(seq, mode="exact" if arm == "new-exact" else "similarity")
    r["detected"] = int(
        r["period"] > 0 and r["coverage"] >= MIN_COVERAGE and float(r["n_copies"]) >= MIN_COPIES
    )
    r["detected_strict"] = int(
        r["period"] > 0 and r["coverage"] >= 0.3 and float(r["n_copies"]) >= 3.0
    )
    return r


def period_verdict(rep, true_p):
    if rep <= 0 or true_p <= 0:
        return "none"
    if abs(rep - true_p) <= 1:
        return "exact"
    for k in range(2, 9):
        if abs(rep - k * true_p) <= 1:
            return "multiple"
        if true_p % k == 0 and abs(rep - true_p // k) <= 1:
            return "divisor"
    return "wrong"


def run(fa, truth):
    seqs = dict(NEW.read_fasta(fa))
    rows = []
    for i, t in truth.iterrows():
        s = seqs[t.seq_id]
        for arm in ARMS:
            r = call(arm, s)
            rows.append(
                {
                    "seq_id": t.seq_id,
                    "arm": arm,
                    "part": t.part,
                    "donor_class": t.donor_class,
                    "length": t.length,
                    "true_period": t.true_period,
                    "true_copies": t.true_copies,
                    "target_identity": t.target_identity,
                    "true_identity": t.true_identity,
                    "rep_period": r["period"],
                    "rep_score": r["score"],
                    "rep_copies": float(r["n_copies"]),
                    "rep_coverage": r["coverage"],
                    "rep_entropy": r.get("region_entropy", 0.0),
                    "rep_unit_period": r.get("unit_period", 0),
                    "detected": r["detected"],
                    "detected_strict": r["detected_strict"],
                }
            )
        if (i + 1) % 250 == 0:
            print(f"  {i + 1}/{len(truth)}", file=sys.stderr)
    return pd.DataFrame(rows)


def ident_bin(v):
    return np.floor(v * 20) / 20 + 0.025  # 5% bins, plotted at the centre


def floor_at(recall_by_bin, level):
    """Highest identity bin centre below which recall stays under `level`."""
    s = recall_by_bin.sort_index(ascending=False)
    last = None
    for x, r in s.items():
        if r < level:
            return last, x
        last = x
    return last, None


def summarize(calls, out_tsv):
    lines = []
    syn = calls[calls.part == "synthetic"].copy()
    syn["ibin"] = ident_bin(syn.true_identity.astype(float))
    rec = syn.groupby(["arm", "ibin"]).detected.mean().unstack(0)

    for arm in ARMS:
        r = rec[arm].dropna()
        a50, b50 = floor_at(r, 0.5)
        a90, b90 = floor_at(r, 0.9)
        lines.append(
            {
                "metric": "divergence_floor",
                "arm": arm,
                "key": "recall_50",
                "value": b50 if b50 is not None else "never_below",
            }
        )
        lines.append(
            {
                "metric": "divergence_floor",
                "arm": arm,
                "key": "recall_90",
                "value": b90 if b90 is not None else "never_below",
            }
        )
        lines.append(
            {
                "metric": "recall_overall",
                "arm": arm,
                "key": "synthetic",
                "value": round(syn[syn.arm == arm].detected.mean(), 4),
            }
        )
        lines.append(
            {
                "metric": "recall_overall",
                "arm": arm,
                "key": "synthetic_strict_thresholds",
                "value": round(syn[syn.arm == arm].detected_strict.mean(), 4),
            }
        )
        pf = syn[syn.arm == arm]
        pf_ok = (pf.rep_period - pf.true_period).abs() <= 1
        lines.append(
            {
                "metric": "period_found",
                "arm": arm,
                "key": "synthetic_any_coverage",
                "value": round(float(pf_ok.mean()), 4),
            }
        )
        pfb = pf.assign(ok=pf_ok).groupby("ibin").ok.mean()
        a50p, b50p = floor_at(pfb, 0.5)
        lines.append(
            {
                "metric": "period_found_floor",
                "arm": arm,
                "key": "rate_50",
                "value": b50p if b50p is not None else "never_below",
            }
        )

    for arm in ARMS:
        d = syn[(syn.arm == arm) & (syn.detected == 1)]
        v = [period_verdict(p, tp) for p, tp in zip(d.rep_period, d.true_period, strict=False)]
        n = max(len(v), 1)
        for k in ("exact", "multiple", "divisor", "wrong"):
            lines.append(
                {
                    "metric": "period_verdict",
                    "arm": arm,
                    "key": k,
                    "value": round(v.count(k) / n, 4),
                }
            )
        bias = (d.rep_copies - d.true_copies.astype(float)).astype(float)
        lines.append(
            {
                "metric": "copy_bias_mean",
                "arm": arm,
                "key": "synthetic_detected",
                "value": round(float(bias.mean()), 3),
            }
        )
        lines.append(
            {
                "metric": "copy_bias_median",
                "arm": arm,
                "key": "synthetic_detected",
                "value": round(float(bias.median()), 3),
            }
        )
        lines.append(
            {
                "metric": "copy_abs_err_mean",
                "arm": arm,
                "key": "synthetic_detected",
                "value": round(float(bias.abs().mean()), 3),
            }
        )

    neg = calls[calls.part == "negative"]
    for arm in ARMS:
        for kind, g in neg[neg.arm == arm].groupby("donor_class"):
            lines.append(
                {
                    "metric": "false_positive_rate",
                    "arm": arm,
                    "key": kind,
                    "value": round(float(g.detected.mean()), 4),
                }
            )
        lines.append(
            {
                "metric": "false_positive_rate",
                "arm": arm,
                "key": "all",
                "value": round(float(neg[neg.arm == arm].detected.mean()), 4),
            }
        )

    sow = calls[calls.part == "sowgp"]
    for arm in ARMS:
        g = sow[sow.arm == arm]
        lines.append(
            {
                "metric": "sowgp",
                "arm": arm,
                "key": "recall",
                "value": round(float(g.detected.mean()), 4),
            }
        )
        lines.append(
            {
                "metric": "sowgp",
                "arm": arm,
                "key": "period_47_frac",
                "value": round(float((g.rep_period == 47).mean()), 4),
            }
        )
        gd = g[g.detected == 1]
        b = (gd.rep_copies - gd.true_copies.astype(float)).astype(float)
        lines.append(
            {
                "metric": "sowgp",
                "arm": arm,
                "key": "copy_bias_mean_detected",
                "value": round(float(b.mean()), 3),
            }
        )
        lines.append(
            {
                "metric": "sowgp",
                "arm": arm,
                "key": "copy_abs_err_mean_detected",
                "value": round(float(b.abs().mean()), 3),
            }
        )
        lines.append(
            {
                "metric": "sowgp",
                "arm": arm,
                "key": "copy_exact_frac_detected",
                "value": round(float((b == 0).mean()), 4),
            }
        )
        lines.append({"metric": "sowgp", "arm": arm, "key": "n_detected", "value": int(len(gd))})

    rc = calls[calls.part == "real_candidate"]
    for arm in ARMS:
        lines.append(
            {
                "metric": "real_candidate",
                "arm": arm,
                "key": "n_still_called",
                "value": int(rc[rc.arm == arm].detected.sum()),
            }
        )

    m = pd.DataFrame(lines)
    m.to_csv(out_tsv, sep="\t", index=False)
    return m, rec, syn, neg


def figure(rec, syn, neg, out_png):
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), facecolor=SURFACE)
    for a in ax.ravel():
        a.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            a.spines[sp].set_color(GRID)
        a.tick_params(colors=MUTED, labelsize=9)
        a.grid(axis="y", color=GRID, lw=0.8, zorder=0)
        a.set_axisbelow(True)

    # (a) recall vs measured unit identity
    a = ax[0, 0]
    dy = {"old": 10, "new-exact": 0, "new": -10}
    for arm in ARMS:
        r = rec[arm].dropna()
        a.plot(
            r.index * 100, r.values * 100, lw=2, color=ARM_COLOR[arm], marker="o", ms=4, zorder=3
        )
        a.annotate(
            arm,
            (r.index[-1] * 100, r.values[-1] * 100),
            xytext=(6, dy[arm]),
            textcoords="offset points",
            color=ARM_COLOR[arm],
            fontsize=9,
            va="center",
        )
    a.axhline(50, color=MUTED, lw=1, ls=":")
    a.axhline(90, color=MUTED, lw=1, ls=":")
    a.set_xlabel("measured mean unit-to-unit identity (%)", color=INK, fontsize=10)
    a.set_ylabel("detected (%)", color=INK, fontsize=10)
    a.set_title("(a) the divergence floor", color=INK, fontsize=11, loc="left")
    a.set_ylim(-3, 103)

    # (b) copy-number error vs identity
    a = ax[0, 1]
    syn = syn.copy()
    syn["err"] = syn.rep_copies - syn.true_copies.astype(float)
    for arm in ARMS:
        d = syn[(syn.arm == arm) & (syn.detected == 1)]
        g = d.groupby("ibin").err.mean()
        a.plot(g.index * 100, g.values, lw=2, color=ARM_COLOR[arm], marker="o", ms=4, zorder=3)
        if len(g):
            a.annotate(
                arm,
                (g.index[-1] * 100, g.values[-1]),
                xytext=(6, dy[arm] / 8),
                textcoords="offset points",
                color=ARM_COLOR[arm],
                fontsize=9,
                va="center",
            )
    a.axhline(0, color=MUTED, lw=1)
    a.set_xlabel("measured mean unit-to-unit identity (%)", color=INK, fontsize=10)
    a.set_ylabel("reported − true copies", color=INK, fontsize=10)
    a.set_title("(b) copy-number bias, detected arrays only", color=INK, fontsize=11, loc="left")

    # (c) recall by true period, at identity >= 60%
    a = ax[1, 0]
    hi = syn[syn.true_identity.astype(float) >= 0.6]
    piv = hi.groupby(["arm", "true_period"]).detected.mean().unstack(0)
    x = np.arange(len(piv.index))
    w = 0.27
    for k, arm in enumerate(ARMS):
        a.bar(
            x + (k - 1) * w,
            piv[arm].values * 100,
            w - 0.03,
            color=ARM_COLOR[arm],
            label=arm,
            zorder=3,
            edgecolor=SURFACE,
            linewidth=2,
        )
    a.set_xticks(x)
    a.set_xticklabels(piv.index, fontsize=9)
    a.set_xlabel("true period (aa)", color=INK, fontsize=10)
    a.set_ylabel("detected (%)", color=INK, fontsize=10)
    a.set_title("(c) recall by period, unit identity ≥ 60%", color=INK, fontsize=11, loc="left")
    a.set_ylim(0, 105)
    a.legend(
        frameon=False,
        fontsize=9,
        labelcolor=INK,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
    )

    # (d) false positives on the negatives
    a = ax[1, 1]
    kinds = sorted(neg.donor_class.unique())
    x = np.arange(len(kinds))
    for k, arm in enumerate(ARMS):
        vals = [neg[(neg.arm == arm) & (neg.donor_class == c)].detected.mean() * 100 for c in kinds]
        a.bar(
            x + (k - 1) * w,
            vals,
            w - 0.03,
            color=ARM_COLOR[arm],
            label=arm,
            zorder=3,
            edgecolor=SURFACE,
            linewidth=2,
        )
        n_of = {c: int((neg.donor_class == c).sum() / len(ARMS)) for c in kinds}
        for xi, v, c in zip(x + (k - 1) * w, vals, kinds, strict=False):
            a.annotate(
                f"{round(v * n_of[c] / 100)}/{n_of[c]}",
                (xi, max(v, 0)),
                xytext=(0, 4),
                textcoords="offset points",
                ha="center",
                fontsize=7.5,
                color=INK,
                rotation=90,
            )
    a.set_xticks(x)
    a.set_xticklabels(kinds, fontsize=9)
    a.set_ylabel("false-positive rate (%)", color=INK, fontsize=10)
    a.set_ylim(0, 8)
    a.set_title("(d) false positives on non-repeat controls", color=INK, fontsize=11, loc="left")
    a.legend(frameon=False, fontsize=9, labelcolor=INK, ncol=3, loc="upper right")

    fig.suptitle(
        "Tandem-repeat detector: divergence floor, copy-number bias and false positives",
        color=INK,
        fontsize=13,
        x=0.01,
        ha="left",
    )
    fig.text(
        0.01,
        0.005,
        "2880 synthetic arrays (BLOSUM62-informed divergence model, no indels), 800 non-repeat controls. "
        "Detection = period > 0, coverage ≥ 0.25, copies ≥ 2.5.",
        color=MUTED,
        fontsize=8,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    fig.savefig(out_png.replace(".png", ".pdf"), facecolor=SURFACE)
    print(f"wrote {out_png}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--prefix", default=str(HERE / "repeat_benchmark"))
    args = ap.parse_args()
    truth = pd.read_csv(f"{args.prefix}_truth.tsv", sep="\t")
    calls = run(f"{args.prefix}.fa", truth)
    calls.to_csv(f"{args.prefix}_calls.tsv", sep="\t", index=False)
    m, rec, syn, neg = summarize(calls, f"{args.prefix}_metrics.tsv")
    plots = Path(args.prefix).parent / "plots"
    plots.mkdir(exist_ok=True)
    figure(rec, syn, neg, str(plots / Path(args.prefix).name) + ".png")
    with pd.option_context("display.max_rows", 200, "display.width", 140):
        print(m.to_string(index=False))


if __name__ == "__main__":
    main()
