"""Power calibration for the adhesion ablation harness: what could this benchmark detect?

This answers one question, before any component is built or tested:

    How large an improvement would a component have to make before `adhesion_v1` could tell?

If that minimum detectable effect is large, the finding is that the benchmark is too small, and
curating more truth outranks building the harness. That decision is the point of this script.

It deliberately tests NO real component. It plants synthetic signals of known strength, runs them
through the real decision rule, and reports the rejection rate.

Run (a few minutes; HPCC login node is fine -- one core, no SLURM needed):
    python3 analysis/ablation_power/run_power.py --benchmark data/benchmark/adhesion_v1 \
        --out /path/to/NEW_DIR

NOT `python3 -I`, unlike the stdlib-only scripts elsewhere in this repo. Isolated mode implies -s,
which drops the user site-packages -- and on the HPCC numpy, scipy and scikit-learn live in
~/.local/lib/python3.12/site-packages. With -I this script dies on `import numpy`. Because the
result therefore depends on packages outside the repo, their versions are recorded in
baseline.json alongside the benchmark digest.

Writes baseline.json, type_i.tsv, power.tsv, summary.md. The output directory must not exist.
"""

import argparse
import json
import platform
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from benchmark_data import TRAIN_POOLS, assert_no_leakage, load, splits  # noqa: E402
from model import choose_bar  # noqa: E402
from simulate import minimum_detectable_effect, power_curve, type_i_rate  # noqa: E402
from stats import per_fold_report  # noqa: E402

D_GRID = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0)
SENSITIVITY_GRID = (0.05, 0.15, 0.25, 0.35, 0.45, 0.60, 0.80)
NOISE_DIMS = (1, 5, 20)
TARGET_POWER = 0.8


def fail(message):
    print(f"FATAL: {message}", file=sys.stderr)
    raise SystemExit(1)


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write(
                "\t".join("" if row.get(c) is None else str(row.get(c, "")) for c in columns) + "\n"
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--train-pool", choices=TRAIN_POOLS, default="all")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--permutations", type=int, default=400)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--attenuation", type=float, default=1 / 3)
    parser.add_argument("--binary-fpr", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--no-verify", action="store_true", help="skip the SHA256SUMS check (tests only)"
    )
    args = parser.parse_args(argv)

    out = Path(args.out)
    if out.exists():
        fail(f"{out} exists; nothing is overwritten")

    bench = load(args.benchmark, train_pool=args.train_pool, verify=not args.no_verify)
    assert_no_leakage(bench)

    scored_mask = np.array([f in set(bench.scored_folds) for f in bench.fold])
    if not scored_mask.any():
        fail("no scorable fold in this benchmark")

    # ---- the bar, and the split table, both recorded before anything is simulated
    bar_name, bar_scores, family, bar_c = choose_bar(bench, tune=True, seed=args.seed)
    split_rows = []
    for held, train, test, dropped in splits(bench):
        split_rows.append(
            {
                "fold": held,
                "test_positive": int(bench.y[test].sum()),
                "test_negative": int((bench.y[test] == 0).sum()),
                "train_positive": int(bench.y[train].sum()),
                "train_negative": int((bench.y[train] == 0).sum()),
                "dropped_for_leakage": dropped,
            }
        )

    effective = _effective_sample_size(bench, scored_mask)
    baseline = {
        "benchmark_version": bench.version,
        "benchmark_sha256sums_digest": bench.digest,
        "train_pool": bench.train_pool,
        "bar": bar_name,
        "bar_C_per_fold": bar_c,
        "baseline_family": family,
        "splits": split_rows,
        "per_fold": per_fold_report(
            bench.y[scored_mask], bar_scores[scored_mask], np.array(bench.fold)[scored_mask]
        ),
        "pooled_counts": {
            "positives": int(bench.y[scored_mask].sum()),
            "negatives": int((bench.y[scored_mask] == 0).sum()),
            "clusters": int(len(set(np.array(bench.cluster)[scored_mask]))),
        },
        "effective_sample_size": effective,
        "settings": {
            "replicates": args.replicates,
            "permutations": args.permutations,
            "alpha": args.alpha,
            "attenuation": args.attenuation,
            "binary_fpr": args.binary_fpr,
            "seed": args.seed,
            "python": platform.python_version(),
            "packages": _package_versions(),
        },
    }

    # ---- type-I rate, per component dimensionality
    type_i = [
        type_i_rate(
            bench,
            bar_name,
            bar_scores,
            dims,
            c=bar_c,
            replicates=args.replicates,
            n_perm=args.permutations,
            seed=args.seed,
            alpha=args.alpha,
        )
        for dims in NOISE_DIMS
    ]

    # ---- power curves
    strong = max(
        bench.scored_folds,
        key=lambda f: int(np.sum(np.array(bench.fold) == f)),
    )
    curves = []
    curves += power_curve(
        bench,
        bar_name,
        bar_scores,
        "homogeneous",
        D_GRID,
        c=bar_c,
        replicates=args.replicates,
        n_perm=args.permutations,
        seed=args.seed,
        alpha=args.alpha,
    )
    curves += power_curve(
        bench,
        bar_name,
        bar_scores,
        "heterogeneous",
        D_GRID,
        c=bar_c,
        replicates=args.replicates,
        n_perm=args.permutations,
        seed=args.seed,
        alpha=args.alpha,
        strong_fold=strong,
        attenuation=args.attenuation,
    )
    curves += power_curve(
        bench,
        bar_name,
        bar_scores,
        "binary",
        SENSITIVITY_GRID,
        c=bar_c,
        replicates=args.replicates,
        n_perm=args.permutations,
        seed=args.seed,
        alpha=args.alpha,
        fpr=args.binary_fpr,
    )

    mde = {}
    for shape in ("homogeneous", "heterogeneous", "binary"):
        rows = [r for r in curves if r["shape"] == shape]
        found = minimum_detectable_effect(rows, TARGET_POWER)
        mde[shape] = None if found is None else {"effect": found[0], "delta_concordance": found[1]}

    out.mkdir(parents=True)
    (out / "baseline.json").write_text(json.dumps(baseline, indent=2, sort_keys=True) + "\n")
    write_tsv(
        out / "type_i.tsv",
        type_i,
        [
            "dims",
            "replicates",
            "alpha",
            "rejection_rate",
            "mean_delta_under_null",
            "fraction_delta_positive",
        ],
    )
    write_tsv(out / "power.tsv", curves, ["shape", "effect", "replicates", "power", "mean_delta"])
    (out / "summary.md").write_text(_summary(baseline, type_i, mde, strong, args))
    print(f"wrote {out}")
    return 0


def _package_versions():
    """Record the versions this result depends on: they are outside the repo and outside the digest."""
    import numpy
    import scipy
    import sklearn

    return {
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "scikit-learn": sklearn.__version__,
    }


def _effective_sample_size(bench, mask):
    """Kish effective sample size, counting a paralogue family once rather than nine times."""
    out = {}
    cluster = np.array(bench.cluster)[mask]
    y = bench.y[mask]
    for label, name in ((1, "positives"), (0, "negatives")):
        sel = y == label
        _, counts = np.unique(cluster[sel], return_counts=True)
        if counts.size == 0:
            out[name] = None
            continue
        out[name] = {
            "n_proteins": int(sel.sum()),
            "n_clusters": int(counts.size),
            "kish_ess": float(counts.sum() ** 2 / np.sum(counts**2)),
        }
    return out


def _summary(baseline, type_i, mde, strong, args):
    pooled = baseline["pooled_counts"]
    ess = baseline["effective_sample_size"]
    lines = [
        f"# Ablation power calibration ({date.today().isoformat()})",
        "",
        f"Benchmark `{baseline['benchmark_version']}`, training pool `{baseline['train_pool']}`, "
        f"bar `{baseline['bar']}` (best of the declared family by inner cluster-grouped CV).",
        "",
        "## What the benchmark actually weighs",
        "",
        f"- {pooled['positives']} positives and {pooled['negatives']} hard negatives in "
        f"{pooled['clusters']} clusters across the scorable folds.",
    ]
    if ess.get("positives"):
        lines.append(
            f"- **Effective sample size {ess['positives']['kish_ess']:.1f} positives** "
            f"({ess['positives']['n_clusters']} clusters), not {ess['positives']['n_proteins']}: "
            "paralogue families are one observation, not many."
        )
    lines += [
        "",
        "## Does the decision rule mean what it says?",
        "",
        f"Rejection rate under a pure-noise component, nominal alpha = {args.alpha}:",
        "",
        "| component dims | rejection rate | mean delta under the null |",
        "|---|---|---|",
    ]
    for row in type_i:
        lines.append(
            f"| {row['dims']} | {row['rejection_rate']:.3f} | {row['mean_delta_under_null']:+.4f} |"
        )
    lines += [
        "",
        "A rejection rate above alpha means the rule is anti-conservative and every later "
        "component test is easier to pass than it claims. A non-zero mean delta means a component "
        "is rewarded or penalised for its feature count alone, which is why the null for a real "
        "component must be a permutation of that component, never the number zero.",
        "",
        "## Minimum detectable effect at 80% power",
        "",
        "| signal | effect | delta in stratified concordance |",
        "|---|---|---|",
    ]
    names = {
        "homogeneous": "Gaussian, same in every clade (optimistic bound)",
        "heterogeneous": f"Gaussian, full in {strong}, {args.attenuation:.2f}x elsewhere (realistic)",
        "binary": f"domain call at {args.binary_fpr:.0%} FPR (sensitivity needed)",
    }
    for shape, label in names.items():
        found = mde.get(shape)
        if found is None:
            lines.append(f"| {label} | **not reached on the grid** | — |")
        else:
            lines.append(f"| {label} | {found['effect']:.2f} | {found['delta_concordance']:+.3f} |")
    lines += [
        "",
        "The binary row is the actionable one: it says what sensitivity a domain model would need, "
        "at the stated false-positive rate on hard negatives, before this benchmark could see it.",
        "",
        "## How to read a null result",
        "",
        "A component that fails to clear this bar has not been shown to be useless. It has been "
        "shown to be undetectable at this benchmark's size. Report the interval's upper bound, not "
        'the phrase "no improvement".',
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
