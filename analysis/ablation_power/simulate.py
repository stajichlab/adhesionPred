"""Planted-signal simulation: what improvement could this benchmark actually detect?

Three signal shapes, because they answer different questions and only the third is actionable:

  homogeneous    a Gaussian feature with the same effect in every clade. The most favourable case
                 that could ever occur, so its curve is an optimistic bound.
  heterogeneous  full effect in one clade, attenuated elsewhere. Leave-one-clade-out cannot see a
                 signal that does not transfer -- the model scoring a held-out clade never saw it --
                 so this is the realistic case for a domain family with narrow taxonomic range.
  binary         a domain call with a stated sensitivity on positives and false-positive rate on
                 hard negatives. This is the only parameterisation a biologist can act on:
                 "a model that hits X% of adhesins at Y% FPR".

Two things the review showed must be measured rather than assumed:

  TYPE-I RATE. The decision rule's real rejection rate under a pure-noise component, at the nominal
  level. A design whose nominal 5% is really 9.5% turns six component tests into a coin flip.
  DIMENSIONALITY BIAS. Pure noise is not worth zero. Measured on adhesion_v1, one noise feature
  gains a mean +0.007 AP and twenty loses 0.025, so a component is rewarded or penalised by its
  feature count before it says anything about biology. The null must therefore be a permuted
  version of THAT component, not the number zero.
"""

import numpy as np
from model import loco_scores
from stats import paired_permutation_test, stratified_concordance


class SimulationError(Exception):
    pass


def plant_homogeneous(y, effect, rng, dims=1):
    """Gaussian feature with mean separation `effect` (Cohen's d) between classes, every clade."""
    y = np.asarray(y)
    x = rng.standard_normal((y.size, dims))
    x[:, 0] += effect * y
    return x


def plant_heterogeneous(y, fold, effect, rng, strong_fold, attenuation=1 / 3, dims=1):
    """Full effect inside `strong_fold`, attenuated elsewhere: a signal that transfers poorly."""
    y = np.asarray(y)
    fold = np.asarray(fold)
    scale = np.where(fold == strong_fold, 1.0, attenuation)
    x = rng.standard_normal((y.size, dims))
    x[:, 0] += effect * scale * y
    return x


def plant_binary(y, sensitivity, fpr, rng):
    """A domain call: present with probability `sensitivity` on positives, `fpr` on negatives."""
    y = np.asarray(y)
    probability = np.where(y == 1, sensitivity, fpr)
    return (rng.random(y.size) < probability).astype(float)[:, None]


def _delta(bench, bar_name, bar_scores, component, scored_mask, c, seed):
    """Stratified-concordance difference, component-augmented minus bar."""
    with_component, _ = loco_scores(bench, bar_name, component=component, c=c, seed=seed)
    y = bench.y[scored_mask]
    fold = np.array(bench.fold)[scored_mask]
    return (
        stratified_concordance(y, with_component[scored_mask], fold)
        - stratified_concordance(y, bar_scores[scored_mask], fold),
        with_component,
    )


def one_replicate(bench, bar_name, bar_scores, component, c=1.0, n_perm=1000, seed=0, alpha=0.05):
    """Run one simulated component through the real decision rule. Returns (delta, p, rejected)."""
    scored_mask = np.array([f in set(bench.scored_folds) for f in bench.fold])
    delta, with_component = _delta(bench, bar_name, bar_scores, component, scored_mask, c, seed)
    if not np.isfinite(delta):
        return float("nan"), float("nan"), False
    _, p = paired_permutation_test(
        bench.y[scored_mask],
        bar_scores[scored_mask],
        with_component[scored_mask],
        np.array(bench.cluster)[scored_mask],
        np.array(bench.fold)[scored_mask],
        n=n_perm,
        seed=seed,
        alternative="greater",
    )
    return float(delta), float(p), bool(p <= alpha)


# --------------------------------------------------------------------------- parallel execution
#
# The work is embarrassingly parallel: every (shape, effect, replicate) is independent. Measured
# serially on the HPCC login node, one training pool took over 15 minutes; the whole grid at the
# replicate count the result deserves is not practical in one process.
#
# Each task carries its OWN seed and builds its own component from it. That is not only for
# parallelism: the serial version shared one generator across the replicates of a grid point, so a
# result depended on the order they ran in. Per-task seeding makes a run reproducible whatever the
# worker count, which is the property a number quoted in a report needs.

_WORKER = {}


def _init_worker(bench, bar_name, bar_scores, bar_c, n_perm, alpha):
    _WORKER.update(
        bench=bench,
        bar_name=bar_name,
        bar_scores=bar_scores,
        bar_c=bar_c,
        n_perm=n_perm,
        alpha=alpha,
    )


def build_component(bench, task):
    """Make one task's component from its own seed. Deterministic in (shape, effect, replicate)."""
    rng = np.random.default_rng(task["seed"])
    shape = task["shape"]
    if shape == "noise":
        return rng.standard_normal((len(bench), task["dims"]))
    if shape == "homogeneous":
        return plant_homogeneous(bench.y, task["effect"], rng)
    if shape == "heterogeneous":
        if not task.get("strong_fold"):
            raise SimulationError("heterogeneous needs strong_fold")
        return plant_heterogeneous(
            bench.y,
            bench.fold,
            task["effect"],
            rng,
            task["strong_fold"],
            task.get("attenuation", 1 / 3),
        )
    if shape == "binary":
        return plant_binary(bench.y, task["effect"], task.get("fpr", 0.05), rng)
    raise SimulationError(f"unknown shape {shape!r}")


def _run_one(bench, bar_name, bar_scores, bar_c, n_perm, alpha, task):
    component = build_component(bench, task)
    delta, p, rejected = one_replicate(
        bench,
        bar_name,
        bar_scores,
        component,
        c=bar_c,
        n_perm=n_perm,
        seed=task["seed"],
        alpha=alpha,
    )
    return {**task, "delta": delta, "p": p, "rejected": bool(rejected)}


def _run_one_in_worker(task):
    return _run_one(
        _WORKER["bench"],
        _WORKER["bar_name"],
        _WORKER["bar_scores"],
        _WORKER["bar_c"],
        _WORKER["n_perm"],
        _WORKER["alpha"],
        task,
    )


def run_tasks(bench, bar_name, bar_scores, tasks, bar_c=1.0, n_perm=400, alpha=0.05, workers=1):
    """Run every task, in one process when workers == 1 and in a pool otherwise.

    The two paths must agree exactly; a test pins that.
    """
    if workers == 1:
        return [_run_one(bench, bar_name, bar_scores, bar_c, n_perm, alpha, t) for t in tasks]

    from concurrent.futures import ProcessPoolExecutor

    chunk = max(1, len(tasks) // (workers * 8))
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_worker,
        initargs=(bench, bar_name, bar_scores, bar_c, n_perm, alpha),
    ) as pool:
        return list(pool.map(_run_one_in_worker, tasks, chunksize=chunk))


def _summarise(rows, keys):
    deltas = [r["delta"] for r in rows if np.isfinite(r["delta"])]
    out = {k: rows[0][k] for k in keys}
    out.update(
        replicates=len(rows),
        power=sum(r["rejected"] for r in rows) / len(rows),
        mean_delta=float(np.mean(deltas)) if deltas else float("nan"),
    )
    return out


def noise_tasks(dims_list, replicates, seed):
    return [
        {
            "shape": "noise",
            "dims": int(d),
            "effect": 0.0,
            "replicate": r,
            "seed": seed + 1000 * d + r,
        }
        for d in dims_list
        for r in range(replicates)
    ]


def signal_tasks(shape, grid, replicates, seed, **extra):
    tasks = []
    for i, point in enumerate(grid):
        for r in range(replicates):
            tasks.append(
                {
                    "shape": shape,
                    "effect": float(point),
                    "replicate": r,
                    "seed": seed + 100_000 * (i + 1) + r,
                    **extra,
                }
            )
    return tasks


def type_i_rate(
    bench,
    bar_name,
    bar_scores,
    dims,
    replicates=200,
    c=1.0,
    n_perm=400,
    seed=0,
    alpha=0.05,
    workers=1,
):
    """Rejection rate under a pure-noise component of a given dimensionality.

    Under the null this should equal `alpha`. The review measured 9.5% at a nominal 5% for the
    original design, so this number is the check that the decision rule means what it says.
    """
    tasks = noise_tasks([dims], replicates, seed)
    rows = run_tasks(
        bench, bar_name, bar_scores, tasks, bar_c=c, n_perm=n_perm, alpha=alpha, workers=workers
    )
    deltas = [r["delta"] for r in rows if np.isfinite(r["delta"])]
    return {
        "dims": int(dims),
        "replicates": int(replicates),
        "alpha": float(alpha),
        "rejection_rate": sum(r["rejected"] for r in rows) / len(rows),
        "mean_delta_under_null": float(np.mean(deltas)) if deltas else float("nan"),
        "fraction_delta_positive": float(np.mean(np.array(deltas) > 0)) if deltas else float("nan"),
    }


def power_curve(
    bench,
    bar_name,
    bar_scores,
    shape,
    grid,
    replicates=200,
    c=1.0,
    n_perm=400,
    seed=0,
    alpha=0.05,
    strong_fold=None,
    attenuation=1 / 3,
    fpr=0.05,
    workers=1,
):
    """Rejection rate against effect size, for one signal shape.

    `grid` is a list of effect sizes: Cohen's d for the Gaussian shapes, sensitivity for `binary`.
    """
    if shape not in ("homogeneous", "heterogeneous", "binary"):
        raise SimulationError(f"unknown shape {shape!r}")
    if shape == "heterogeneous" and strong_fold is None:
        raise SimulationError("heterogeneous needs strong_fold")
    extra = {"attenuation": attenuation, "fpr": fpr}
    if strong_fold is not None:
        extra["strong_fold"] = strong_fold
    tasks = signal_tasks(shape, grid, replicates, seed, **extra)
    rows = run_tasks(
        bench, bar_name, bar_scores, tasks, bar_c=c, n_perm=n_perm, alpha=alpha, workers=workers
    )
    by_point = {}
    for row in rows:
        by_point.setdefault(row["effect"], []).append(row)
    return [_summarise(by_point[e], ("shape", "effect")) for e in sorted(by_point)]


def minimum_detectable_effect(rows, target_power=0.8):
    """Smallest effect whose measured power reaches `target_power`, linearly interpolated.

    Returns None when the grid never reaches the target -- which is itself the finding.
    """
    ordered = sorted(rows, key=lambda r: r["effect"])
    for previous, current in zip(ordered, ordered[1:], strict=False):
        if previous["power"] < target_power <= current["power"]:
            span = current["power"] - previous["power"]
            if span == 0:
                return current["effect"], current["mean_delta"]
            frac = (target_power - previous["power"]) / span
            effect = previous["effect"] + frac * (current["effect"] - previous["effect"])
            delta = previous["mean_delta"] + frac * (current["mean_delta"] - previous["mean_delta"])
            return float(effect), float(delta)
    if ordered and ordered[0]["power"] >= target_power:
        return float(ordered[0]["effect"]), float(ordered[0]["mean_delta"])
    return None
