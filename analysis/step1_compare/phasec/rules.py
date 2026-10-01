"""The rule candidates R0, R1, R2 (Phase C spec 3.3, rulings C-2, C-3). Needs numpy.

R0: SignalP call is SP.
R1: SP and PredGPI class at or above g.
R2: SP and (PredGPI class at or above g, or ser_thr_frac >= t). R2 is "the rule".

g takes highly_probable, probable, weakly (PredGPI scores 1.0, 0.70, 0.55); t takes 0.20 to 0.40
in steps of 0.05 (owner decision 2026-10-01: at t = 0.10 R2 equals R0 on the real data, because
96.7% of the SP proteins have ser_thr_frac >= 0.10). g and t are fitted by Youden's J (recall minus FPR) on the training rows of
one outer training set. A rule has no fitted model, so J on the training rows equals J on the
pooled inner out-of-fold rows. Ties in J: the first setting in grid order wins (g from
highly_probable down, then t from 0.40 down), that is the stricter rule.
"""

import numpy as np

GPI_RANK = {"highly_probable": 3, "probable": 2, "weakly": 1, "none": 0, "too_short": 0}
G_VALUES = ("highly_probable", "probable", "weakly")
T_VALUES = (0.20, 0.25, 0.30, 0.35, 0.40)  # literals: 0.2 + 0.05 * 2 != 0.30
RULES = ("R0", "R1", "R2")


def gpi_rank(calls) -> np.ndarray:
    try:
        return np.array([GPI_RANK[c] for c in calls], dtype=np.int8)
    except KeyError as exc:
        raise ValueError(f"unknown gpi_call {exc.args[0]!r}") from exc


def rule_call(rule: str, sp, rank, st, g: str | None = None, t: float | None = None):
    """Boolean call of `rule` for arrays sp (bool), rank (gpi_rank), st (ser_thr_frac)."""
    sp = np.asarray(sp, dtype=bool)
    if rule == "R0":
        return sp.copy()
    gpi = np.asarray(rank) >= GPI_RANK[g]
    if rule == "R1":
        return sp & gpi
    if rule == "R2":
        return sp & (gpi | (np.asarray(st, dtype=np.float64) >= t))
    raise ValueError(f"unknown rule {rule!r}")


def youden(y, call) -> float:
    """Recall minus FPR; NaN when y has no positives or no negatives."""
    y = np.asarray(y, dtype=bool)
    call = np.asarray(call, dtype=bool)
    p, n = y.sum(), (~y).sum()
    if p == 0 or n == 0:
        return float("nan")
    return float((y & call).sum() / p - (~y & call).sum() / n)


def fit_rule(rule: str, y, sp, rank, st) -> dict:
    """g and t with the highest Youden's J on these rows. R0 has no parameter."""
    if rule == "R0":
        return {"g": None, "t": None, "j": youden(y, rule_call("R0", sp, rank, st))}
    best = None
    t_grid = T_VALUES[::-1] if rule == "R2" else (None,)
    for g in G_VALUES:
        for t in t_grid:
            j = youden(y, rule_call(rule, sp, rank, st, g, t))
            if np.isnan(j):
                raise ValueError("Youden's J is not defined: the rows lack a class")
            if best is None or j > best["j"]:
                best = {"g": g, "t": t, "j": j}
    return best
