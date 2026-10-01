"""Fit and score the Phase C candidates on one outer training set (Phase C spec 3.3).

Nested protocol (one rule for every fitted quantity). For one outer training set:
1. Inner folds: StratifiedGroupKFold(3, shuffle, seed) over the training rows, groups =
   MMseqs2 clusters. Each inner test fold must hold both classes, else ModelError.
2. Logistic-regression candidates (B0, B1, M8, M35, M8-C, M35-C; pipeline StandardScaler +
   LogisticRegression(class_weight="balanced")): for each C in C_GRID, the inner out-of-fold
   decision values; C = the value with the highest inner out-of-fold PR-AUC (ties: the smaller
   C). H: the same over (ESM variant, C); ties: the first variant in M8, M35, M8-C, M35-C order.
3. With the chosen setting: the decision threshold maximises Youden's J on the inner
   out-of-fold values (ties: the higher threshold); Platt scaling (a, b) is fitted on the same
   values (ruling C-10). The final pipeline is refitted on all training rows.
4. Rules: g and t by Youden's J on the training rows (rules.fit_rule).
Only the training rows' labels are read. Scores are decision values; `prob` is the Platt
probability; `call` is score >= threshold.
"""

import warnings

import numpy as np
import rules
import universe
from metrics import pr_auc
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

CANDIDATES = ("B0", "B1", "R0", "R1", "R2", "M8", "M35", "M8-C", "M35-C", "H")
LR_CANDIDATES = ("B0", "B1", "M8", "M35", "M8-C", "M35-C")
ML_CANDIDATES = ("M8", "M35", "M8-C", "M35-C", "H")
H_VARIANTS = ("M8", "M35", "M8-C", "M35-C")
C_GRID = (0.001, 0.003, 0.01, 0.1, 1.0, 10.0)  # ruling C-12: S1 fold 0 chose 0.01, the old edge
INNER_FOLDS = 3
MAX_ITER = 5000
SCORE_CHUNK = 20000


class ModelError(ValueError):
    """An outer training set cannot be fitted (for example an inner fold without positives)."""


def inner_folds(y, groups, seed: int, n_folds: int = INNER_FOLDS):
    y = np.asarray(y, dtype=bool)
    groups = np.asarray(groups).astype(str)
    n_pos, n_neg = len(set(groups[y])), len(set(groups[~y]))
    if n_pos < n_folds or n_neg < n_folds:
        raise ModelError(
            f"{n_pos} positive and {n_neg} negative clusters; the inner {n_folds}-fold split "
            "needs at least one cluster of each class per fold"
        )
    cv = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    folds = list(cv.split(np.zeros(len(y)), y, groups))
    for k, (_, test) in enumerate(folds):
        if y[test].all() or not y[test].any():
            raise ModelError(
                f"inner fold {k} has no {'negatives' if y[test].all() else 'positives'}"
            )
    return folds


def make_lr(C: float):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(C=C, class_weight="balanced", max_iter=MAX_ITER),
    )


def _fit(model, X, y, notes: dict):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(X, y)
    notes["convergence_warnings"] = notes.get("convergence_warnings", 0) + sum(
        issubclass(w.category, ConvergenceWarning) for w in caught
    )
    return model


def oof_scores(X, y, folds, C: float, notes: dict) -> np.ndarray:
    out = np.empty(len(y))
    for train, test in folds:
        out[test] = _fit(make_lr(C), X[train], y[train], notes).decision_function(X[test])
    return out


def youden_threshold(y, score) -> float:
    """The score with the highest recall - FPR when calls are score >= threshold."""
    y = np.asarray(y, dtype=bool)
    score = np.asarray(score, dtype=np.float64)
    order = np.argsort(-score, kind="stable")
    s = score[order]
    ends = np.r_[np.nonzero(s[1:] != s[:-1])[0], len(s) - 1]
    tp = np.cumsum(y[order])[ends] / y.sum()
    fp = np.cumsum(~y[order])[ends] / (~y).sum()
    return float(s[ends][int(np.argmax(tp - fp))])


def fit_platt(y, score) -> tuple[float, float]:
    """Platt scaling: p = 1 / (1 + exp(-(a * score + b))), fitted without class weights."""
    lr = LogisticRegression(C=1e6, max_iter=MAX_ITER)
    lr.fit(np.asarray(score, dtype=np.float64)[:, None], np.asarray(y, dtype=bool))
    return float(lr.coef_[0, 0]), float(lr.intercept_[0])


def platt(score, a: float, b: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(a * np.asarray(score, dtype=np.float64) + b)))


def _score(model, u, candidate, idx, h_variant) -> np.ndarray:
    out = np.empty(len(idx))
    for start in range(0, len(idx), SCORE_CHUNK):
        part = idx[start : start + SCORE_CHUNK]
        out[start : start + len(part)] = model.decision_function(
            universe.features(u, candidate, part, h_variant)
        )
    return out


def fit_unit(u, train_idx, y, groups, score_idx, seed: int, candidates=CANDIDATES):
    """Fit every candidate on the training rows; score the rows score_idx.

    Return (params, scores): params[candidate] is a JSON-ready dict; scores[candidate] is a
    dict with `score` and `prob` (arrays, or None for rules) and `call` (bool array)."""
    train_idx = np.asarray(train_idx, dtype=np.int64)
    score_idx = np.asarray(score_idx, dtype=np.int64)
    y = np.asarray(y, dtype=bool)
    folds = inner_folds(y, groups, seed)
    params, scores = {}, {}
    base = {"n_train": int(len(y)), "n_train_pos": int(y.sum())}
    for cand in candidates:
        if cand in rules.RULES:
            fit = rules.fit_rule(cand, y, u.sp[train_idx], u.rank[train_idx], u.st[train_idx])
            call = rules.rule_call(
                cand, u.sp[score_idx], u.rank[score_idx], u.st[score_idx], fit["g"], fit["t"]
            )
            params[cand] = {**base, **fit}
            scores[cand] = {"score": None, "prob": None, "call": call}
            continue
        notes: dict = {}
        settings = (
            [(v, c) for v in H_VARIANTS for c in C_GRID]
            if cand == "H"
            else [(None, c) for c in C_GRID]
        )
        best, inner = None, {}
        cache = {}
        for variant, C in settings:
            if variant not in cache:
                cache[variant] = universe.features(u, cand, train_idx, variant)
            oof = oof_scores(cache[variant], y, folds, C, notes)
            ap = float(pr_auc(np.ones(len(y)), y, oof)[0])
            inner[f"{variant}:{C}" if variant else f"{C}"] = ap
            if best is None or ap > best[0]:
                best = (ap, variant, C, oof)
        ap, variant, C, oof = best
        threshold = youden_threshold(y, oof)
        a, b = fit_platt(y, oof)
        model = _fit(make_lr(C), cache[variant], y, notes)
        s = _score(model, u, cand, score_idx, variant)
        params[cand] = {
            **base,
            "C": C,
            "h_variant": variant,
            "threshold": threshold,
            "platt_a": a,
            "platt_b": b,
            "inner_pr_auc": inner,
            "convergence_warnings": notes.get("convergence_warnings", 0),
        }
        scores[cand] = {"score": s, "prob": platt(s, a, b), "call": s >= threshold}
    return params, scores


# --- work units for 10_fit_and_score.py (module level, so worker processes can import them) ---

_WORKER: dict = {}


def init_worker(phaseb: str) -> None:
    """Process-pool initializer: load the universe once per worker (10 verified the files)."""
    _WORKER["u"] = universe.load(phaseb, verify=False)


def run_task(task: dict, u=None):
    """Fit one outer training set. task: key, train, y, groups, score, seed, candidates."""
    u = u if u is not None else _WORKER["u"]
    try:
        params, scores = fit_unit(
            u,
            u.rows(task["train"]),
            task["y"],
            task["groups"],
            u.rows(task["score"]),
            task["seed"],
            task["candidates"],
        )
    except ModelError as exc:
        raise ModelError(f"{task['key']}: {exc}") from exc
    return task["key"], params, scores
