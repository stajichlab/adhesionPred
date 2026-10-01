"""Weighted classification metrics for Phase C (Phase C spec 4). Needs numpy.

Every function takes a weight matrix W of shape (B, n) (or a vector of length n, read as one
row) and returns one value per row of W. A row of W holds the number of times each protein is
in one bootstrap resample; a row of ones gives the point estimate. A metric that is not defined
for a row (no positives, no negatives, no positive calls) is NaN for that row. Non-finite
scores and scores whose shape differs from the label shape raise ValueError.

Definitions (thresholds are the distinct scores; a protein is called when score >= threshold):
- recall = TP / P; precision = TP / (TP + FP); FPR = FP / N (binary calls).
- ROC-AUC: area under the ROC step curve with ties counted as half (trapezoid per tie group).
- PR-AUC: average precision, sum over thresholds of (R_k - R_k-1) * P_k (sklearn definition).
- precision at recall r: the highest precision at any threshold with recall >= r.
- recall at FPR f: the highest recall at any threshold with FPR <= f (the empty call set, with
  recall 0 and FPR 0, counts as a threshold).
- FPR of a subset at recall r: the FPR inside `neg_mask` at the first (highest) threshold whose
  recall reaches r.
"""

import numpy as np

TOL = 1e-12


def as_weights(W, n: int) -> np.ndarray:
    W = np.asarray(W, dtype=np.float64)
    if W.ndim == 1:
        W = W[None, :]
    if W.ndim != 2 or W.shape[1] != n:
        raise ValueError(f"weights must have shape (B, {n}), got {W.shape}")
    return W


def _div(a, b) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    out = np.full(np.broadcast(a, b).shape, np.nan)
    ok = b > 0
    np.divide(a, b, out=out, where=ok)
    return out


def _bool(x) -> np.ndarray:
    return np.asarray(x, dtype=bool)


def recall(W, y, call) -> np.ndarray:
    y, call = _bool(y), _bool(call)
    W = as_weights(W, len(y))
    return _div(W @ (y & call), W @ y)


def precision(W, y, call) -> np.ndarray:
    y, call = _bool(y), _bool(call)
    W = as_weights(W, len(y))
    return _div(W @ (y & call), W @ call)


def fpr(W, y, call) -> np.ndarray:
    y, call = _bool(y), _bool(call)
    W = as_weights(W, len(y))
    return _div(W @ (~y & call), W @ ~y)


def curve(W, y, score, extra=None):
    """Cumulative weighted TP and FP at each distinct threshold, from the highest score down.

    Return (tp, fp, p, n, extra_cum): tp and fp have shape (B, G) for G distinct scores; p and n
    have shape (B,); extra_cum is the cumulative weight of the rows in the boolean mask `extra`
    (shape (B, G)) or None."""
    y = _bool(y)
    score = np.asarray(score, dtype=np.float64)
    if score.shape != y.shape:
        raise ValueError(f"score shape {score.shape} differs from label shape {y.shape}")
    if not np.isfinite(score).all():
        raise ValueError("scores must be finite")
    W = as_weights(W, len(y))
    order = np.argsort(-score, kind="stable")
    s = score[order]
    ends = np.r_[np.nonzero(s[1:] != s[:-1])[0], len(s) - 1] if len(s) else np.array([], int)
    Ws = W[:, order]
    tp = np.cumsum(Ws * y[order], axis=1)[:, ends]
    fp = np.cumsum(Ws * ~y[order], axis=1)[:, ends]
    ex = None
    if extra is not None:
        ex = np.cumsum(Ws * _bool(extra)[order], axis=1)[:, ends]
    return tp, fp, W @ y, W @ ~y, ex


def _roc(tp, fp, p, n) -> np.ndarray:
    tp0 = np.c_[np.zeros(len(p)), tp]
    fp0 = np.c_[np.zeros(len(p)), fp]
    area = np.sum((fp0[:, 1:] - fp0[:, :-1]) * (tp0[:, 1:] + tp0[:, :-1]) / 2.0, axis=1)
    return _div(area, p * n)


def _prec(tp, fp, empty) -> np.ndarray:
    called = tp + fp
    return np.where(called > 0, tp / np.where(called > 0, called, 1.0), empty)


def _ap(tp, fp, p) -> np.ndarray:
    dtp = np.diff(np.c_[np.zeros(len(p)), tp], axis=1)
    return _div(np.sum(dtp * _prec(tp, fp, 0.0), axis=1), p)


def _prec_at_recall(tp, fp, p, level) -> np.ndarray:
    level = np.broadcast_to(np.asarray(level, dtype=np.float64), p.shape)[:, None]
    prec = _prec(tp, fp, np.nan)
    ok = (_div(tp, p[:, None]) >= level - TOL) & ~np.isnan(prec)
    best = np.where(ok, prec, -np.inf).max(axis=1, initial=-np.inf)
    best[~np.isfinite(best) | ~(p > 0)] = np.nan
    return best


def _recall_at_fpr(tp, fp, p, n, level) -> np.ndarray:
    level = np.broadcast_to(np.asarray(level, dtype=np.float64), p.shape)[:, None]
    rec = np.c_[np.zeros(len(p)), _div(tp, p[:, None])]
    rate = np.c_[np.zeros(len(p)), _div(fp, n[:, None])]
    best = np.where(rate <= level + TOL, rec, -np.inf).max(axis=1)
    best[~((p > 0) & (n > 0))] = np.nan
    return best


def roc_auc(W, y, score) -> np.ndarray:
    tp, fp, p, n, _ = curve(W, y, score)
    return _roc(tp, fp, p, n)


def pr_auc(W, y, score) -> np.ndarray:
    tp, fp, p, _, _ = curve(W, y, score)
    return _ap(tp, fp, p)


def precision_at_recall(W, y, score, level) -> np.ndarray:
    tp, fp, p, _, _ = curve(W, y, score)
    return _prec_at_recall(tp, fp, p, level)


def recall_at_fpr(W, y, score, level) -> np.ndarray:
    tp, fp, p, n, _ = curve(W, y, score)
    return _recall_at_fpr(tp, fp, p, n, level)


def scored_summary(W, y, score, recall_levels=(0.8, 0.9), fpr_level=0.01) -> dict:
    """ROC-AUC, PR-AUC, precision at each recall level, recall at fpr_level; one curve."""
    tp, fp, p, n, _ = curve(W, y, score)
    out = {"roc_auc": _roc(tp, fp, p, n), "pr_auc": _ap(tp, fp, p)}
    for r in recall_levels:
        out[f"precision_at_recall_{r}"] = _prec_at_recall(tp, fp, p, r)
    out[f"recall_at_fpr_{fpr_level}"] = _recall_at_fpr(tp, fp, p, n, fpr_level)
    return out


def fpr_at_recall(W, y, score, level, neg_mask) -> np.ndarray:
    """FPR inside the rows of `neg_mask` at the first threshold where recall >= level."""
    neg_mask = _bool(neg_mask) & ~_bool(y)
    tp, _, p, _, sub = curve(W, y, score, extra=neg_mask)
    W = as_weights(W, len(neg_mask))
    total = W @ neg_mask
    level = np.broadcast_to(np.asarray(level, dtype=np.float64), p.shape)[:, None]
    reached = _div(tp, p[:, None]) >= level - TOL
    out = np.full(len(p), np.nan)
    has = reached.any(axis=1) & (p > 0) & (total > 0)
    first = reached.argmax(axis=1)
    rows = np.nonzero(has)[0]
    out[rows] = sub[rows, first[rows]] / total[rows]
    return out


def brier(W, y, prob) -> np.ndarray:
    y = _bool(y).astype(np.float64)
    prob = np.asarray(prob, dtype=np.float64)
    if prob.shape != y.shape:
        raise ValueError(f"prob shape {prob.shape} differs from label shape {y.shape}")
    W = as_weights(W, len(y))
    return _div(W @ (prob - y) ** 2, W.sum(axis=1))


def precision_at_prevalence(rec, rate, prevalence: float) -> np.ndarray:
    """Precision at an assumed prevalence (spec 4, ruling C-11):
    recall * pi / (recall * pi + FPR * (1 - pi))."""
    rec = np.asarray(rec, dtype=np.float64)
    rate = np.asarray(rate, dtype=np.float64)
    return _div(rec * prevalence, rec * prevalence + rate * (1.0 - prevalence))
