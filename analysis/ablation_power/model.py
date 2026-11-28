"""The leave-one-clade-out fit used by the power study.

The BAR is the best of a declared family -- length-only, composition-only, composition+length --
chosen by inner cluster-grouped cross-validation inside the training part, never by looking at the
test fold. It is not "composition + length" by default: measured on adhesion_v1 that 21-feature
model scores a pooled AP of 0.551 while log length ALONE scores 0.690, so taking it as the bar
would let a component win by stabilising a noise-limited fit rather than by knowing anything about
adhesion.

Everything that learns -- the scaler and the regression -- is refit inside every fold, and, when
the bootstrap refits, inside every resample.
"""

import numpy as np
from benchmark_data import FEATURE_SETS, splits

C_GRID = (0.01, 0.1, 1.0, 10.0)
INNER_SPLITS = 5


class ModelError(Exception):
    pass


def _pipeline(c):
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    return Pipeline(
        [
            ("scale", StandardScaler()),
            ("lr", LogisticRegression(C=c, max_iter=5000, solver="lbfgs")),
        ]
    )


def _design(bench, feature_set, component=None):
    cols = list(FEATURE_SETS[feature_set])
    X = bench.features[:, cols]
    if component is not None:
        component = np.asarray(component, dtype=float)
        if component.ndim == 1:
            component = component[:, None]
        if component.shape[0] != X.shape[0]:
            raise ModelError(f"component has {component.shape[0]} rows, benchmark has {X.shape[0]}")
        X = np.hstack([X, component])
    return X


def _inner_choose_c(X, y, groups, seed):
    """Pick C by cluster-grouped CV inside the training part only."""
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold

    uniq = np.unique(groups)
    n_splits = min(INNER_SPLITS, uniq.size)
    if n_splits < 2 or len(np.unique(y)) < 2:
        return 1.0
    best_c, best_score = 1.0, -np.inf
    for c in C_GRID:
        scores = []
        for tr, te in GroupKFold(n_splits=n_splits).split(X, y, groups):
            if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                continue
            model = _pipeline(c).fit(X[tr], y[tr])
            scores.append(roc_auc_score(y[te], model.predict_proba(X[te])[:, 1]))
        if scores and np.mean(scores) > best_score:
            best_c, best_score = c, float(np.mean(scores))
    return best_c


def loco_scores(bench, feature_set, component=None, c=1.0, tune=False, seed=0):
    """Out-of-fold scores for every protein in a scorable fold, in benchmark order.

    Proteins outside the scorable folds get NaN: they can be trained on but are never scored.

    `c` may be a scalar or a {fold: C} mapping. The mapping matters: the bar and the
    component-augmented model MUST be fitted at the same regularisation, or their difference
    measures the change in C rather than the component. Fitting the bar with a tuned C and the
    component with a fixed one made a pure-noise component appear to cost 0.15 concordance.
    """
    X = _design(bench, feature_set, component)
    y = bench.y
    cluster = np.array(bench.cluster)
    out = np.full(len(bench), np.nan, dtype=float)
    chosen = {}
    for held, train, test, _dropped in splits(bench):
        if train.size == 0 or len(np.unique(y[train])) < 2:
            continue
        if tune:
            use_c = _inner_choose_c(X[train], y[train], cluster[train], seed)
        elif isinstance(c, dict):
            use_c = c.get(held, 1.0)
        else:
            use_c = c
        chosen[held] = use_c
        model = _pipeline(use_c).fit(X[train], y[train])
        out[test] = model.predict_proba(X[test])[:, 1]
    return out, chosen


def choose_bar(bench, tune=True, seed=0):
    """Evaluate the declared baseline family and return (name, scores, table, c_per_fold).

    The winner is the bar every component must beat. It is chosen once, on the real benchmark, by
    inner cluster-grouped CV inside each training part -- never by looking at a test fold.

    The per-fold C it settled on is returned so that every later fit, bar and component alike, can
    be pinned to it. Holding C fixed is what makes a difference in score attributable to the
    component rather than to the regularisation.
    """
    from stats import stratified_concordance

    scored = np.array([f in set(bench.scored_folds) for f in bench.fold])
    table = {}
    best_name, best_value, best_scores, best_c = None, -np.inf, None, {}
    for name in FEATURE_SETS:
        scores, chosen = loco_scores(bench, name, tune=tune, seed=seed)
        value = stratified_concordance(
            bench.y[scored], scores[scored], np.array(bench.fold)[scored]
        )
        table[name] = {"stratified_concordance": float(value), "C_per_fold": dict(chosen)}
        if value > best_value:
            best_name, best_value, best_scores, best_c = name, value, scores, dict(chosen)
    return best_name, best_scores, table, best_c
