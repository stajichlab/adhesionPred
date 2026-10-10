#!/usr/bin/python3.12
"""Statistics of the secretome discovery pass (B4 of docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md, section 5.2).

Input per family: one row per genus: genus, group ("target" or "background"), pool (pool size), hits (family hits in the pool).
Presence of a genus = hits > 0. Copies per 1,000 = 1000 * hits / pool.
Fisher p and q are a ranking statistic under an independence assumption that does not hold (shared ancestry); see the spec.
"""

import io

import numpy as np
from scipy import stats as sst
from scipy.optimize import minimize


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if p.size == 0:
        return []
    return [float(x) for x in sst.false_discovery_control(p, method="bh")]


def _logit_p(rows):
    """Likelihood-ratio p for the group term in a logistic regression of presence on group and log pool size."""
    y = np.array([1.0 if r["hits"] > 0 else 0.0 for r in rows])
    g = np.array([1.0 if r["group"] == "target" else 0.0 for r in rows])
    z = np.log(np.array([r["pool"] for r in rows], dtype=float))
    z = z - z.mean()
    if y.sum() in (0, len(y)):
        return float("nan")

    def nll(beta, X):
        eta = X @ beta
        return float(np.sum(np.logaddexp(0, eta) - y * eta) + 1e-6 * np.sum(beta[1:] ** 2))

    X1 = np.column_stack([np.ones_like(y), z])
    X2 = np.column_stack([np.ones_like(y), z, g])
    f1 = minimize(nll, np.zeros(2), args=(X1,), method="BFGS").fun
    f2 = minimize(nll, np.zeros(3), args=(X2,), method="BFGS").fun
    return float(sst.chi2.sf(max(0.0, 2 * (f1 - f2)), 1))


def family_stats(rows):
    t = [r for r in rows if r["group"] == "target"]
    b = [r for r in rows if r["group"] == "background"]
    tp = sum(1 for r in t if r["hits"] > 0)
    bp = sum(1 for r in b if r["hits"] > 0)
    table = [[tp, len(t) - tp], [bp, len(b) - bp]]
    odds, p = sst.fisher_exact(table, alternative="two-sided")
    tc = [1000 * r["hits"] / r["pool"] for r in t]
    bc = [1000 * r["hits"] / r["pool"] for r in b]
    try:
        mw = float(sst.mannwhitneyu(tc, bc, alternative="two-sided", method="asymptotic").pvalue)
    except ValueError:
        mw = float("nan")
    return {
        "n_target": len(t),
        "n_background": len(b),
        "t_present": tp,
        "b_present": bp,
        "t_frac": tp / len(t) if t else float("nan"),
        "b_frac": bp / len(b) if b else float("nan"),
        "odds_ratio": float(odds),
        "fisher_p": float(p),
        "logit_p": _logit_p(rows),
        "t_copies_median": float(np.median(tc)) if tc else float("nan"),
        "b_copies_median": float(np.median(bc)) if bc else float("nan"),
        "mwu_p": mw,
    }


def report(families, min_genera=3):
    """families: {name: rows}. BH runs over all families with at least one present genus; the reporting rule limits what is listed."""
    stats = {}
    for name, rows in families.items():
        if not any(r["hits"] > 0 for r in rows):
            continue
        stats[name] = family_stats(rows)
    names = sorted(stats)
    q = bh([stats[n]["fisher_p"] for n in names])
    enriched, depleted = [], []
    for n, qi in zip(names, q, strict=True):
        s = dict(stats[n], family=n, fisher_q=qi)
        if s["t_frac"] > s["b_frac"] and s["t_present"] >= min_genera:
            enriched.append(s)
        elif s["t_frac"] < s["b_frac"]:
            depleted.append(s)
    key = lambda s: s["fisher_p"]  # noqa: E731
    return {
        "enriched": sorted(enriched, key=key),
        "depleted": sorted(depleted, key=key),
        "n_tested": len(names),
    }


def fitch_gains_losses(newick, present, all_tips):
    """Gains and losses of presence on a rooted tree: Fitch two-pass reconstruction with ties resolved to the parent state
    (a root tie resolves to absent), then 0->1 and 1->0 changes counted on the branches. Unrooted input is the caller's problem."""
    from Bio import Phylo

    tree = Phylo.read(io.StringIO(newick), "newick")
    tips = set(all_tips)

    def down(node):
        if node.is_terminal():
            node.sets = {1} if node.name in present else {0}
            return
        for c in node.clades:
            down(c)
        inter = set.intersection(*[c.sets for c in node.clades])
        node.sets = inter if inter else set.union(*[c.sets for c in node.clades])

    def up(node, parent_state):
        if parent_state is None:
            node.state = 0 if 0 in node.sets else 1
        else:
            node.state = parent_state if parent_state in node.sets else sorted(node.sets)[0]
        for c in node.clades:
            up(c, node.state)

    down(tree.root)
    up(tree.root, None)
    gains = losses = 0
    stack = [tree.root]
    while stack:
        n = stack.pop()
        for c in n.clades:
            if n.state == 0 and c.state == 1:
                gains += 1
            elif n.state == 1 and c.state == 0:
                losses += 1
            stack.append(c)
    assert tips >= {n.name for n in tree.get_terminals()}
    return gains, losses


__all__ = ["bh", "family_stats", "report", "fitch_gains_losses"]
