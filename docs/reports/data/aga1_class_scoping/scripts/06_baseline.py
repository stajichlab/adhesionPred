#!/usr/bin/env python3
"""Descriptive baseline: how well do local sequence features alone separate class positives from hard negatives?

Needs numpy + scikit-learn (repo environment python). No tuning: logistic regression, L2 C=1, standardised features,
balanced class weights, fixed seeds. Training = positives + hard negatives of the training folds. Background proteins
(assumed negatives) are NEVER trained on; they are scored as extra test negatives.

Schemes
  leave-clade-out (LCO): hold out one fine clade (the clades with positives); train on the others.
  leave-cluster-out (LCL): 5 folds of whole MMseqs2 clusters (30 percent id, 50 percent coverage), repeated with 20 seeds.
Positive sets: tierA (independent repeat evidence) and all (tier A + B family-inference).
Intervals: percentile cluster bootstrap (resample whole clusters of the evaluated proteins, 2000 draws, seed 1).
An interval is NOT given when the evaluated set has fewer than 5 positives, fewer than 3 positive clusters,
fewer than 10 negatives or fewer than 3 negative clusters.

Outputs: baseline_lco.tsv, baseline_lcl.tsv, single_feature_auc.tsv, repeat_call_vs_features.tsv, oof_scores.tsv
"""

import math
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module  # noqa: E402

from common import OUT, write_tsv  # noqa: E402

inv = import_module("05_inventory")

FEATSETS = {
    "comp": ["frac_S", "frac_T", "frac_N", "frac_Q", "frac_G", "frac_P", "cys_frac"],
    "comp_len_struct": ["frac_S", "frac_T", "frac_N", "frac_Q", "frac_G", "frac_P", "cys_frac", "log_len", "max_ST30", "lc_frac", "mean_entropy30", "nterm_frac", "log_run", "sequon_per100"],
    "comp_len_struct_repeat": ["frac_S", "frac_T", "frac_N", "frac_Q", "frac_G", "frac_P", "cys_frac", "log_len", "max_ST30", "lc_frac", "mean_entropy30", "nterm_frac", "log_run", "sequon_per100", "r14_z_seq", "r14_coverage", "r02_coverage", "r14_copies", "r02_copies"],
}  # fmt: skip
NEG_SETS = ["hn_all", "hn_st30", "hn_stfrac25", "hn_repeat", "bg"]


def auc(pos, neg):
    """Mann-Whitney AUC with ties counted as 0.5."""
    if not len(pos) or not len(neg):
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv))
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    r = ranks[: len(pos)].sum()
    return (r - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def cluster_boot(score, is_pos, cluster, n=2000, seed=1):
    """Percentile CI of the AUC by resampling clusters; None when too few positives."""
    score, is_pos, cluster = map(np.asarray, (score, is_pos, cluster))
    npos = int(is_pos.sum())
    pos_cl = len(set(cluster[is_pos]))
    nneg = int((~is_pos).sum())
    neg_cl = len(set(cluster[~is_pos]))
    if npos < 5 or pos_cl < 3 or nneg < 10 or neg_cl < 3:
        return None
    groups = defaultdict(list)
    for i, c in enumerate(cluster):
        groups[c].append(i)
    keys = list(groups)
    rng = random.Random(seed)
    vals = []
    for _ in range(n):
        idx = [i for k in (rng.choice(keys) for _ in keys) for i in groups[k]]
        p = score[idx][is_pos[idx]]
        q = score[idx][~is_pos[idx]]
        if len(p) and len(q):
            vals.append(auc(p, q))
    vals.sort()
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]


def wilson(k, n):
    if n == 0:
        return float("nan"), float("nan")
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def matrix(P, accs, fs):
    rows = []
    for a in accs:
        p = P[a]
        v = {k: float(p[k]) if p.get(k, "") != "" else 0.0 for k in ("frac_S", "frac_T", "frac_N", "frac_Q", "frac_G", "frac_P", "max_ST30", "lc_frac", "mean_entropy30", "nterm_frac", "sequon_per100", "r14_z_seq", "r14_coverage", "r02_coverage", "r14_copies", "r02_copies")}  # fmt: skip
        v["cys_frac"] = float(p["n_cys"]) / float(p["length"])
        v["log_len"] = math.log10(float(p["length"]))
        v["log_run"] = math.log10(float(p["longest_ST_run"]) + 1)
        rows.append([v[k] for k in fs])
    return np.array(rows)


def fit_score(P, train, test, fs):
    """train: [(acc, y)], test: [acc]. Returns scores of test."""
    Xtr = matrix(P, [a for a, _ in train], fs)
    ytr = np.array([y for _, y in train])
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit(
        sc.transform(Xtr), ytr
    )
    return m.decision_function(sc.transform(matrix(P, test, fs)))


def main():
    P = inv.load()
    S = {a: set(inv.subsets(p)) for a, p in P.items()}
    posset = {"tierA": "pos_tierA", "all": "pos_all"}
    hn = [a for a in P if "hn_all" in S[a]]
    bg = [a for a in P if "bg" in S[a]]
    oof_rows = []

    def evaluate(rows, tag, scheme, fold, fs, ps, test_pos, score):
        for ns in NEG_SETS:
            negs = [a for a in score if ns in S[a]]
            pos = list(test_pos)
            if not pos or not negs:
                continue
            a_ = auc(np.array([score[a] for a in pos]), np.array([score[a] for a in negs]))
            ids = pos + negs
            ci = cluster_boot(
                [score[a] for a in ids],
                [a in set(pos) for a in ids],
                [P[a]["cluster"] for a in ids],
            )
            rows.append({"scheme": scheme, "fold": fold, "featset": fs, "pos_set": ps, "neg_set": ns, "n_pos": len(pos), "n_pos_clusters": len({P[a]["cluster"] for a in pos}), "n_neg": len(negs), "auc": round(a_, 3), "ci_lo": "" if ci is None else round(ci[0], 3), "ci_hi": "" if ci is None else round(ci[1], 3), "interval": "cluster bootstrap" if ci else "n too small"})  # fmt: skip

    # ---- leave-clade-out
    lco = []
    clades = sorted({P[a]["clade_fine"] for a in P if "pos_all" in S[a]})
    for fs, cols in FEATSETS.items():
        for ps, sub in posset.items():
            pos_all = [a for a in P if sub in S[a]]
            for clade in clades:
                test_pos = [a for a in pos_all if P[a]["clade_fine"] == clade]
                if not test_pos:
                    continue
                train = [(a, 1) for a in pos_all if P[a]["clade_fine"] != clade] + [
                    (a, 0) for a in hn if P[a]["clade_fine"] != clade
                ]
                test = test_pos + [a for a in hn + bg if P[a]["clade_fine"] == clade]
                if sum(y for _, y in train) < 2:
                    continue
                sc = dict(zip(test, fit_score(P, train, test, cols)))  # noqa: B905
                evaluate(lco, "", "leave-clade-out", clade, fs, ps, test_pos, sc)
                if fs == "comp_len_struct" and ps == "all":
                    for a in test:
                        oof_rows.append(
                            {"accession": a, "scheme": "LCO", "score": round(float(sc[a]), 4)}
                        )
    write_tsv(
        OUT / "baseline_lco.tsv",
        lco,
        [
            "scheme",
            "fold",
            "featset",
            "pos_set",
            "neg_set",
            "n_pos",
            "n_pos_clusters",
            "n_neg",
            "auc",
            "ci_lo",
            "ci_hi",
            "interval",
        ],
    )

    # ---- leave-cluster-out, 5 folds x 20 repeats
    lcl = []
    for fs, cols in FEATSETS.items():
        for ps, sub in posset.items():
            pos_all = [a for a in P if sub in S[a]]
            per_rep = defaultdict(list)
            for seed in range(20):
                rng = random.Random(seed)
                pcl = sorted({P[a]["cluster"] for a in pos_all})
                ocl = sorted({P[a]["cluster"] for a in hn + bg} - set(pcl))
                rng.shuffle(pcl)
                rng.shuffle(ocl)
                fold_of = {c: i % 5 for i, c in enumerate(pcl)}
                fold_of.update({c: i % 5 for i, c in enumerate(ocl)})
                score = {}
                for k in range(5):
                    train = [(a, 1) for a in pos_all if fold_of[P[a]["cluster"]] != k] + [
                        (a, 0) for a in hn if fold_of[P[a]["cluster"]] != k
                    ]
                    test = [a for a in pos_all + hn + bg if fold_of[P[a]["cluster"]] == k]
                    score.update(dict(zip(test, fit_score(P, train, test, cols))))  # noqa: B905
                for ns in NEG_SETS:
                    negs = [a for a in score if ns in S[a]]
                    per_rep[ns].append(
                        auc(
                            np.array([score[a] for a in pos_all]),
                            np.array([score[a] for a in negs]),
                        )
                    )
                if seed == 0:
                    evaluate(
                        lcl,
                        "",
                        "leave-cluster-out (seed 0 pooled OOF)",
                        "all folds",
                        fs,
                        ps,
                        pos_all,
                        score,
                    )
                    if fs == "comp_len_struct" and ps == "all":
                        for a in score:
                            oof_rows.append(
                                {
                                    "accession": a,
                                    "scheme": "LCL_seed0",
                                    "score": round(float(score[a]), 4),
                                }
                            )
            for ns in NEG_SETS:
                v = per_rep[ns]
                lcl.append({"scheme": "leave-cluster-out (20 seeds)", "fold": "all folds", "featset": fs, "pos_set": ps, "neg_set": ns, "n_pos": len(pos_all), "n_pos_clusters": len({P[a]["cluster"] for a in pos_all}), "n_neg": len([a for a in P if ns in S[a]]), "auc": round(float(np.mean(v)), 3), "ci_lo": round(min(v), 3), "ci_hi": round(max(v), 3), "interval": "min-max over 20 seeds (split variability only)"})  # fmt: skip
    write_tsv(
        OUT / "baseline_lcl.tsv",
        lcl,
        [
            "scheme",
            "fold",
            "featset",
            "pos_set",
            "neg_set",
            "n_pos",
            "n_pos_clusters",
            "n_neg",
            "auc",
            "ci_lo",
            "ci_hi",
            "interval",
        ],
    )
    write_tsv(OUT / "oof_scores.tsv", oof_rows, ["accession", "scheme", "score"])

    # ---- single features (no training): AUC of the raw value, higher = more positive-like
    single = []
    feats = ["length", "frac_S", "frac_T", "frac_ST", "frac_STP", "frac_STNQGP", "frac_N", "frac_Q", "frac_G", "frac_P", "n_cys", "max_ST30", "lc_frac", "mean_entropy30", "nterm_frac", "longest_ST_run", "sequon_per100", "r14_z_seq", "r14_coverage", "r02_coverage"]  # fmt: skip
    for ps, sub in posset.items():
        pos = [a for a in P if sub in S[a]]
        for f in feats:
            for ns in NEG_SETS:
                negs = [a for a in P if ns in S[a]]
                v = lambda a, f=f: float(P[a][f]) if P[a][f] != "" else 0.0  # noqa: E731
                ids = pos + negs
                ci = cluster_boot(
                    [v(a) for a in ids],
                    [a in set(pos) for a in ids],
                    [P[a]["cluster"] for a in ids],
                )
                single.append({"feature": f, "pos_set": ps, "neg_set": ns, "n_pos": len(pos), "n_neg": len(negs), "auc_raw": round(auc(np.array([v(a) for a in pos]), np.array([v(a) for a in negs])), 3), "ci_lo": "" if ci is None else round(ci[0], 3), "ci_hi": "" if ci is None else round(ci[1], 3)})  # fmt: skip
    write_tsv(
        OUT / "single_feature_auc.tsv",
        single,
        ["feature", "pos_set", "neg_set", "n_pos", "n_neg", "auc_raw", "ci_lo", "ci_hi"],
    )

    # ---- existing repeat call on the same proteins
    rc = []
    for ps, sub in posset.items():
        pos = [a for a in P if sub in S[a]]
        for rule in ("call_repeat", "call_repeat_z35"):
            k = sum(int(P[a][rule]) for a in pos)
            lo, hi = wilson(k, len(pos))
            rc.append({"rule": rule, "set": ps, "kind": "positives called (sensitivity)", "n": len(pos), "called": k, "rate": round(k / len(pos), 3), "wilson_lo": round(lo, 3), "wilson_hi": round(hi, 3)})  # fmt: skip
            for ns in NEG_SETS:
                negs = [a for a in P if ns in S[a]]
                k = sum(int(P[a][rule]) for a in negs)
                lo, hi = wilson(k, len(negs))
                rc.append({"rule": rule, "set": ns, "kind": "negatives called (false-positive rate)", "n": len(negs), "called": k, "rate": round(k / len(negs), 3), "wilson_lo": round(lo, 3), "wilson_hi": round(hi, 3)})  # fmt: skip
    write_tsv(
        OUT / "repeat_call_vs_features.tsv",
        rc,
        ["rule", "set", "kind", "n", "called", "rate", "wilson_lo", "wilson_hi"],
    )
    print("done")


if __name__ == "__main__":
    main()
