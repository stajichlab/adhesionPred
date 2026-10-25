"""Features-only classification study for the AGA1 class (scikit-learn, numpy, pandas, matplotlib; no sequence is written).

Sub-commands
  features  build the per-protein feature table (numbers only, keyed by accession; no residues, no repeat-unit text)
  analyse   PCA, grouped cross-validation of three models, hand rules, tier B and excluded-role scores, clade hold-out,
            ablation, permutation importance, length checks; writes tables and a PNG into --out-dir
  apply     fit one pipeline on ALL tier A plus the negatives and score a NEW feature table

Labels (owner decisions 2026-10-24/25)
  positives  = role class_positive and pos_tier A (never B)
  negatives  = role hard_negative + background_assumed_negative
  extras     = tier B, excluded_*, adhesin_other_mechanism, kiss_row_positive, uncurated_repeat_bearing: NEVER fitted on, never in
               headline metrics, scored descriptively by models that did not see their cluster
  SED1 and DAN4 are negatives in the headline; a set-aside run leaves them out entirely.

Every imputation and scaling step lives inside the scikit-learn Pipeline, so it is refit inside each fold. Folds are grouped by
MMseqs cluster (a cluster-less protein is its own cluster), so no cluster spans train and test. Thresholds for fixed false-positive
rates are read from an INNER grouped cross-validation on the training part only.

Fails loudly (FATAL:, exit 1) on missing columns or an existing output file.
"""

import argparse
import hashlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

ST_SET = frozenset("ST")
AA20 = "ACDEFGHIKLMNPQRSTVWY"
CHARGED, HYDROPHOBIC = "DEKR", "AVILMFW"
GPI_ORD = {"none": 0.0, "weak": 1.0, "probable": 2.0, "high": 3.0}
SEQUON = re.compile(r"(?=N[^P][ST])")
NEW_COL_SED = ("Q01589", "P47179")  # SED1, DAN4

GROUPS = {
    "composition": [
        "frac_S", "frac_T", "frac_ST", "frac_P", "frac_N", "frac_Q", "frac_G", "frac_A", "frac_charged", "frac_C",
        "frac_hydrophobic",
    ],
    "length": ["length"],
    "predictors": ["sp_yes", "sp_prob", "gpi_yes", "gpi_missing", "gpi_conf_ord"],
    "disorder": ["disorder_frac", "disorder_longest_run", "disorder_mean", "st_in_disordered"],
    "repeats": ["rep_coverage", "rep_n_copies", "rep_z_seq", "rep_missing"],
    "stalk": ["stalk_run_len", "max_st_window_frac", "ctail_hydro"],
    "glyco": ["sequon_per100"],
    "complexity": ["aa_entropy"],
}  # fmt: skip
FEATURES = [f for g in GROUPS.values() for f in g]
LOG_COLS = {"length", "disorder_longest_run", "stalk_run_len", "rep_n_copies"}
META = ["accession", "gene", "role", "pos_tier", "clade_coarse", "cluster_id"]
NEG_ROLES = ("hard_negative", "background_assumed_negative")
FPRS = (0.01, 0.035, 0.05)
MODELS = ("logreg", "rf", "gboost")
N_SPLITS, N_INNER = 5, 3


class FeatureError(Exception):
    pass


def group_of(feature):
    for g, fs in GROUPS.items():
        if feature in fs:
            return g
    raise KeyError(feature)


# ---------------------------------------------------------------- features
def sequon_per100(residues):
    """N-x-S/T (x not Pro) sites per 100 residues; overlapping sites each count; NA for an empty sequence"""
    if not residues:
        return float("nan")
    return 100.0 * len(SEQUON.findall(residues)) / len(residues)


def composition(residues):
    n = len(residues)
    if n == 0:
        raise FeatureError("empty sequence")
    f = lambda letters: sum(residues.count(a) for a in letters) / n  # noqa: E731
    out = {f"frac_{a}": f(a) for a in "STPNQGAC"}
    out["frac_ST"] = out["frac_S"] + out["frac_T"]
    out["frac_charged"] = f(CHARGED)
    out["frac_hydrophobic"] = f(HYDROPHOBIC)
    return out


def aa_entropy(residues):
    """Shannon entropy (bits) of the 20-letter composition; low values mean low complexity"""
    counts = np.array([residues.count(a) for a in AA20], float)
    tot = counts.sum()
    if tot == 0:
        return float("nan")
    p = counts[counts > 0] / tot
    return float(-(p * np.log2(p)).sum())


def st_in_disordered(residues, scores, threshold=0.5):
    """Ser/Thr fraction among residues with disorder score >= threshold (NaN if none are disordered)"""
    if len(residues) != len(scores):
        raise FeatureError("residues and scores differ in length")
    dis = [a for a, s in zip(residues, scores, strict=False) if s >= threshold]
    return sum(a in ST_SET for a in dis) / len(dis) if dis else float("nan")


def sequence_features(residues, scores):
    out = composition(residues)
    out["sequon_per100"] = sequon_per100(residues)
    out["aa_entropy"] = aa_entropy(residues)
    out["st_in_disordered"] = st_in_disordered(residues, scores)
    return out


def gpi_encoding(gpi_pred, conf):
    """(gpi_yes, gpi_missing, gpi_conf_ord). 'NA' (protein <= 40 aa) stays missing with an indicator, never 0."""
    if gpi_pred == "NA":
        return float("nan"), 1.0, float("nan")
    if gpi_pred not in ("yes", "no"):
        raise FeatureError(f"unexpected gpi_pred {gpi_pred!r}")
    if gpi_pred == "no":
        return 0.0, 0.0, 0.0
    if conf not in GPI_ORD:
        raise FeatureError(f"unexpected gpi_confidence_class {conf!r} for a GPI call")
    return 1.0, 0.0, GPI_ORD[conf]


def to_float(x):
    return float("nan") if x in ("NA", "", None) else float(x)


def build_features(proteins, preds, stalk, rep, raw):
    """proteins/preds/stalk/rep: DataFrames; raw: accession -> (residues, scores). Returns the feature table."""
    acc = list(proteins["accession"])
    for name, df, key in (("predictions", preds, "accession"), ("stalk", stalk, "accession")):
        if set(df[key]) != set(acc) or df[key].duplicated().any():
            raise FeatureError(f"{name} table accessions differ from proteins.tsv")
    if set(raw) != set(acc):
        raise FeatureError("raw AIUPred accessions differ from proteins.tsv")
    p = preds.set_index("accession")
    s = stalk.set_index("accession")
    r = rep.set_index("protein")
    rows = []
    for a in acc:
        res, sc = raw[a]
        row = {"accession": a, "length": float(len(res))}
        if int(p.at[a, "disorder_n_res"]) != len(res) or int(s.at[a, "n_res"]) != len(res):
            raise FeatureError(f"{a}: residue count differs between tables")
        row.update(sequence_features("".join(res), sc))
        row["sp_yes"] = 1.0 if p.at[a, "signal_peptide_pred"] == "yes" else 0.0
        row["sp_prob"] = to_float(p.at[a, "sp_prob"])
        row["gpi_yes"], row["gpi_missing"], row["gpi_conf_ord"] = gpi_encoding(
            p.at[a, "gpi_pred"], p.at[a, "gpi_confidence_class"]
        )
        for c in ("disorder_frac", "disorder_longest_run", "disorder_mean"):
            row[c] = to_float(p.at[a, c])
        if a in r.index:
            row["rep_coverage"] = float(r.at[a, "rep_coverage"])
            row["rep_n_copies"] = float(r.at[a, "rep_n_copies"])
            row["rep_z_seq"] = float(r.at[a, "rep_z_seq"])
            row["rep_missing"] = 0.0
        else:  # detector skips proteins < 80 aa: missing, never 0
            row.update(rep_coverage=np.nan, rep_n_copies=np.nan, rep_z_seq=np.nan, rep_missing=1.0)
        for c in ("stalk_run_len", "max_st_window_frac", "ctail_hydro"):
            row[c] = to_float(s.at[a, c])
        rows.append(row)
    feat = pd.DataFrame(rows)
    meta = proteins[["accession", "gene", "role", "pos_tier", "clade_coarse", "cluster_id"]].copy()
    meta["cluster_id"] = [
        c if c else f"SINGLETON_{a}"
        for a, c in zip(meta["accession"], meta["cluster_id"], strict=False)
    ]
    return meta.merge(feat, on="accession")[META + FEATURES]


def rule_flags(df):
    """hand rules on a feature table: broad = SP and GPI and longest disordered run >= 200; narrow = broad and project repeat call"""
    broad = (
        (df["sp_yes"] == 1) & (df["gpi_yes"] == 1) & (df["disorder_longest_run"] >= 200)
    ).to_numpy()
    rep = (
        (df["rep_coverage"] >= 0.25) & (df["rep_n_copies"] >= 2.5)
    ).to_numpy()  # NaN compares False
    return broad, broad & rep


# ---------------------------------------------------------------- labels
def label_sets(df, drop=()):
    """masks: eligible (tier A + negatives, minus drop), y over eligible, extras (everything else not dropped)"""
    pos = ((df["role"] == "class_positive") & (df["pos_tier"] == "A")).to_numpy()
    neg = df["role"].isin(NEG_ROLES).to_numpy()
    keep = ~df["accession"].isin(drop).to_numpy()
    elig = (pos | neg) & keep
    return elig, pos, neg, ~(pos | neg) & keep


# ---------------------------------------------------------------- models
def _log_cols(X, idx):
    X = np.array(X, float, copy=True)
    X[:, idx] = np.log1p(np.clip(X[:, idx], 0, None))
    return X


def build_pipeline(model, cols, seed=0):
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import FunctionTransformer, StandardScaler

    idx = [i for i, c in enumerate(cols) if c in LOG_COLS]
    if model == "logreg":
        clf = LogisticRegression(C=0.1, class_weight="balanced", max_iter=5000)
    elif model == "rf":
        clf = RandomForestClassifier(
            n_estimators=300,
            min_samples_leaf=2,
            class_weight="balanced_subsample",
            random_state=seed,
        )
    elif model == "gboost":
        clf = GradientBoostingClassifier(
            n_estimators=100, max_depth=2, learning_rate=0.1, subsample=0.8, random_state=seed
        )
    else:
        raise FeatureError(f"unknown model {model!r}")
    return Pipeline(
        [
            ("log", FunctionTransformer(_log_cols, kw_args={"idx": idx}, validate=False)),
            ("imp", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("sc", StandardScaler()),
            ("clf", clf),
        ]
    )


def fit_scores(model, cols, seed, Xtr, ytr, *Xtests):
    pipe = build_pipeline(model, cols, seed).fit(Xtr, ytr)
    return [pipe.predict_proba(X)[:, 1] for X in Xtests]


def grouped_splits(y, groups, n_splits, seed):
    from sklearn.model_selection import StratifiedGroupKFold

    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(cv.split(np.zeros(len(y)), y, groups))


def inner_thresholds(model, cols, seed, X, y, groups, fprs, n_inner=N_INNER):
    """score thresholds giving the target false-positive rates on INNER out-of-fold negatives of the training part"""
    oof = np.full(len(y), np.nan)
    for tr, te in grouped_splits(y, groups, n_inner, seed):
        oof[te] = fit_scores(model, cols, seed, X[tr], y[tr], X[te])[0]
    neg = oof[y == 0]
    return [float(np.quantile(neg, 1 - f, method="lower")) for f in fprs]


def cv_one_seed(model, cols, seed, X, y, g, XX, gx, fprs, perm_groups=None, n_splits=N_SPLITS):
    """one grouped CV repetition. Returns out-of-fold scores, calls at each fpr, and extras bookkeeping."""
    n, m = len(y), len(XX)
    oof = np.full(n, np.nan)
    called = np.zeros((len(fprs), n), bool)
    ex_s, ex_n = np.zeros(m), np.zeros(m)
    ex_c = np.zeros((len(fprs), m))
    rng = np.random.default_rng(seed)
    perm = {k: np.full(n, np.nan) for k in (perm_groups or {})}
    for tr, te in grouped_splits(y, g, n_splits, seed):
        if set(g[tr]) & set(g[te]):
            raise AssertionError("a cluster spans train and test")
        pipe = build_pipeline(model, cols, seed).fit(X[tr], y[tr])
        oof[te] = pipe.predict_proba(X[te])[:, 1]
        thr = inner_thresholds(model, cols, seed, X[tr], y[tr], g[tr], fprs) if fprs else []
        for k, t in enumerate(thr):
            called[k, te] = oof[te] > t
        if m:
            unseen = ~np.isin(gx, list(set(g[tr])))
            if unseen.any():
                s = pipe.predict_proba(XX[unseen])[:, 1]
                ex_s[unseen] += s
                ex_n[unseen] += 1
                for k, t in enumerate(thr):
                    ex_c[k, unseen] += s > t
        for name, idx in (perm_groups or {}).items():
            Xp = X[te].copy()
            Xp[:, idx] = Xp[rng.permutation(len(te))][:, idx]
            perm[name][te] = pipe.predict_proba(Xp)[:, 1]
    return {"oof": oof, "called": called, "ex_s": ex_s, "ex_n": ex_n, "ex_c": ex_c, "perm": perm}


def run_cv(model, cols, X, y, g, XX, gx, fprs=FPRS, seeds=None, perm_groups=None, n_jobs=1):
    from joblib import Parallel, delayed

    seeds = range(10) if seeds is None else seeds

    outs = Parallel(n_jobs=n_jobs)(
        delayed(cv_one_seed)(model, cols, s, X, y, g, XX, gx, fprs, perm_groups) for s in seeds
    )
    res = {
        "oof_seeds": np.stack([o["oof"] for o in outs]),
        "called": np.stack([o["called"] for o in outs]),  # seeds x fprs x n
    }
    res["oof"] = res["oof_seeds"].mean(0)
    res["call_rate"] = res["called"].mean(0)  # fprs x n
    ex_n = np.sum([o["ex_n"] for o in outs], 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        res["ex_score"] = np.sum([o["ex_s"] for o in outs], 0) / ex_n
        res["ex_call"] = np.sum([o["ex_c"] for o in outs], 0) / ex_n
    res["perm"] = {k: np.stack([o["perm"][k] for o in outs]).mean(0) for k in (perm_groups or {})}
    return res


# ---------------------------------------------------------------- statistics
def auc(score, y):
    y = np.asarray(y)
    pos = int((y == 1).sum())
    neg = len(y) - pos
    if pos == 0 or neg == 0:
        return float("nan")
    r = pd.Series(np.asarray(score, float)).rank(method="average").to_numpy()
    return float((r[y == 1].sum() - pos * (pos + 1) / 2) / (pos * neg))


def cluster_bootstrap(stat, groups, n_boot=2000, seed=0):
    """(point, lo, hi): percentile 95% interval of stat(index array) over cluster resamples"""
    groups = np.asarray(groups)
    uniq, inv = np.unique(groups, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(len(uniq))]
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(uniq), len(uniq))
        v = stat(np.concatenate([members[i] for i in pick]))
        if not np.isnan(v):
            vals.append(v)
    point = stat(np.arange(len(groups)))
    if not vals:
        return point, float("nan"), float("nan")
    return point, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def fmt(t, d=3):
    return f"{t[0]:.{d}f} [{t[1]:.{d}f}, {t[2]:.{d}f}]"


def recall_stat(rate, y):
    return lambda i: float(rate[i][y[i] == 1].mean()) if (y[i] == 1).any() else float("nan")


def fpr_stat(rate, y):
    return lambda i: float(rate[i][y[i] == 0].mean()) if (y[i] == 0).any() else float("nan")


def rule_table(name, flag, y, g, n_boot):
    f = flag.astype(float)
    rec = cluster_bootstrap(recall_stat(f, y), g, n_boot)
    fp = cluster_bootstrap(fpr_stat(f, y), g, n_boot)
    return {"method": name, "recall": rec, "fpr": fp}


# ---------------------------------------------------------------- PCA
def pca_table(df, elig):
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler

    X = _log_cols(
        df[FEATURES].to_numpy(float), [i for i, c in enumerate(FEATURES) if c in LOG_COLS]
    )
    imp = SimpleImputer(strategy="median", keep_empty_features=True).fit(X[elig])
    sc = StandardScaler().fit(imp.transform(X[elig]))
    Z = sc.transform(imp.transform(X))
    pca = PCA(n_components=10, random_state=0).fit(Z[elig])
    coords = pca.transform(Z)
    cd = pd.DataFrame(coords[:, :6], columns=[f"PC{i + 1}" for i in range(6)])
    cd.insert(0, "accession", df["accession"].to_numpy())
    cd.insert(1, "role", df["role"].to_numpy())
    cd.insert(2, "pos_tier", df["pos_tier"].to_numpy())
    cd.insert(3, "fitted_on", np.where(elig, "eligible", "projected_only"))
    var = pd.DataFrame(
        {
            "pc": [f"PC{i + 1}" for i in range(10)],
            "variance_explained": pca.explained_variance_ratio_,
            "cumulative": np.cumsum(pca.explained_variance_ratio_),
        }
    )
    load = pd.DataFrame(
        pca.components_[:6].T, index=FEATURES, columns=[f"PC{i + 1}" for i in range(6)]
    )
    load.insert(0, "group", [group_of(f) for f in FEATURES])
    load.index.name = "feature"
    return cd, var, load.reset_index()


def group_loadings(load):
    pcs = [c for c in load.columns if c.startswith("PC")]
    return (load.groupby("group")[pcs].agg(lambda s: float((s**2).sum()))).reset_index()


def pca_plot(cd, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {
        "hard_negative": "#9aa5b1", "background_assumed_negative": "#cfd6dd", "class_positive_A": "#d62728",
        "class_positive_B": "#ff9896", "adhesin_other_mechanism": "#2ca02c",
        "excluded_indirect_regulator": "#9467bd", "excluded_surface_other_adhesion_phenotype": "#8c564b",
        "kiss_row_positive": "#17becf", "uncurated_repeat_bearing": "#bcbd22",
    }  # fmt: skip
    lab = np.where(
        cd["role"] == "class_positive", "class_positive_" + cd["pos_tier"].astype(str), cd["role"]
    )
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (a, b) in zip(axes, (("PC1", "PC2"), ("PC1", "PC3")), strict=False):
        for k in colors:
            sel = lab == k
            if sel.any():
                ax.scatter(
                    cd.loc[sel, a], cd.loc[sel, b], s=34 if "positive" in k else 14, c=colors[k], label=f"{k} ({sel.sum()})",
                    edgecolors="k" if k == "class_positive_A" else "none", linewidths=0.6, alpha=0.9,
                )  # fmt: skip
        ax.set_xlabel(a)
        ax.set_ylabel(b)
    axes[0].legend(fontsize=6.5, loc="best")
    fig.suptitle(
        "PCA of the AGA1-class feature table (fitted on tier A + negatives; other roles projected)"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- analysis
def _arrays(df, cols, mask):
    return df.loc[mask, cols].to_numpy(float)


def analyse(df, out_dir, seeds=None, n_boot=2000, n_jobs=-1, abl_seeds=None):
    seeds = range(10) if seeds is None else seeds
    abl_seeds = range(3) if abl_seeds is None else abl_seeds
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    targets = [
        "pca_coords.tsv", "pca_variance.tsv", "pca_loadings.tsv", "pca_group_loadings.tsv", "pca_roles.png",
        "headline_metrics.tsv", "headline_setaside_sed1_dan4.tsv", "tierB_recall.tsv", "paired_vs_rules.tsv", "extra_roles_scores.tsv", "clade_holdout.tsv",
        "matched_fpr_descriptive.tsv", "ablation.tsv", "permutation_importance.tsv", "length_check.tsv", "single_feature_auc.tsv", "oof_scores.tsv",
    ]  # fmt: skip
    for t in targets:
        if (out / t).exists():
            raise FeatureError(f"refusing to overwrite {out / t}")
    elig, pos, neg, extra = label_sets(df)
    y_all = pos.astype(int)

    # PCA
    cd, var, load = pca_table(df, elig)
    cd.to_csv(out / "pca_coords.tsv", sep="\t", index=False, float_format="%.5f")
    var.to_csv(out / "pca_variance.tsv", sep="\t", index=False, float_format="%.5f")
    load.to_csv(out / "pca_loadings.tsv", sep="\t", index=False, float_format="%.4f")
    group_loadings(load).to_csv(
        out / "pca_group_loadings.tsv", sep="\t", index=False, float_format="%.4f"
    )
    pca_plot(cd, out / "pca_roles.png")

    broad, narrow = rule_flags(df)
    g_e = df.loc[elig, "cluster_id"].to_numpy()
    y_e = y_all[elig]
    X_e = _arrays(df, FEATURES, elig)
    X_x = _arrays(df, FEATURES, extra)
    g_x = df.loc[extra, "cluster_id"].to_numpy()

    results = {}
    for mdl in MODELS:
        print(f"[cv] {mdl}", flush=True)
        results[mdl] = run_cv(mdl, FEATURES, X_e, y_e, g_e, X_x, g_x, FPRS, seeds, n_jobs=n_jobs)

    def headline(res_by_model, y, g, flags_e, label_extra, tierB_mask_x, boot):
        rows = []
        for mdl, r in res_by_model.items():
            row = {
                "method": mdl,
                "auc": cluster_bootstrap(lambda i, r=r: auc(r["oof"][i], y[i]), g, boot),
            }
            sd = [auc(s, y) for s in r["oof_seeds"]]
            row["auc_seed_range"] = f"{min(sd):.3f}-{max(sd):.3f}"
            for k, f in enumerate(FPRS):
                row[f"recall@{f}"] = cluster_bootstrap(recall_stat(r["call_rate"][k], y), g, boot)
                row[f"fpr@{f}"] = cluster_bootstrap(fpr_stat(r["call_rate"][k], y), g, boot)
                tb = r["ex_call"][k][tierB_mask_x]
                row[f"tierB_recall@{f}"] = float(np.nanmean(tb))
            rows.append(row)
        for name, fl in flags_e.items():
            t = rule_table(name, fl, y, g, boot)
            rows.append({"method": name, "recall_rule": t["recall"], "fpr_rule": t["fpr"]})
        return rows

    def flat(rows):
        recs = []
        for r in rows:
            d = {}
            for k, v in r.items():
                d[k] = fmt(v) if isinstance(v, tuple) and len(v) == 3 else v
            recs.append(d)
        return pd.DataFrame(recs)

    flags_e = {"rule_broad": broad[elig], "rule_narrow": narrow[elig]}
    tierB_x = (
        (df.loc[extra, "role"] == "class_positive") & (df.loc[extra, "pos_tier"] == "B")
    ).to_numpy()
    rows = headline(results, y_e, g_e, flags_e, "main", tierB_x, n_boot)
    # rule recall on tier B
    tb_all = ((df["role"] == "class_positive") & (df["pos_tier"] == "B")).to_numpy()
    for r in rows:
        if r["method"] == "rule_broad":
            r["tierB_recall"] = float(broad[tb_all].mean())
        if r["method"] == "rule_narrow":
            r["tierB_recall"] = float(narrow[tb_all].mean())
    # tier B recall with cluster-bootstrap intervals (never trained on; scored by models that did not see the cluster)
    gB = g_x[tierB_x]
    onesB = np.ones(int(tierB_x.sum()), int)
    brows = []
    for mdl, r in results.items():
        for k, f in enumerate(FPRS):
            rate = np.nan_to_num(r["ex_call"][k][tierB_x], nan=0.0)
            brows.append(
                {
                    "method": mdl,
                    "fpr_target": f,
                    "tierB_recall": fmt(cluster_bootstrap(recall_stat(rate, onesB), gB, n_boot)),
                }
            )
    for name, fl in (("rule_broad", broad), ("rule_narrow", narrow)):
        brows.append({"method": name, "fpr_target": "rule", "tierB_recall": fmt(cluster_bootstrap(recall_stat(fl[extra][tierB_x].astype(float), onesB), gB, n_boot))})  # fmt: skip
    pd.DataFrame(brows).to_csv(out / "tierB_recall.tsv", sep="\t", index=False)
    # paired comparison with the hand rules (same clusters resampled): recall and false-positive-rate differences, model minus rule
    prow = []
    for mdl, r in results.items():
        for k, f in enumerate(FPRS):
            for name, fl in flags_e.items():
                fl = fl.astype(float)
                d_rec = cluster_bootstrap(
                    lambda i, r=r, k=k, fl=fl: recall_stat(r["call_rate"][k], y_e)(i)
                    - recall_stat(fl, y_e)(i),
                    g_e,
                    n_boot,
                )
                d_fpr = cluster_bootstrap(
                    lambda i, r=r, k=k, fl=fl: fpr_stat(r["call_rate"][k], y_e)(i)
                    - fpr_stat(fl, y_e)(i),
                    g_e,
                    n_boot,
                )
                prow.append({"model": mdl, "fpr_target": f, "rule": name, "recall_diff_model_minus_rule": fmt(d_rec), "fpr_diff_model_minus_rule": fmt(d_fpr)})  # fmt: skip
    pd.DataFrame(prow).to_csv(out / "paired_vs_rules.tsv", sep="\t", index=False)
    hm = flat(rows)
    hm.to_csv(out / "headline_metrics.tsv", sep="\t", index=False)

    # recall of model at the hand rules' own false-positive rate, read from the pooled OOF scores (optimistic, descriptive)
    matched = []
    for name, fl in flags_e.items():
        rate = fl[y_e == 0].mean()
        for mdl, r in results.items():
            t = np.quantile(r["oof"][y_e == 0], 1 - rate, method="lower")
            matched.append(
                {
                    "rule": name,
                    "rule_fpr": rate,
                    "model": mdl,
                    "model_recall_at_that_fpr_pooled_oof": float((r["oof"][y_e == 1] > t).mean()),
                    "rule_recall": float(fl[y_e == 1].mean()),
                }
            )
    pd.DataFrame(matched).to_csv(
        out / "matched_fpr_descriptive.tsv", sep="\t", index=False, float_format="%.3f"
    )

    # set-aside run (SED1 and DAN4 left out entirely)
    elig2, pos2, _, extra2 = label_sets(df, drop=NEW_COL_SED)
    y2 = pos2.astype(int)[elig2]
    g2 = df.loc[elig2, "cluster_id"].to_numpy()
    res2 = {}
    for mdl in MODELS:
        print(f"[cv set-aside] {mdl}", flush=True)
        res2[mdl] = run_cv(
            mdl, FEATURES, _arrays(df, FEATURES, elig2), y2, g2, _arrays(df, FEATURES, extra2),
            df.loc[extra2, "cluster_id"].to_numpy(), FPRS, seeds, n_jobs=n_jobs,
        )  # fmt: skip
    tb2 = (
        (df.loc[extra2, "role"] == "class_positive") & (df.loc[extra2, "pos_tier"] == "B")
    ).to_numpy()
    flat(
        headline(
            res2,
            y2,
            g2,
            {"rule_broad": broad[elig2], "rule_narrow": narrow[elig2]},
            "setaside",
            tb2,
            n_boot,
        )
    ).to_csv(out / "headline_setaside_sed1_dan4.tsv", sep="\t", index=False)

    # extras: descriptive scores by role
    xr = df.loc[extra].reset_index(drop=True)
    recs = []
    for role, idx in xr.groupby(xr["role"] + "|" + xr["pos_tier"].astype(str)).groups.items():
        idx = np.asarray(idx)
        bx, nx = broad[extra][idx], narrow[extra][idx]
        for mdl, r in results.items():
            rec = {"role|tier": role, "n": len(idx), "model": mdl, "median_score": float(np.nanmedian(r["ex_score"][idx])),
                   "rule_broad_rate": float(bx.mean()), "rule_narrow_rate": float(nx.mean())}  # fmt: skip
            for k, f in enumerate(FPRS):
                rec[f"called@{f}"] = float(np.nanmean(r["ex_call"][k][idx]))
            recs.append(rec)
    pd.DataFrame(recs).to_csv(
        out / "extra_roles_scores.tsv", sep="\t", index=False, float_format="%.3f"
    )

    # tier B recall with cluster interval, per model (extra to the headline)
    # clade hold-out
    clade = df["clade_coarse"].to_numpy()
    rows = []
    for held in ("Saccharomycotina", "non-Saccharomycotina"):
        test_m = (
            (clade == "Saccharomycotina")
            if held == "Saccharomycotina"
            else (clade != "Saccharomycotina")
        )
        tr_i = np.flatnonzero(elig & ~test_m)
        te_i = np.flatnonzero(elig & test_m)
        ytr, yte = y_all[tr_i], y_all[te_i]
        gtr, gte = df["cluster_id"].to_numpy()[tr_i], df["cluster_id"].to_numpy()[te_i]
        Xtr, Xte = df[FEATURES].to_numpy(float)[tr_i], df[FEATURES].to_numpy(float)[te_i]
        for mdl in MODELS:
            sc, cl = [], []
            for s in range(5):
                thr = inner_thresholds(mdl, FEATURES, s, Xtr, ytr, gtr, FPRS)
                v = fit_scores(mdl, FEATURES, s, Xtr, ytr, Xte)[0]
                sc.append(v)
                cl.append([v > t for t in thr])
            sc, cl = np.mean(sc, 0), np.mean(np.array(cl), 0)
            r = {"held_out": held, "model": mdl, "train_pos": int(ytr.sum()), "train_pos_clusters": len(set(gtr[ytr == 1])),
                 "test_pos": int(yte.sum()), "test_pos_clusters": len(set(gte[yte == 1])), "test_neg": int((yte == 0).sum()),
                 "auc": fmt(cluster_bootstrap(lambda i: auc(sc[i], yte[i]), gte, n_boot))}  # fmt: skip  # noqa: B023
            for k, f in enumerate(FPRS):
                r[f"recall@{f}"] = fmt(cluster_bootstrap(recall_stat(cl[k], yte), gte, n_boot))
                r[f"fpr@{f}"] = fmt(cluster_bootstrap(fpr_stat(cl[k], yte), gte, n_boot))
            rows.append(r)
        for name, fl in (("rule_broad", broad), ("rule_narrow", narrow)):
            t = rule_table(name, fl[te_i], yte, gte, n_boot)
            rows.append({"held_out": held, "model": name, "test_pos": int(yte.sum()), "test_pos_clusters": len(set(gte[yte == 1])), "test_neg": int((yte == 0).sum()), "recall@rule": fmt(t["recall"]), "fpr@rule": fmt(t["fpr"])})  # fmt: skip
    pd.DataFrame(rows).to_csv(out / "clade_holdout.tsv", sep="\t", index=False)

    # ablation by group
    abl = []
    for mdl in MODELS:
        full = run_cv(mdl, FEATURES, X_e, y_e, g_e, X_x[:0], g_x[:0], (), abl_seeds, n_jobs=n_jobs)
        arms = [("all", list(FEATURES))]
        arms += [(f"drop_{k}", [f for f in FEATURES if f not in v]) for k, v in GROUPS.items()]
        arms += [(f"only_{k}", list(v)) for k, v in GROUPS.items()]
        for arm, cols in arms:
            print(f"[ablation] {mdl} {arm}", flush=True)
            Xa = df.loc[elig, cols].to_numpy(float)
            r = (
                full
                if arm == "all"
                else run_cv(mdl, cols, Xa, y_e, g_e, Xa[:0], g_e[:0], (), abl_seeds, n_jobs=n_jobs)
            )
            t = np.quantile(r["oof"][y_e == 0], 0.95, method="lower")
            abl.append({"model": mdl, "arm": arm, "n_features": len(cols),
                        "auc": fmt(cluster_bootstrap(lambda i, r=r: auc(r["oof"][i], y_e[i]), g_e, 500)),
                        "recall_at_fpr0.05_pooled_oof_threshold": float((r["oof"][y_e == 1] > t).mean())})  # fmt: skip
    pd.DataFrame(abl).to_csv(out / "ablation.tsv", sep="\t", index=False)

    # permutation importance by group
    pg = {k: [FEATURES.index(f) for f in v] for k, v in GROUPS.items()}
    perm_rows = []
    for mdl in MODELS:
        print(f"[perm] {mdl}", flush=True)
        r = run_cv(mdl, FEATURES, X_e, y_e, g_e, X_x[:0], g_x[:0], (), abl_seeds, pg, n_jobs=n_jobs)
        for k in GROUPS:
            drop = cluster_bootstrap(
                lambda i, r=r, k=k: auc(r["oof"][i], y_e[i]) - auc(r["perm"][k][i], y_e[i]),
                g_e,
                500,
            )
            perm_rows.append(
                {"model": mdl, "group": k, "auc_drop": fmt(drop), "auc_drop_point": drop[0]}
            )
    pd.DataFrame(perm_rows).to_csv(out / "permutation_importance.tsv", sep="\t", index=False)

    # length checks
    L = df.loc[elig, "length"].to_numpy()
    big = L >= 500
    Xb = X_e[big]
    rows = []
    sc_len = np.log(L)
    sc_run = df.loc[elig, "disorder_longest_run"].to_numpy()
    resb = {}
    for mdl in MODELS:
        print(f"[length-only subset] {mdl}", flush=True)
        resb[mdl] = run_cv(
            mdl, FEATURES, Xb, y_e[big], g_e[big], Xb[:0], g_e[:0], FPRS, seeds, n_jobs=n_jobs
        )
    for subset, m in (("all", np.ones(len(y_e), bool)), (">=500 aa", big)):
        ys, gs = y_e[m], g_e[m]
        row_base = {"subset": subset, "n_pos": int(ys.sum()), "pos_clusters": len(set(gs[ys == 1])), "n_neg": int((ys == 0).sum()),
                    "median_len_pos": float(np.median(L[m][ys == 1])), "median_len_neg": float(np.median(L[m][ys == 0]))}  # fmt: skip
        for mdl in MODELS:
            sc = results[mdl]["oof"][m] if subset == "all" else resb[mdl]["oof"]
            lab = f"{mdl} (trained on all, scored here)" if subset != "all" else mdl
            for base_name, bs in (
                ("length_alone", sc_len[m]),
                ("longest_disordered_run_alone", sc_run[m]),
            ):
                d = cluster_bootstrap(
                    lambda i, sc=sc, bs=bs, ys=ys: auc(sc[i], ys[i]) - auc(bs[i], ys[i]), gs, n_boot
                )
                rows.append({**row_base, "model": mdl, "compared_with": base_name, "model_auc": fmt(cluster_bootstrap(lambda i, sc=sc, ys=ys: auc(sc[i], ys[i]), gs, 500)),
                             "baseline_auc": auc(bs, ys), "auc_difference_model_minus_baseline": fmt(d)})  # fmt: skip
            if subset != "all":
                sc2 = results[mdl]["oof"][m]
                rows.append({**row_base, "model": lab, "compared_with": "-", "model_auc": fmt(cluster_bootstrap(lambda i, sc2=sc2, ys=ys: auc(sc2[i], ys[i]), gs, 500))})  # fmt: skip
    pd.DataFrame(rows).to_csv(out / "length_check.tsv", sep="\t", index=False)

    # single-feature AUCs, orientation-free, and length proxy correlation among negatives
    srows = []
    Ln = df.loc[elig & neg, "length"]
    for f in FEATURES:
        v = df.loc[elig, f].to_numpy(float)
        ok = ~np.isnan(v)
        a_all = auc(v[ok], y_e[ok])
        okb = ok & big
        a_big = auc(v[okb], y_e[okb])
        rho = df.loc[elig & neg, f].corr(Ln, method="spearman")
        srows.append({"feature": f, "group": group_of(f), "auc_all": a_all, "auc_ge500aa": a_big, "spearman_with_length_in_negatives": rho,
                      "n_missing": int((~ok).sum())})  # fmt: skip
    pd.DataFrame(srows).to_csv(
        out / "single_feature_auc.tsv", sep="\t", index=False, float_format="%.3f"
    )

    # out-of-fold scores per protein (extras: mean over models that did not see their cluster)
    o = df[["accession", "gene", "role", "pos_tier", "clade_coarse", "cluster_id"]].copy()
    o["set"] = np.where(elig, "train_eligible", "extra_never_fitted")
    o["rule_broad"] = broad.astype(int)
    o["rule_narrow"] = narrow.astype(int)
    for mdl, r in results.items():
        col = np.full(len(df), np.nan)
        col[np.flatnonzero(elig)] = r["oof"]
        col[np.flatnonzero(extra)] = r["ex_score"]
        o[f"score_{mdl}"] = col
        for k, f in enumerate(FPRS):
            c = np.full(len(df), np.nan)
            c[np.flatnonzero(elig)] = r["call_rate"][k]
            c[np.flatnonzero(extra)] = r["ex_call"][k]
            o[f"callrate_{mdl}_fpr{f}"] = c
    o.to_csv(out / "oof_scores.tsv", sep="\t", index=False, float_format="%.4f")
    return results


# ---------------------------------------------------------------- apply
def score_new(train, new, model="rf", seed=0):
    """fit on ALL tier A plus the negatives of `train`; score `new` (a feature table). Fails on a missing feature column."""
    miss = [c for c in FEATURES if c not in new.columns]
    if miss:
        raise FeatureError(f"new table lacks feature columns: {', '.join(miss)}")
    if "accession" not in new.columns:
        raise FeatureError("new table lacks the accession column")
    miss = [c for c in FEATURES + ["role", "pos_tier"] if c not in train.columns]
    if miss:
        raise FeatureError(f"training table lacks columns: {', '.join(miss)}")
    elig, pos, _, _ = label_sets(train)
    pipe = build_pipeline(model, FEATURES, seed).fit(
        train.loc[elig, FEATURES].to_numpy(float), pos[elig].astype(int)
    )
    s = pipe.predict_proba(new[FEATURES].to_numpy(float))[:, 1]
    return pd.DataFrame({"accession": new["accession"].to_numpy(), f"score_{model}": s})


# ---------------------------------------------------------------- CLI
def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _refuse(path):
    if Path(path).exists():
        raise FeatureError(f"refusing to overwrite existing {path}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("features")
    f.add_argument("--proteins", required=True)
    f.add_argument("--predictions", required=True)
    f.add_argument("--stalk", required=True)
    f.add_argument("--rep14", required=True)
    f.add_argument("--raw", required=True, help="raw AIUPred file with residues (never committed)")
    f.add_argument("--out", required=True)
    a = sub.add_parser("analyse")
    a.add_argument("--features", required=True)
    a.add_argument("--out-dir", required=True)
    a.add_argument("--seeds", type=int, default=10)
    a.add_argument("--n-boot", type=int, default=2000)
    a.add_argument("--n-jobs", type=int, default=-1)
    p = sub.add_parser("apply")
    p.add_argument("--train-features", required=True)
    p.add_argument("--new-features", required=True)
    p.add_argument("--model", choices=MODELS, default="rf")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "features":
            _refuse(args.out)
            from stalk_features import StalkError, read_raw

            try:
                raw = read_raw(args.raw)
            except StalkError as e:
                raise FeatureError(str(e)) from None
            rd = lambda p_, **k: pd.read_csv(p_, sep="\t", dtype=str, keep_default_na=False, **k)  # noqa: E731
            rep = pd.read_csv(args.rep14, sep="\t").drop(columns=["rep_unit"], errors="ignore")
            df = build_features(rd(args.proteins), rd(args.predictions), rd(args.stalk), rep, raw)
            df.to_csv(args.out, sep="\t", index=False, float_format="%.6g")
            print(f"wrote {args.out}: {len(df)} proteins, {len(FEATURES)} features")
        elif args.cmd == "analyse":
            df = pd.read_csv(
                args.features,
                sep="\t",
                keep_default_na=False,
                na_values=[""],
                dtype={"pos_tier": str},
            )
            df["pos_tier"] = df["pos_tier"].fillna("")
            analyse(df, args.out_dir, range(args.seeds), args.n_boot, args.n_jobs)
        else:
            _refuse(args.out)

            def rt(p_):
                return pd.read_csv(
                    p_, sep="\t", keep_default_na=False, na_values=[""], dtype={"pos_tier": str}
                )

            train = rt(args.train_features)
            train["pos_tier"] = train["pos_tier"].fillna("")
            res = score_new(train, rt(args.new_features), args.model, args.seed)
            res.to_csv(args.out, sep="\t", index=False, float_format="%.5f")
            print(f"wrote {args.out}: {len(res)} proteins")
    except (FeatureError, FileNotFoundError, pd.errors.ParserError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
