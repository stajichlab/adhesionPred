"""Compare CV schemes and feature sets for the adhesion classifier."""

import importlib.util
import os
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location(
    "features", str(REPO) + "/src/surface_glyco/features.py"
)
_f = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(_f)
extract_sequence_features = _f.extract_sequence_features

meta = pd.read_csv("train_meta.csv")
clu = pd.read_csv("clu_cluster.tsv", sep="\t", header=None, names=["rep", "id"]).drop_duplicates(
    "id"
)
meta = meta.merge(clu, on="id", how="left")
y = meta.label.values
groups = meta.rep.astype("category").cat.codes.values

X = {
    "length_only": np.log(meta.seq.str.len().values)[:, None],
    "aa_comp+length": np.array([extract_sequence_features(s) for s in meta.seq]),
    "esm2_8M_legacy_pool": np.load("train_legacy.npy"),
    "esm2_8M_masked_pool": np.load("train_masked.npy"),
}


def clf():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=1.0))


def run(Xf, splits):
    p = np.full(len(y), np.nan)
    for tr, te in splits:
        m = clf().fit(Xf[tr], y[tr])
        p[te] = m.predict_proba(Xf[te])[:, 1]
    ok = ~np.isnan(p)
    pred = p[ok] > 0.5
    return {
        "roc_auc": roc_auc_score(y[ok], p[ok]),
        "pr_auc": average_precision_score(y[ok], p[ok]),
        "recall": (pred & (y[ok] == 1)).sum() / (y[ok] == 1).sum(),
        "fpr": (pred & (y[ok] == 0)).sum() / (y[ok] == 0).sum(),
    }


# leave-family-out: hold out one positive family + one negative species per fold
fam = {"ALS1_homologs": 0, "FLO_homologs": 1, "FLO11_Scer": 2}
neg = {
    "FungiDB-68_CalbicansSC5314_AnnotatedProteins": 0,
    "FungiDB-68_Spombe972h_AnnotatedProteins": 1,
    "FungiDB-68_CneoformansJEC21_AnnotatedProteins": 2,
}
fold = meta.source.map({**fam, **neg}).fillna(-1).values  # Calb/Scer singletons always train
lfo = [(np.where(fold != k)[0], np.where(fold == k)[0]) for k in range(3)]

schemes = {
    "random_5fold (current)": list(
        StratifiedKFold(5, shuffle=True, random_state=0).split(X["length_only"], y)
    ),
    "homology_grouped_5fold": list(
        StratifiedGroupKFold(5, shuffle=True, random_state=0).split(X["length_only"], y, groups)
    ),
    "leave_family_out": lfo,
}
out = []
for fn, Xf in X.items():
    for sn, sp in schemes.items():
        out.append(dict(features=fn, cv=sn, **run(Xf, sp)))
res = pd.DataFrame(out)
print(res.round(3).to_string(index=False))
res.round(4).to_csv("cv_results.csv", index=False)

# Length-matched hard negatives: can the model rank positives above LONG negatives?
long_neg = (y == 0) & (meta.seq.str.len().values > 900)
print("\nlong negatives (>900 aa):", long_neg.sum())


if os.path.exists("scer_legacy.npy"):
    shipped = pickle.load(open(str(REPO) + "/models/adhesion_model_esm2_t6_8M_UR50D.pkl", "rb"))
    sc = pd.read_csv("scer_meta.csv")
    sc["len"] = sc.seq.str.len()
    sc["p_shipped"] = shipped.predict_proba(np.load("scer_legacy.npy"))[:, 1]
    sc["p_retrained"] = (
        clf().fit(X["esm2_8M_masked_pool"], y).predict_proba(np.load("scer_masked.npy"))[:, 1]
    )
    sc["p_aacomp"] = (
        clf()
        .fit(X["aa_comp+length"], y)
        .predict_proba(np.array([extract_sequence_features(s) for s in sc.seq]))[:, 1]
    )
    adh = {
        "YAR050W": "FLO1",
        "YHR211W": "FLO5",
        "YAL063C": "FLO9",
        "YKR102W": "FLO10",
        "YIR019C": "FLO11",
        "YNR044W": "AGA1",
        "YGL032C": "AGA2",
        "YJR004C": "SAG1",
        "YCR089W": "FIG2",
    }
    hard = {
        "YMR307W": "GAS1",
        "YGR189C": "CRH1",
        "YKL096W": "CWP1",
        "YKL096W-A": "CWP2",
        "YDR077W": "SED1",
        "YER011W": "TIR1",
        "YJR150C": "DAN1",
        "YKL164C": "PIR1",
        "YJL159W": "HSP150",
        "YBR078W": "ECM33",
        "YLR120C": "YPS1",
        "YLR300W": "EXG1",
        "YGR014W": "MSB2",
        "YDR420W": "HKR1",
        "YLR110C": "CCW12",
        "YJL158C": "CIS3",
        "YIL011W": "TIR3",
        "YDR134C": "CCW22",
        "YOR009W": "TIR4",
        "YLR194C": "NCW2",
    }
    sc["id0"] = sc.id.str.split().str[0]
    for col in ["p_shipped", "p_retrained", "p_aacomp"]:
        called = sc[col] > 0.5
        print(
            f"\n[{col}] S288C ORFs={len(sc)} called={called.sum()} ({called.mean():.2%}); "
            f"known adhesins called {sc[sc.id0.isin(adh)][col].gt(0.5).sum()}/{len(adh)}; "
            f"hard negatives called {sc[sc.id0.isin(hard)][col].gt(0.5).sum()}/{len(hard)}; "
            f"median len called={sc[called].len.median()} vs all={sc.len.median()}; "
            f"spearman(score,len)={sc[col].corr(np.log(sc.len), method='spearman'):.2f}"
        )
    t = sc[sc.id0.isin(adh) | sc.id0.isin(hard)].copy()
    t["gene"] = t.id0.map({**adh, **hard})
    t["class"] = np.where(t.id0.isin(adh), "adhesin", "hard_neg")
    print(
        t[["gene", "class", "len", "p_shipped", "p_retrained", "p_aacomp"]]
        .sort_values(["class", "gene"])
        .round(3)
        .to_string(index=False)
    )
    sc.drop(columns="seq").to_csv("scer_scores.csv", index=False)
