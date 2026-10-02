"""Score S288C with the shipped 8M pickle under both poolings (issue: legacy pickles, new pooling).

Reads scer_legacy.npy / scer_masked.npy / scer_meta.csv written by 01_embed_esm2_8M.py
and the shipped model in models/. Reports how much the pooling change moves the scores.
Usage (in the directory holding the .npy files): python 03_pooling_mismatch.py
"""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
shipped = pickle.load(open(REPO / "models" / "adhesion_model_esm2_t6_8M_UR50D.pkl", "rb"))
meta = pd.read_csv("scer_meta.csv")
p_old = shipped.predict_proba(np.load("scer_legacy.npy"))[:, 1]
p_new = shipped.predict_proba(np.load("scer_masked.npy"))[:, 1]

call_old, call_new = p_old > 0.5, p_new > 0.5
rho = spearmanr(p_old, p_new)
out = {
    "n_proteins": len(meta),
    "spearman_rho": round(float(rho.statistic), 4),
    "max_abs_diff": round(float(np.abs(p_old - p_new).max()), 4),
    "median_abs_diff": round(float(np.median(np.abs(p_old - p_new))), 4),
    "calls_legacy_pooling": int(call_old.sum()),
    "calls_masked_pooling": int(call_new.sum()),
    "calls_in_both": int((call_old & call_new).sum()),
    "calls_only_legacy": int((call_old & ~call_new).sum()),
    "calls_only_masked": int((~call_old & call_new).sum()),
}
for k, v in out.items():
    print(f"{k}\t{v}")
pd.DataFrame(
    {
        "id": meta.id,
        "len": meta.seq.str.len(),
        "p_legacy_pooling": p_old,
        "p_masked_pooling": p_new,
    }
).to_csv("pooling_mismatch_scores.tsv", sep="\t", index=False)
pd.Series(out).to_csv("pooling_mismatch_summary.tsv", sep="\t", header=False)
