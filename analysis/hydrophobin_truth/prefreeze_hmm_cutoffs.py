#!/usr/bin/python3.12
"""Cutoffs of the per-fold HMMs (and the all-cluster HMM) from the pre-freeze tuning scores. No positive score is used.

For the HMM the search option is default filters. The option rule would compare in-sample recovery of training clusters, which
ties by construction, and a tie goes to default filters. Both options' cutoffs are stored for the report.
Usage: prefreeze_hmm_cutoffs.py --dir analysis/hydrophobin_truth
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


relaxed = load("relaxed")
pd = load("prefreeze_data")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    sc = d / "scores_pre_freeze_hmm"
    ctx = pd.load_tuning(d)
    names = sorted({p.name.split(".")[0] for p in sc.glob("fold_*.tblout")})
    out = {"per10k_limit": 5, "hard_negative_rate": 0.05, "floor": relaxed.FLOOR, "folds": {}}
    for n in names:
        rec = {}
        for opt in ("default", "nobias"):
            cands = []
            for p in ctx["tuning"]:
                scores = relaxed.parse_tblout(open(sc / f"{n}.{p}.{opt}.tblout"))
                cands += pd.candidates(p, scores, ctx["meta"][p])
            hn = pd.hard_negative_scores(
                relaxed.parse_tblout(open(sc / f"{n}.hn_tuning.{opt}.tblout")),
                ctx["hn_meta"],
                ctx["groups"],
            )
            cut = relaxed.lowest_cutoff(cands, ctx["sizes"], hn, 5, 0.05, ignore_groups={"HsbA"})
            rec[opt] = {
                "cutoff": cut,
                "extra_calls": {
                    p: sum(1 for c in cands if c["proteome"] == p and c["score"] >= cut)
                    for p in ctx["tuning"]
                },
                "hn_rates": {
                    g: round(sum(1 for s in v if s >= cut) / len(v), 4) for g, v in hn.items()
                },
            }
        rec["chosen_option"] = "default"
        out["folds"][n] = rec
    json.dump(out, open(d / "prefreeze_hmm.json", "w"), indent=1)
    for n, r in out["folds"].items():
        print(n, "default", r["default"]["cutoff"], "nobias", r["nobias"]["cutoff"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
