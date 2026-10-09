#!/usr/bin/python3.12
"""Cutoffs of the relaxed Pfam level from the pre-freeze scores (tuning data only), per search option, and the option rule.

The cutoff rule uses no positive score. The T2 positives are scored only by choose_option, to choose between the two
search options (recovered Pfam-missed clusters at each option's cutoff). Writes prefreeze_relaxed.json.
Usage: prefreeze_cutoffs.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
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
PER10K = 5
HN_RATE = 0.05
OPTIONS = ("default", "nobias")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    sc = d / "scores_pre_freeze"
    ctx = pd.load_tuning(d)
    tuning, sizes, meta, groups, hn_meta, r0ref = (
        ctx[k] for k in ("tuning", "sizes", "meta", "groups", "hn_meta", "r0ref")
    )
    t2 = pd.read_fasta(sc / "t2_positives.faa")
    folds = list(csv.DictReader(open(d / "folds.tsv"), delimiter="\t"))
    missed_clusters = {
        f["held_out_cluster"]: set(f["test_members"].split(","))
        for f in folds
        if f["pfam_missed"] == "1"
    }
    pfam_missed = set(
        "A0A0A2KI06 A0A6V8R0V1 G9MKB2 I1RDP9 I1RDQ0 O94217 Q00367 Q0KKA0 Q7Z9L4".split()
    )
    out = {
        "per10k_limit": PER10K,
        "hard_negative_rate": HN_RATE,
        "floor": relaxed.FLOOR,
        "tuning_proteomes": tuning,
        "options": {},
    }
    recovered = {}
    for opt in OPTIONS:
        cands = []
        for p in tuning:
            scores = relaxed.parse_tblout(open(sc / f"{p}.{opt}.tblout"))
            cands += pd.candidates(p, scores, meta[p])
        hn_scores = pd.hard_negative_scores(
            relaxed.parse_tblout(open(sc / f"hn_tuning.{opt}.tblout")), hn_meta, groups
        )
        cutoff = relaxed.lowest_cutoff(
            cands, sizes, hn_scores, PER10K, HN_RATE, ignore_groups={"HsbA"}
        )
        per_proteome = {
            p: sum(1 for c in cands if c["proteome"] == p and c["score"] >= cutoff) for p in tuning
        }
        hn_rates = {
            g: round(sum(1 for s in v if s >= cutoff) / len(v), 4) for g, v in hn_scores.items()
        }
        t2s = relaxed.parse_tblout(open(sc / f"t2_positives.{opt}.tblout"))
        called_t2 = {
            i
            for i, s in t2s.items()
            if s >= cutoff and r0ref.get(i, "") == "called" and pd.n_cys(t2[i]) >= 8
        }
        rec = {c for c, members in missed_clusters.items() if members & called_t2 & pfam_missed}
        recovered[opt] = rec
        out["options"][opt] = {
            "cutoff": cutoff,
            "extra_calls_per_tuning_proteome": per_proteome,
            "tuning_proteome_sizes": sizes,
            "hard_negative_call_rates": hn_rates,
            "n_candidates_at_floor": len(cands),
            "pfam_missed_clusters_recovered": sorted(rec),
            "pfam_missed_proteins_called": sorted(called_t2 & pfam_missed),
        }
    out["chosen_option"] = relaxed.choose_option(recovered)
    out["note"] = (
        "T2 positive scores were used only to choose between the two options; no positive score set a cutoff."
    )
    json.dump(out, open(d / "prefreeze_relaxed.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "options"}, indent=1))
    for opt, v in out["options"].items():
        print(
            opt,
            "cutoff",
            v["cutoff"],
            "| recovered Pfam-missed clusters",
            len(v["pfam_missed_clusters_recovered"]),
            "| extra calls",
            v["extra_calls_per_tuning_proteome"],
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
