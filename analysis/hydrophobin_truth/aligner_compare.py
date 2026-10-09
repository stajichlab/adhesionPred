#!/usr/bin/python3.12
"""Exploratory aligner comparison (after L5; does not change the frozen v1).

Per aligner and per search option (default filters, --nobias): per-fold cutoffs from the tuning data by the frozen rule, then
leave-cluster-out recall of T2 proteins and Pfam-missed clusters, test hard-negative call rates (all-cluster HMM), recovery of the
Jensen proteins that Pfam misses, and alignment diagnostics. The numbers come from held-out members, so they are exploratory:
a change of aligner would be a new freeze judged on data not used here.
Usage: aligner_compare.py --dir analysis/hydrophobin_truth --out docs/reports/data/sorting_hat/hydrophobin_ext
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
lco = load("lco")
ev = load("evaluate")

ARMS = {
    "mafft_linsi": {
        "models": "hmm_models",
        "tuning": "scores_pre_freeze_hmm",
        "heldout": "scores_post_freeze",
    },
    "famsa": {
        "models": "aligner_compare/hmm_famsa",
        "tuning": "aligner_compare/scores_famsa",
        "heldout": "aligner_compare/scores_famsa",
    },
    "muscle5": {
        "models": "aligner_compare/hmm_muscle5",
        "tuning": "aligner_compare/scores_muscle5",
        "heldout": "aligner_compare/scores_muscle5",
    },
}


def diagnostics(model_dir):
    rows = []
    for f in sorted(Path(model_dir).glob("fold_*.aln.fa")):
        aln = lco.read_aligned_fasta(f.read_text())
        n, w = len(aln), len(aln[0])
        slots_ok = lco.cys_slots_check(aln)
        ordinal = lco.ordinal_cys_check(aln)
        same_col = lco.cys_check(aln)
        rows.append(
            {
                "fold": f.name[:-7],
                "n": n,
                "width": w,
                "slots_ok": slots_ok,
                "ordinal": ordinal,
                "same_column": same_col,
            }
        )
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    ctx = pd.load_tuning(d)
    truth = {r["accession"]: r for r in ev.read_tsv(d / "truth_all.tsv")}
    seqs = {i: r["sequence"] for i, r in truth.items()}
    t2 = {i for i, r in truth.items() if r["tier"] == "T2"}
    r0ref = ctx["r0ref"]
    ncys = {i: s.upper().count("C") for i, s in seqs.items()}
    clusters = {}
    for r in ev.read_tsv(d / "clusters_positives.tsv"):
        if r["tier"] == "T2":
            clusters.setdefault(r["cluster"], set()).add(r["id"])
    folds = {int(r["fold"]): r for r in ev.read_tsv(d / "folds.tsv")}
    missed_clusters = [c for c, m in clusters.items() if m & ev.PFAM_MISSED]
    sc_post = d / "scores_post_freeze"
    hn_groups = {r["id"]: r["group"] for r in ev.read_tsv(sc_post / "fasta/hn_test_groups.tsv")}
    hn_seq = pd.read_fasta(sc_post / "fasta/hn_test.faa")
    hn_n = {i: pd.n_cys(s) for i, s in hn_seq.items()}
    hn_r0 = {i: r0ref.get("HN_" + i.rsplit("-", 1)[0], "") for i in hn_seq}
    lp_seq = pd.read_fasta(sc_post / "fasta/lp_scored.faa")
    lp_n = {i: pd.n_cys(s) for i, s in lp_seq.items()}
    lp_r0 = {i: r0ref.get(i, "") for i in lp_seq}
    lp_missed = {"LP_AFLA_014260", "LP_AO090012000143", "LP_ATEG_08089"}
    results = []
    for arm, cfg in ARMS.items():
        mdir, tdir, hdir = d / cfg["models"], d / cfg["tuning"], d / cfg["heldout"]
        if not (mdir / "fold_all.hmm").exists() or not list(tdir.glob("fold_all.*tblout")):
            print(arm, "not available yet; skipped")
            continue
        diag = diagnostics(mdir)
        for opt in ("default", "nobias"):
            called, cutoffs = set(), {}
            for i, f in folds.items():
                name = f"fold_{i:02d}"
                cands = []
                for p in ctx["tuning"]:
                    cands += pd.candidates(
                        p,
                        relaxed.parse_tblout(open(tdir / f"{name}.{p}.{opt}.tblout")),
                        ctx["meta"][p],
                    )
                hn = pd.hard_negative_scores(
                    relaxed.parse_tblout(open(tdir / f"{name}.hn_tuning.{opt}.tblout")),
                    ctx["hn_meta"],
                    ctx["groups"],
                )
                cut = relaxed.lowest_cutoff(
                    cands, ctx["sizes"], hn, 5, 0.05, ignore_groups={"HsbA"}
                )
                cutoffs[name] = cut
                s = ev.parse_tblout(hdir / f"hmm.{name}.test.{opt}.tblout")
                called |= ev.called_set(s, cut, r0ref, ncys) & set(f["test_members"].split(","))
            cands = []
            for p in ctx["tuning"]:
                cands += pd.candidates(
                    p,
                    relaxed.parse_tblout(open(tdir / f"fold_all.{p}.{opt}.tblout")),
                    ctx["meta"][p],
                )
            hn = pd.hard_negative_scores(
                relaxed.parse_tblout(open(tdir / f"fold_all.hn_tuning.{opt}.tblout")),
                ctx["hn_meta"],
                ctx["groups"],
            )
            cut_all = relaxed.lowest_cutoff(
                cands, ctx["sizes"], hn, 5, 0.05, ignore_groups={"HsbA"}
            )
            hn_scores = ev.parse_tblout(hdir / f"hmm.fold_all.hn_test.{opt}.tblout")
            rates = {}
            for g in sorted(set(hn_groups.values())):
                ids = [i for i, gg in hn_groups.items() if gg == g]
                rates[g] = round(
                    len(ev.called_set(hn_scores, cut_all, hn_r0, hn_n) & set(ids)) / len(ids), 4
                )
            lp_called = ev.called_set(
                ev.parse_tblout(hdir / f"hmm.fold_all.lp_scored.{opt}.tblout"), cut_all, lp_r0, lp_n
            )
            rec = ev.recovered_clusters(called, clusters, ev.PFAM_MISSED)
            results.append(
                {
                    "aligner": arm,
                    "option": opt,
                    "t2_called": len(called & t2),
                    "t2_n": len(t2),
                    "t2_recall": round(len(called & t2) / len(t2), 3),
                    "pfam_missed_clusters_recovered": len(rec),
                    "of": len(missed_clusters),
                    "pfam_missed_proteins_called": len(called & ev.PFAM_MISSED),
                    "lp_missed_recovered": len(lp_called & lp_missed),
                    "hn_test_max_rate": max(rates.values()),
                    "hn_test_rates": json.dumps(rates),
                    "cutoff_all_cluster": cut_all,
                    "median_fold_cutoff": sorted(cutoffs.values())[len(cutoffs) // 2],
                    "folds": len(diag),
                    "folds_slots_ok": sum(1 for r in diag if r["slots_ok"]),
                    "folds_ordinal_ok": sum(1 for r in diag if r["ordinal"]),
                    "folds_same_column_ok": sum(1 for r in diag if r["same_column"]),
                    "median_width": sorted(r["width"] for r in diag)[len(diag) // 2],
                }
            )
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "aligner_comparison.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    for r in results:
        print({k: v for k, v in r.items() if k not in ("hn_test_rates",)})
    return 0


if __name__ == "__main__":
    sys.exit(main())
