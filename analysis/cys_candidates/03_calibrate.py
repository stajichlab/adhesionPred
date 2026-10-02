#!/usr/bin/env python3.12
"""Measure the feature values that the default thresholds are chosen from.

Prints three tables to stdout (TSV): (1) controls, (2) distribution of mature length,
Cys fraction and max Cys in a window over the SP-called proteins of each proteome,
(3) how many SP-called proteins pass at each tried (L, K) pair, with W fixed.
Controls are UniProt sequences with no SignalP run, so their features are computed on the
full sequence; PRA3 is also computed on the annotated ortholog after SP removal.
Read-only. Python 3.12, standard library only.

Usage: 03_calibrate.py --manifest M.tsv --control-fasta-dir DIR --pra3-id CIMG_02492-t26_1-p1
"""

import argparse
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cys_candidates as cc  # noqa: E402

CONTROLS = {
    "PRA3_Q2TVJ9_uniprot_153aa": "PRA3_Q2TVJ9.fasta",
    "Ag2_PRA_Q6QJA6": "Ag2_PRA_Q6QJA6.fasta",
    "PRA2_Q6K1L8": "PRA2_Q6K1L8.fasta",
    "RodA_P41746": "RodA_P41746.fasta",
    "CalA_Q4WXJ1": "CalA_Q4WXJ1.fasta",
    "SOWgp_Q8NK60": "SOWgp_Q8NK60.fasta",
    "CTS1_Q1E3R8": "CTS1_Q1E3R8.fasta",
}


def quantiles(values):
    v = sorted(values)
    q = statistics.quantiles(v, n=100, method="inclusive")
    return {"median": q[49], "p90": q[89], "p95": q[94], "p99": q[98], "max": v[-1]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--control-fasta-dir", required=True)
    ap.add_argument("--pra3-id", required=True, help="short id of the annotated PRA3 ortholog")
    ap.add_argument("--window", type=int, default=cc.DEFAULT_WINDOW)
    a = ap.parse_args()
    try:
        jobs = cc.read_manifest(a.manifest)
        w = a.window
        print("# table 1: controls (window W=%d)" % w)
        print("control\tlength\tmature_length\tcys\tcys_frac\tmax_cys_window\tcc\tcxc\tpest")
        for name, fn in CONTROLS.items():
            ((_, seq),) = cc.read_fasta(Path(a.control_fasta_dir) / fn)
            f = cc.features(seq, None, w)
            print(
                f"{name}\t{len(seq)}\t{f['mature_length']}\t{f['cys_count']}\t"
                f"{f['cys_frac']:.4f}\t{f['max_cys_window']}\t{f['cc_pairs']}\t{f['cxc']}\t"
                f"{f['pest_frac']:.4f}"
            )
        pools = {}
        for name, fasta, sp_path, _dom in jobs:
            recs = cc.read_fasta(fasta)
            sp = cc.read_signalp(sp_path)
            feats = []
            for header, seq in recs:
                call, _prob, cs = sp[header]
                if call != "SP":
                    continue
                f = cc.features(seq, cs, w)
                f["id"] = cc.short_id(header)
                f["length"] = len(seq)
                feats.append(f)
                if f["id"] == a.pra3_id:
                    print(
                        f"PRA3_annotated_{name}_{a.pra3_id}\t{len(seq)}\t{f['mature_length']}\t"
                        f"{f['cys_count']}\t{f['cys_frac']:.4f}\t{f['max_cys_window']}\t"
                        f"{f['cc_pairs']}\t{f['cxc']}\t{f['pest_frac']:.4f}"
                    )
            pools[name] = feats
        print("\n# table 2: SP-called proteins per proteome")
        print("proteome\tn_sp\tstat\tmature_length\tcys_frac\tmax_cys_window")
        for name, feats in pools.items():
            ml = quantiles([f["mature_length"] for f in feats])
            cf = quantiles([f["cys_frac"] for f in feats])
            mw = quantiles([f["max_cys_window"] for f in feats])
            for k in ml:
                print(f"{name}\t{len(feats)}\t{k}\t{ml[k]:.4g}\t{cf[k]:.4g}\t{mw[k]:.4g}")
        print("\n# table 3: SP-called proteins passing (mature length <= L and max Cys >= K)")
        Ls = (150, 200, 250, 300, 400)
        Ks = (6, 7, 8, 9, 10, 12)
        print("proteome\tL\t" + "\t".join(f"K>={k}" for k in Ks))
        for name, feats in pools.items():
            for L in Ls:
                counts = [
                    sum(1 for f in feats if f["mature_length"] <= L and f["max_cys_window"] >= k)
                    for k in Ks
                ]
                print(f"{name}\t{L}\t" + "\t".join(map(str, counts)))
        print("\n# table 4: window size W. max Cys in a window for the controls, and SP-called")
        print("# proteins with mature length <= 300 and max Cys >= 8 (all seven proteomes summed)")
        Ws = (30, 40, 60, 80, 100)
        ctrl = {}
        for name, fn in CONTROLS.items():
            ((_, seq),) = cc.read_fasta(Path(a.control_fasta_dir) / fn)
            ctrl[name] = seq
        print("W\t" + "\t".join(ctrl) + "\tn_pass_L300_K8")
        for ww in Ws:
            vals = [str(cc.max_cys_window(seq, ww)[0]) for seq in ctrl.values()]
            n_pass = 0
            for _name, fasta, sp_path, _dom in jobs:
                sp = cc.read_signalp(sp_path)
                for header, seq in cc.read_fasta(fasta):
                    call, _prob, cs = sp[header]
                    if call == "SP":
                        f = cc.features(seq, cs, ww)
                        n_pass += f["mature_length"] <= 300 and f["max_cys_window"] >= 8
            print(f"{ww}\t" + "\t".join(vals) + f"\t{n_pass}")
    except cc.Stop as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
