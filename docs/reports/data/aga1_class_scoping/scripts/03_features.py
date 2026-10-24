#!/usr/bin/env python3
"""Per-protein features computable locally with no external tool (stdlib only), plus the two repeat detectors' numbers.

Writes features.tsv (numbers only; no sequence). Features:
  length; frac_S/T/N/Q/G/P; frac_ST; frac_STP; frac_STNQGP; n_cys
  max_ST30      highest S+T fraction in any 30-residue window
  lc_frac       fraction of residues covered by a 30-residue window with Shannon entropy < 2.5 bits (crude low-complexity)
  mean_entropy30  mean window entropy (bits)
  nterm_before_ST  residues before the first 30-residue window with S+T fraction >= 0.5 (= length if none)
  nterm_frac    that number divided by length
  longest_ST_run  length of the maximal-scoring segment with S/T = +1, other = -1 (a run at least 50 percent S/T)
  sequon_per100 N-X-S/T (X not P) sequons per 100 residues (crude N-glycosylation proxy, NOT a predictor)
  r02_*, r14_*  period, coverage, copies (and z_seq for 14) from analysis/cocci_repeats detectors 02 and 14
  call_repeat   existing `tandem_repeat_protein` rule: either detector period>0, coverage>=0.25, copies>=2.5
  call_repeat_z35  the same OR detector-14 z_seq >= 3.5 (the exploratory rule of zseq_rule_scan.tsv; not a shipped rule)
"""

import math
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import CACHE, OUT, all_sequences, read_tsv, write_tsv  # noqa: E402

W = 30
FEATS = [
    "accession", "length", "frac_S", "frac_T", "frac_N", "frac_Q", "frac_G", "frac_P", "frac_ST", "frac_STP",
    "frac_STNQGP", "n_cys", "max_ST30", "lc_frac", "mean_entropy30", "nterm_before_ST", "nterm_frac",
    "longest_ST_run", "sequon_per100", "r02_period", "r02_coverage", "r02_copies", "r14_period", "r14_coverage",
    "r14_copies", "r14_z_seq", "call_repeat", "call_repeat_z35",
]  # fmt: skip


def entropy(s):
    n = len(s)
    return -sum(c / n * math.log2(c / n) for c in Counter(s).values())


def kadane_run(seq):
    best = cur = 0
    for ch in seq:
        cur = max(0, cur + (1 if ch in "ST" else -1))
        best = max(best, cur)
    return best


def feats(seq):
    seq = seq.rstrip("*")
    n = len(seq)
    f = {"length": n}
    for aa in "STNQGP":
        f[f"frac_{aa}"] = seq.count(aa) / n
    f["frac_ST"] = f["frac_S"] + f["frac_T"]
    f["frac_STP"] = f["frac_ST"] + f["frac_P"]
    f["frac_STNQGP"] = sum(f[f"frac_{a}"] for a in "STNQGP")
    f["n_cys"] = seq.count("C")
    if n >= W:
        st = [1 if c in "ST" else 0 for c in seq]
        acc = [0]
        for v in st:
            acc.append(acc[-1] + v)
        wst = [(acc[i + W] - acc[i]) / W for i in range(n - W + 1)]
        ent = [entropy(seq[i : i + W]) for i in range(n - W + 1)]
        f["max_ST30"] = max(wst)
        f["mean_entropy30"] = sum(ent) / len(ent)
        covered = [False] * n
        for i, e in enumerate(ent):
            if e < 2.5:
                for j in range(i, i + W):
                    covered[j] = True
        f["lc_frac"] = sum(covered) / n
        first = next((i for i, v in enumerate(wst) if v >= 0.5), None)
        f["nterm_before_ST"] = n if first is None else first
    else:
        f["max_ST30"] = f["frac_ST"]
        f["mean_entropy30"] = entropy(seq)
        f["lc_frac"] = 0.0
        f["nterm_before_ST"] = n
    f["nterm_frac"] = f["nterm_before_ST"] / n
    f["longest_ST_run"] = kadane_run(seq)
    f["sequon_per100"] = 100 * len(re.findall(r"(?=N[^P][ST])", seq)) / n
    return f


def detector(path):
    out = {}
    for r in read_tsv(path):
        out[r["protein"]] = r
    return out


def called(r):
    return (
        bool(r)
        and int(r["rep_period"]) > 0
        and float(r["rep_coverage"]) >= 0.25
        and float(r["rep_n_copies"]) >= 2.5
    )


def main():
    seqs, _ = all_sequences()
    r02 = detector(CACHE / "detectors/repeat02.tsv")
    r14 = detector(CACHE / "detectors/repeat14.tsv")
    rows = []
    for p in read_tsv(OUT / "proteins.tsv"):
        a = p["accession"]
        f = feats(seqs[a])
        f["accession"] = a
        a2, a4 = r02.get(a), r14.get(a)
        for tag, r in (("r02", a2), ("r14", a4)):
            f[f"{tag}_period"] = int(r["rep_period"]) if r else 0
            f[f"{tag}_coverage"] = float(r["rep_coverage"]) if r else 0.0
            f[f"{tag}_copies"] = float(r["rep_n_copies"]) if r else 0.0
        f["r14_z_seq"] = float(a4["rep_z_seq"]) if a4 else ""
        f["call_repeat"] = int(called(a2) or called(a4))
        f["call_repeat_z35"] = int(
            bool(f["call_repeat"]) or (a4 is not None and float(a4["rep_z_seq"]) >= 3.5)
        )
        rows.append({k: (round(v, 4) if isinstance(v, float) else v) for k, v in f.items()})
    write_tsv(OUT / "features.tsv", rows, FEATS)
    print(len(rows), "rows")


if __name__ == "__main__":
    main()
