#!/usr/bin/python3.12
"""Compare the old (02) and new (14) repeat calls on the real Coccidioides proteins.

Inputs: `repeat_profile_{longread,reference}.tsv` (from 02) and
`repeat_general_{longread,reference}.tsv` (from 17, which runs 14 on the same proteins).

A call is a protein passing the thresholds 03_repeat_surface_candidates.py uses:
coverage >= 0.25 and copies >= 2.5. The report also uses the stricter pair from the
02 summary line (coverage >= 0.3, copies >= 3), because the 2026-09-27 report's "58 with a
substantial repeat" is that number.

More calls is not better. Low-complexity regions are the false-positive mode, so every
gained call is printed with the Wootton-Federhen entropy of its repeat region, the fraction
of the region taken by its three commonest residues, and the region sequence, so it can be
read by eye.

Outputs: `repeat_real_compare.tsv` (per-protein, both detectors, with a gain/loss/shared
label) and a printed spot-check block.

Usage: ./18_repeat_real_compare.py [--top 25]
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
COV, COP = 0.25, 2.5
STRICT_COV, STRICT_COP = 0.3, 3.0


def load(paths):
    d = pd.concat([pd.read_csv(p, sep="\t") for p in paths], ignore_index=True)
    d["key"] = d.strain.astype(str) + "|" + d.protein.astype(str)
    return d


def called(d, cov, cop):
    return set(d[(d.rep_coverage >= cov) & (d.rep_n_copies >= cop) & (d.rep_period > 0)].key)


def top3_frac(s):
    if not s:
        return 0.0
    c = Counter(s)
    return sum(sorted(c.values(), reverse=True)[:3]) / len(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--out", default=str(HERE / "repeat_real_compare.tsv"))
    args = ap.parse_args()

    old = load([HERE / "repeat_profile_longread.tsv", HERE / "repeat_profile_reference.tsv"])
    new = load([HERE / "repeat_general_longread.tsv", HERE / "repeat_general_reference.tsv"])
    print(f"old rows {len(old)}, new rows {len(new)}", file=sys.stderr)
    shared_keys = set(old.key) & set(new.key)
    print(f"proteins in both: {len(shared_keys)}", file=sys.stderr)

    for name, (cov, cop) in {
        "03 thresholds (cov>=0.25, copies>=2.5)": (COV, COP),
        "02 summary thresholds (cov>=0.3, copies>=3)": (STRICT_COV, STRICT_COP),
    }.items():
        o, n = called(old, cov, cop), called(new, cov, cop)
        o &= shared_keys
        n &= shared_keys
        print(
            f"\n{name}\n  old {len(o)}  new {len(n)}  shared {len(o & n)}  "
            f"gained {len(n - o)}  lost {len(o - n)}"
        )

    o, n = called(old, COV, COP) & shared_keys, called(new, COV, COP) & shared_keys
    # suffix every column by hand: pandas only suffixes names the two frames share, and
    # the new detector adds columns the old one has no counterpart for.
    o2 = old.rename(columns={c: f"{c}_old" for c in old.columns if c != "key"})
    n2 = new.rename(columns={c: f"{c}_new" for c in new.columns if c != "key"})
    m = o2.merge(n2, on="key")
    m["call"] = [
        "shared" if (k in o and k in n) else "gained" if k in n else "lost" if k in o else "neither"
        for k in m.key
    ]
    m["region_top3"] = [round(top3_frac(str(u)), 3) for u in m.rep_unit_new.fillna("")]
    keep = [
        "strain_old",
        "protein_old",
        "length_old",
        "call",
        "rep_period_old",
        "rep_n_copies_old",
        "rep_coverage_old",
        "rep_period_new",
        "rep_n_copies_new",
        "rep_coverage_new",
        "rep_score_new",
        "rep_region_entropy_new",
        "rep_unit_period_new",
        "region_top3",
        "pct_ser_thr_new",
        "pct_pro_new",
        "pct_cys_new",
        "top3_frac_new",
        "rep_unit_new",
    ]
    out = m[m.call != "neither"][keep].sort_values(
        ["call", "rep_coverage_new"], ascending=[True, False]
    )
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {len(out)} rows to {args.out}", file=sys.stderr)

    g = m[m.call == "gained"].copy()
    lc = g[(g.rep_region_entropy_new < 2.5) | (g.region_top3 > 0.6)]
    print(
        f"\ngained calls: {len(g)}; of these {len(lc)} look low-complexity "
        f"(region entropy < 2.5 bits or top-3 residues > 60% of the unit)"
    )
    print(f"\nSPOT CHECK: {min(args.top, len(g))} gained calls, highest new coverage first")
    hdr = f"{'protein':<24}{'len':>5}{'per':>5}{'cop':>5}{'cov':>6}{'H':>6}{'t3':>6}  unit"
    print(hdr)
    for _, r in g.sort_values("rep_coverage_new", ascending=False).head(args.top).iterrows():
        print(
            f"{str(r.protein_old)[:23]:<24}{r.length_old:>5}{r.rep_period_new:>5}"
            f"{r.rep_n_copies_new:>5}{r.rep_coverage_new:>6.2f}"
            f"{r.rep_region_entropy_new:>6.2f}{r.region_top3:>6.2f}  {str(r.rep_unit_new)[:70]}"
        )

    lost = m[m.call == "lost"]
    print(f"\nLOST calls: {len(lost)}")
    print(hdr.replace("unit", "old unit"))
    for _, r in lost.sort_values("rep_coverage_old", ascending=False).head(args.top).iterrows():
        print(
            f"{str(r.protein_old)[:23]:<24}{r.length_old:>5}{r.rep_period_old:>5}"
            f"{r.rep_n_copies_old:>5.1f}{r.rep_coverage_old:>6.2f}"
            f"{r.rep_region_entropy_new:>6.2f}{r.region_top3:>6.2f}  {str(r.rep_unit_old)[:70]}"
        )


if __name__ == "__main__":
    main()
