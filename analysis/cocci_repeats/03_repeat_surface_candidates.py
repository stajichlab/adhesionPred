#!/usr/bin/env python3
"""Join repeat structure + SignalP + composition into class-2a candidates.

Class 2a is repeat/avidity surface proteins (docs/TOOL-ARCHITECTURE.md). A candidate needs
all three: a tandem repeat array, a secretion signal, and the composition of one of the two
known sub-classes (Ser/Thr-rich FLO/ALS type, or Pro/Cys-rich SOWgp/BAD1 type).
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_signalp(d):
    """protein id -> (prediction, SP probability) from SignalP 6 prediction_results.txt"""
    out = {}
    f = Path(d) / "prediction_results.txt"
    if not f.exists():
        return out
    for line in open(f):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        out[p[0].split()[0]] = (p[1], p[2])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--profiles", nargs="+", required=True)
    ap.add_argument("--signalp-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-coverage", type=float, default=0.25)
    ap.add_argument("--min-copies", type=float, default=2.5)
    args = ap.parse_args()

    sp = {}
    for d in sorted(Path(args.signalp_dir).iterdir()):
        if d.is_dir():
            sp[d.name] = read_signalp(d)
    print(
        f"SignalP loaded for {len(sp)} proteomes: "
        + ", ".join(f"{k}={sum(1 for v in s.values() if v[0] != 'OTHER')}" for k, s in sp.items())
    )

    rows = []
    for prof in args.profiles:
        for r in csv.DictReader(open(prof), delimiter="\t"):
            strain = r["strain"]
            call = sp.get(strain, {}).get(r["protein"], ("", ""))
            r["signalp"] = call[0] or "not_run"
            r["signalp_prob"] = call[1]
            r["secreted"] = "yes" if call[0] and call[0] != "OTHER" else "no"
            st, pro, cys = float(r["pct_ser_thr"]), float(r["pct_pro"]), float(r["pct_cys"])
            # KNOWN DEFECT, documented 2026-09-30. Behaviour deliberately UNCHANGED here
            # because the 41- and 58-candidate tables and the retrain test in
            # analysis/model_review/repeat_structure_transfer_2a_retrain.py all depend on
            # these labels. Do not "fix" this without re-running those.
            #
            # `pro + cys > 15` is a SUM, so a high proline fraction alone satisfies it and
            # cysteine is never required. Consequences, both measured:
            #
            #  * FALSE POSITIVES: 6 of the 21 candidates labelled ProCys_rich have < 2% Cys.
            #    CIMG_04070 (the PTGIPTEWP family) has 0.0% Cys and 18.5% Pro -> 18.5 > 15
            #    -> labelled "ProCys_rich(SOWgp_BAD1_type)" despite containing no cysteine
            #    and having no relationship to SOWgp or BAD1.
            #  * FALSE NEGATIVES: BAD1 itself (UniProt A4D962, recomputed from sequence:
            #    %Ser+Thr 8.5, %Pro 3.8, %Cys 8.0, Pro+Cys 11.9) scores BELOW 15 and is
            #    classified "other". **The class named after BAD1 does not contain BAD1.**
            #    Three candidates with the same profile (Cys 5.1-6.6%, Pro+Cys 10.6-12.1)
            #    sit in "other" with it: CPOS3700_005551, VFC140_000791, CPOS1038_001811.
            #
            # The ordering also matters: SerThr is tested first, so a protein with
            # Ser+Thr just over 25% never reaches the Pro/Cys test. CIMG_04070 sits 1.8
            # points under that cutoff and 3.5 points over the Pro/Cys one.
            #
            # A corrected rule would test the two axes separately rather than summing them,
            # e.g. Pro-rich (pro > 15), Cys-rich (cys > 5), Pro/Cys-rich (both), and would
            # be calibrated so that SOWgp and BAD1 both land in the class named for them.
            # See REPORT_2026-09-30_ptgiptewp_family.md section 2.
            r["comp_class"] = (
                "SerThr_rich(FLO_ALS_type)"
                if st > 25
                else "ProCys_rich(SOWgp_BAD1_type)"
                if pro + cys > 15
                else "other"
            )
            rows.append(r)

    cand = [
        r
        for r in rows
        if float(r["rep_coverage"]) >= args.min_coverage
        and float(r["rep_n_copies"]) >= args.min_copies
        and r["secreted"] == "yes"
    ]
    print(
        f"\n{len(rows)} proteins; {len(cand)} class-2a candidates "
        f"(repeat coverage >={args.min_coverage}, >={args.min_copies} copies, secreted)"
    )

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(cand, key=lambda r: (r["strain"], -float(r["rep_coverage"]))))

    by_class = defaultdict(int)
    for r in cand:
        by_class[r["comp_class"]] += 1
    print("composition classes:", dict(by_class))
    print(
        f"\n{'strain':<32}{'protein':<20}{'len':>5}{'per':>5}{'cop':>6}{'cov':>6}"
        f"{'S+T':>6}{'P+C':>6}  class"
    )
    for r in sorted(cand, key=lambda r: (r["strain"], -float(r["rep_coverage"]))):
        pc = float(r["pct_pro"]) + float(r["pct_cys"])
        print(
            f"{r['strain'][:31]:<32}{r['protein'][:19]:<20}{r['length']:>5}{r['rep_period']:>5}"
            f"{float(r['rep_n_copies']):>6.1f}{float(r['rep_coverage']):>6.2f}"
            f"{float(r['pct_ser_thr']):>6.1f}{pc:>6.1f}  {r['comp_class']}"
        )


if __name__ == "__main__":
    main()
