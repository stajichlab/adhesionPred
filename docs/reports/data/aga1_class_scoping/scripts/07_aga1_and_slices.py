#!/usr/bin/env python3
"""AGA1 next to named Ser/Thr-rich proteins, matched operating point of the feature model vs the existing repeat call,
and the narrow-class slices. Reads the committed tables written by scripts 01-06 plus data/sorting_hat/family_table.tsv.
Writes aga1_vs_named.tsv, matched_operating_point.tsv, narrow_slices.tsv. Stdlib only.

Narrow-slice rule (provisional, owner wording "at least 5 positives in at least 2 clades"): reported in two readings,
  strict = at least 5 positives IN EACH of at least 2 clades; weak = at least 5 positives in total over at least 2 clades.
The report states what the data say; it does not decide classes.
"""

import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module  # noqa: E402

from common import OUT, ROOT, read_tsv, write_tsv  # noqa: E402

inv = import_module("05_inventory")
NAMED = [  # (label, accession)
    ("AGA1", "P32323"), ("SED1", "Q01589"), ("DAN4", "P47179"), ("MSB2", "P32334"), ("HKR1", "P41809"),
    ("EGT2", "P42835"), ("HPF1", "Q05164"), ("FLO1", "P32768"), ("FLO5", "P38894"), ("FLO9", "P39712"),
    ("FLO10", "P36170"), ("FLO11", "P08640"), ("SAG1", "P20840"), ("HSP150", "P32478"), ("PIR1", "Q03178"),
    ("TIR1", "P10863"),
]  # fmt: skip


def wilson(k, n):
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(c - h, 3), round(c + h, 3)


def main():
    P = inv.load()
    S = {a: set(inv.subsets(p)) for a, p in P.items()}
    oof = defaultdict(dict)
    for r in read_tsv(OUT / "oof_scores.tsv"):
        oof[r["scheme"]][r["accession"]] = float(r["score"])
    hn_ids = [a for a in P if "hn_all" in S[a]]

    # ---- AGA1 vs named
    rows = []
    for label, a in NAMED:
        p = P[a]
        row = {"label": label, "accession": a, "role": p["role"], "curated_class": p["curated_class"], "evidence": p["evidence"],
               "uniprot_gpi_feature": p["uniprot_gpi_feature"], "uniprot_signal_feature": p["uniprot_signal_feature"]}  # fmt: skip
        for k in ("length", "frac_S", "frac_T", "frac_ST", "frac_STNQGP", "max_ST30", "longest_ST_run", "nterm_before_ST", "lc_frac", "n_cys", "sequon_per100", "r02_coverage", "r14_coverage", "r14_z_seq", "call_repeat", "call_repeat_z35"):  # fmt: skip
            row[k] = p[k]
        lco = oof["LCO"].get(a)
        lcl = oof["LCL_seed0"].get(a)
        row["lco_score"] = "" if lco is None else lco
        # rank among hard negatives of the SAME clade fold (Saccharomycotina) for LCO; among all hard negatives for LCL
        sac_hn = [x for x in hn_ids if P[x]["clade_fine"] == p["clade_fine"] and x in oof["LCO"]]
        row["lco_hn_in_clade_scoring_higher"] = "" if lco is None else f"{sum(oof['LCO'][x] > lco for x in sac_hn if x != a)} of {len([x for x in sac_hn if x != a])}"  # fmt: skip
        row["lcl_score"] = "" if lcl is None else lcl
        row["lcl_hn_scoring_higher"] = "" if lcl is None else f"{sum(oof['LCL_seed0'][x] > lcl for x in hn_ids if x != a)} of {len([x for x in hn_ids if x != a])}"  # fmt: skip
        rows.append(row)
    cols = list(rows[0])
    write_tsv(OUT / "aga1_vs_named.tsv", rows, cols)

    # ---- matched operating point (scores pooled over 5 fold models of leave-cluster-out seed 0; threshold chosen on the same scores)
    sc = oof["LCL_seed0"]
    rule = "call_repeat"
    fp_call = sum(int(P[a][rule]) for a in hn_ids)
    hn_sorted = sorted((sc[a] for a in hn_ids), reverse=True)
    thr = (
        hn_sorted[fp_call] if fp_call < len(hn_sorted) else -1e9
    )  # strictly above thr => at most fp_call false positives
    mp = []
    for ps in ("pos_tierA", "pos_all"):
        pos = [a for a in P if ps in S[a]]
        k_model = sum(sc[a] > thr for a in pos)
        k_call = sum(int(P[a][rule]) for a in pos)
        for who, k, fp in (("existing repeat call", k_call, fp_call), ("feature model (comp_len_struct, LCL seed 0 OOF), threshold set so hard-negative false positives <= the call's", k_model, sum(sc[a] > thr for a in hn_ids))):  # fmt: skip
            lo, hi = wilson(k, len(pos))
            mp.append({"positives": ps, "method": who, "n_pos": len(pos), "pos_called": k, "sensitivity": round(k / len(pos), 3), "wilson_lo": lo, "wilson_hi": hi, "hn_all_false_pos": fp, "hn_all_n": len(hn_ids)})  # fmt: skip
    write_tsv(OUT / "matched_operating_point.tsv", mp, list(mp[0]))

    # ---- narrow slices
    fam = {r["pfam_acc"]: r for r in read_tsv(ROOT / "data/sorting_hat/family_table.tsv")}
    pos = [a for a in P if "pos_all" in S[a]]
    by = defaultdict(list)
    for a in pos:
        by[P[a]["narrow_class"] or "(no domain/name rule)"].append(a)
    out = []
    for cls, ids in sorted(by.items()):
        A = [a for a in ids if P[a]["pos_tier"] == "A"]
        clA = Counter(P[a]["clade_fine"] for a in A)
        clAll = Counter(P[a]["clade_fine"] for a in ids)
        pf = Counter()
        for a in ids:
            for x in {y for y in P[a]["pfam_curated"].replace(",", ";").split(";") if y}:
                pf[x] += 1
        pftxt = "; ".join(f"{x} {fam[x]['name']} ({pf[x]}/{len(ids)}; active={fam[x]['active']}, tier={fam[x]['specificity_tier']}, narrow_ok={fam[x]['narrow_ok']})" for x in sorted(pf) if x in fam and pf[x] / len(ids) >= 0.5) or "none at 50% of members"  # fmt: skip
        strict = lambda c: sum(1 for v in c.values() if v >= 5) >= 2  # noqa: E731
        weak = lambda c: sum(c.values()) >= 5 and len(c) >= 2  # noqa: E731
        famAll = Counter(P[a]["clade_family"] or "(unknown)" for a in ids)
        clades_all = set(clAll)
        conf = [a for a in hn_ids if P[a]["clade_fine"] in clades_all]
        conf_st = [a for a in conf if "hn_st30" in S[a]]
        med = statistics.median(sc[a] for a in ids)
        high = sorted((a for a in conf if sc[a] >= med), key=lambda a: -sc[a])
        top = ", ".join(f"{P[a]['gene'] or a}({P[a]['hn_stratum']})" for a in high[:3])
        called = sum(int(P[a]["call_repeat"]) for a in ids)
        out.append({
            "narrow_class": cls, "positives_all": len(ids), "positives_tierA": len(A),
            "clades_tierA": "; ".join(f"{k}:{v}" for k, v in sorted(clA.items())) or "-",
            "clades_all": "; ".join(f"{k}:{v}" for k, v in sorted(clAll.items())),
            "families_all": "; ".join(f"{k}:{v}" for k, v in sorted(famAll.items())),
            "rule_strict_all_at_family_level": "yes" if strict(famAll) else "no",
            "clusters_all": len({P[a]["cluster"] for a in ids}), "clusters_tierA": len({P[a]["cluster"] for a in A}),
            "with_uniprot_gpi_feature": sum(P[a]["uniprot_gpi_feature"] == "yes" for a in ids),
            "pfam_narrow_call_available": pftxt,
            "rule_strict_tierA": "yes" if strict(clA) else "no", "rule_strict_all": "yes" if strict(clAll) else "no",
            "rule_weak_all": "yes" if weak(clAll) else "no",
            "repeat_call_sensitivity": f"{called}/{len(ids)}",
            "hard_negatives_same_clades": len(conf), "st_rich_hard_negatives_same_clades": len(conf_st),
            "hn_scoring_at_or_above_slice_median": len(high), "top_scoring_hard_negatives": top,
        })  # fmt: skip
    write_tsv(OUT / "narrow_slices.tsv", out, list(out[0]))
    print(len(out), "slices")


if __name__ == "__main__":
    main()
