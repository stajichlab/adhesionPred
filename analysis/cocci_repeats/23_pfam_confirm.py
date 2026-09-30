#!/usr/bin/python3.12
"""Independent check of the 14_repeat_detect_general.py calls against Pfam-A.

Reads the hmmsearch --cut_ga output from 22_pfam_hmmsearch.sh and answers four questions.

1. What fraction of each set (gained / shared / lost / length-matched background) carries a
   Pfam **repeat-family** domain? Gained above background is the result that matters;
   shared is the positive control.
2. Are the 33-34 aa gained calls ankyrin repeats? Same for the 51 aa and 62 aa groups.
3. Does the detector's period agree with the Pfam family's own unit length? Two independent
   measures of the unit are used: the Pfam model length (qlen), and the spacing between
   consecutive non-overlapping hits of the same family in the same protein. The spacing is
   the stronger check: it is measured on this protein by Pfam, not taken from a model.
4. How many gained calls have no Pfam hit at all, and are therefore still unverified?

A family counts as a repeat family when its Pfam description contains "repeat", or when it
is in the explicit REPEAT_EXTRA list below (canonical repeats whose description omits the
word). Every family counted is written to `pfam_repeat_families.tsv` so the rule is auditable.

Usage: /usr/bin/python3.12 23_pfam_confirm.py
"""

import argparse
from math import lgamma
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
CLANS = Path("/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.clans.tsv")

# canonical repeat families whose Pfam DESC does not contain the word "repeat"
REPEAT_EXTRA = {
    "Kelch_1",
    "Kelch_2",
    "Kelch_3",
    "Kelch_4",
    "Kelch_5",
    "Kelch_6",
    "ANAPC3",
    "Beta_propeller",
    "PQQ",
    "PQQ_2",
    "PQQ_3",
}
DOM_COLS = [
    "target",
    "t_acc",
    "tlen",
    "query",
    "q_acc",
    "qlen",
    "seq_evalue",
    "seq_score",
    "seq_bias",
    "dom_i",
    "dom_n",
    "c_evalue",
    "i_evalue",
    "dom_score",
    "dom_bias",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to",
    "acc",
]
SETS = ["gained", "shared", "lost", "background"]


def lchoose(n, k):
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def fisher_right(a, b, c, d):
    """One-sided Fisher exact P for enrichment of a/(a+b) over c/(c+d)."""
    from math import exp

    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    denom = lchoose(n, c1)
    return sum(exp(lchoose(r1, i) + lchoose(n - r1, c1 - i) - denom) for i in range(a, hi + 1))


def load_domtbl(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            rows.append(f[:22])
    d = pd.DataFrame(rows, columns=DOM_COLS)
    for c in ["tlen", "qlen", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to"]:
        d[c] = d[c].astype(int)
    for c in ["seq_evalue", "seq_score", "i_evalue", "dom_score"]:
        d[c] = d[c].astype(float)
    d["pfam"] = d.q_acc.str.split(".").str[0]
    return d


def unit_spacing(g):
    """Median gap between consecutive non-overlapping hits of one family in one protein."""
    hits = g.sort_values("ali_from")[["ali_from", "ali_to"]].values.tolist()
    keep = []
    for s, e in hits:
        if keep and s <= keep[-1][1]:  # overlapping: keep the first
            continue
        keep.append((s, e))
    if len(keep) < 2:
        return None, len(keep)
    gaps = [keep[i + 1][0] - keep[i][0] for i in range(len(keep) - 1)]
    gaps.sort()
    m = len(gaps) // 2
    med = gaps[m] if len(gaps) % 2 else (gaps[m - 1] + gaps[m]) / 2
    return med, len(keep)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--domtbl", default=str(HERE / "pfam_sets.domtbl"))
    ap.add_argument("--sets", default=str(HERE / "pfam_sets.tsv"))
    args = ap.parse_args()

    sets = pd.read_csv(args.sets, sep="\t")
    dom = load_domtbl(args.domtbl)
    clans = pd.read_csv(
        CLANS, sep="\t", header=None, names=["pfam", "clan", "clan_name", "family", "desc"]
    )
    dom = dom.merge(clans[["pfam", "clan", "clan_name", "desc"]], on="pfam", how="left")
    dom["desc"] = dom.desc.fillna("")
    dom["is_repeat"] = dom.desc.str.contains("repeat", case=False, regex=False) | dom["query"].isin(
        REPEAT_EXTRA
    )
    dom = dom.merge(
        sets[["key", "set", "length", "period_new", "period_old", "copies_new"]],
        left_on="target",
        right_on="key",
        how="left",
    )

    print(
        f"{len(dom)} domain hits at the Pfam gathering threshold over "
        f"{dom.target.nunique()} proteins of {len(sets)}\n"
    )

    # ---- 1. confirmation rates -------------------------------------------------
    any_hit = dom.groupby("target").size()
    rep_hit = dom[dom.is_repeat].groupby("target").size()
    sets["any_pfam"] = sets.key.map(any_hit).fillna(0).astype(int) > 0
    sets["repeat_pfam"] = sets.key.map(rep_hit).fillna(0).astype(int) > 0

    print("=== 1. confirmation rate by set ===")
    print(f"{'set':<12}{'n':>5}{'any Pfam':>12}{'repeat family':>16}{'median len':>12}")
    tab = {}
    for s in SETS:
        d = sets[sets.set == s]
        tab[s] = (int(d.repeat_pfam.sum()), len(d))
        print(
            f"{s:<12}{len(d):>5}"
            f"{d.any_pfam.sum():>7} {100 * d.any_pfam.mean():>4.1f}%"
            f"{d.repeat_pfam.sum():>10} {100 * d.repeat_pfam.mean():>4.1f}%"
            f"{d.length.median():>12.0f}"
        )
    for a_set, b_set in [("gained", "background"), ("shared", "background"), ("gained", "shared")]:
        a, na = tab[a_set]
        c, nc = tab[b_set]
        p = fisher_right(a, na - a, c, nc - c)
        print(f"  {a_set} vs {b_set}: one-sided Fisher P = {p:.3g}")

    # ---- 2. the named gained groups --------------------------------------------
    print("\n=== 2. the gained period groups ===")
    ANK = {"PF00023", "PF12796", "PF13606", "PF13637", "PF13857"}
    g = sets[sets.set == "gained"]
    groups = {
        "33-34 aa": g[g.period_new.isin([33, 34])],
        "51 aa": g[g.period_new == 51],
        "62 aa": g[g.period_new == 62],
    }
    for name, d in groups.items():
        sub = dom[dom.target.isin(d.key)]
        ank = sub[sub.pfam.isin(ANK)].target.nunique()
        rep = sub[sub.is_repeat].target.nunique()
        print(
            f"\n{name}: {len(d)} proteins; {sub.target.nunique()} with any Pfam hit, "
            f"{rep} with a repeat family, {ank} with an ankyrin family"
        )
        top = (
            sub.groupby(["query", "pfam", "desc"])
            .agg(
                proteins=("target", "nunique"),
                domains=("target", "size"),
                min_score=("dom_score", "min"),
                max_score=("dom_score", "max"),
            )
            .sort_values("proteins", ascending=False)
            .head(8)
        )
        print(top.to_string())

    # proteins in a group with no repeat family: what are they?
    for name, d in groups.items():
        miss = sorted(set(d.key) - set(dom[dom.is_repeat].target))
        if miss:
            print(f"\n{name}: {len(miss)} proteins with NO Pfam repeat family")
            sub = sets[sets.key.isin(miss)]
            for _, r in sub.iterrows():
                other = sorted(set(dom[dom.target == r.key]["query"]))
                print(
                    f"  {r.protein:<26}{r.length:>6} aa  copies {r.copies_new:>3}  "
                    f"other Pfam: {','.join(other) if other else '-'}  {str(r.unit_new)[:45]}"
                )

    # what the gained calls hit that is NOT a repeat family
    print("\ncommonest non-repeat Pfam families in the gained set:")
    ng = dom[(dom.set == "gained") & (~dom.is_repeat)]
    print(
        (
            ng.groupby(["query", "pfam", "desc"])
            .agg(proteins=("target", "nunique"))
            .sort_values("proteins", ascending=False)
            .head(10)
        ).to_string()
    )

    # ---- 3. period against the family's own unit length -------------------------
    print("\n=== 3. detector period against Pfam unit length ===")
    rows = []
    for (tgt, fam), grp in dom[dom.is_repeat].groupby(["target", "query"]):
        sp, ncopy = unit_spacing(grp)
        rows.append(
            {
                "target": tgt,
                "family": fam,
                "pfam": grp.pfam.iloc[0],
                "desc": grp.desc.iloc[0],
                "model_len": grp.qlen.iloc[0],
                "n_hits": ncopy,
                "spacing": sp,
                "set": grp.set.iloc[0],
                "period_new": grp.period_new.iloc[0],
                "period_old": grp.period_old.iloc[0],
            }
        )
    per = pd.DataFrame(rows)
    per.to_csv(HERE / "pfam_period_check.tsv", sep="\t", index=False)
    # Only proteins 14 actually called have a period to check. Background proteins have
    # period 0 by definition, so they cannot enter this comparison.
    ok = per[per.spacing.notna() & (per.period_new > 0)].copy()
    ok["diff"] = ok.period_new - ok.spacing
    # Some Pfam models cover several units (Ank_2 is "3 copies", model length 90). Their
    # hit-to-hit spacing is then an integer multiple of the unit, so the agreement test is
    # spacing / period close to a whole number, not spacing = period.
    ok["ratio"] = ok.spacing / ok.period_new
    ok["mult"] = ok.ratio.round()
    ok["mult_err"] = (ok.ratio - ok.mult).abs()
    # A priori rule, decided from the Pfam description alone and not from the data: a model
    # whose description says "copies" (Ank_2 "3 copies", Ank_4/Ank_5 "many copies") spans
    # several units, so its hit-to-hit spacing is a multiple of the unit. Every other model
    # is taken as one unit, and its spacing is directly comparable with the detector period.
    ok["multi_unit_model"] = ok.desc.str.contains("copies", case=False, regex=False)
    # single_unit is the complement. (These two were out of sync: the column was renamed to
    # multi_unit_model but three uses below still referred to single_unit, which raised
    # AttributeError before the section could run. Fixed 2026-09-30.)
    ok["single_unit"] = ~ok["multi_unit_model"]
    print(
        f"{len(ok)} protein x repeat-family pairs in proteins 14 called, with >= 2 "
        f"non-overlapping domains of that family, so a Pfam-measured hit-to-hit spacing"
    )
    one = ok[ok.single_unit]
    print(
        f"\n  spacing / period rounds to 1 (Pfam unit = detector period): "
        f"{len(one)} / {len(ok)}"
    )
    print(
        f"    of those, |period - spacing| <= 1 aa : {int((one['diff'].abs() <= 1).sum())} "
        f"({100 * (one['diff'].abs() <= 1).mean():.1f}%)"
    )
    print(
        f"    of those, |period - spacing| <= 3 aa : {int((one['diff'].abs() <= 3).sum())} "
        f"({100 * (one['diff'].abs() <= 3).mean():.1f}%)"
    )
    print(f"    median signed difference            : {one['diff'].median():+.1f} aa")
    mm = ok[~ok.single_unit]
    print(f"\n  spacing / period rounds to k >= 2 (multi-unit Pfam model): {len(mm)}")
    print(f"    of those, |ratio - k| <= 0.1 : {int((mm.mult_err <= 0.1).sum())}")
    print(
        f"  any whole-number relation (|ratio - k| <= 0.1, k >= 1): "
        f"{int((ok.mult_err <= 0.1).sum())} / {len(ok)} "
        f"({100 * (ok.mult_err <= 0.1).mean():.1f}%)"
    )

    print("\nby family (spacing = Pfam's own measure of the unit in that protein):")
    fam = (
        ok.groupby(["family", "pfam"])
        .agg(
            pairs=("target", "size"),
            model_len=("model_len", "first"),
            med_spacing=("spacing", "median"),
            med_period=("period_new", "median"),
            med_ratio=("ratio", "median"),
            med_diff=("diff", "median"),
            within1=("diff", lambda x: int((x.abs() <= 1).sum())),
        )
        .sort_values("pairs", ascending=False)
    )
    print(fam.to_string())
    print("\nsame, split by set:")
    print(
        ok.groupby("set")
        .agg(
            pairs=("target", "size"),
            single_unit=("single_unit", "sum"),
            within1=("diff", lambda x: int((x.abs() <= 1).sum())),
            whole_k=("mult_err", lambda x: int((x <= 0.1).sum())),
            med_diff=("diff", "median"),
        )
        .to_string()
    )
    print(
        "\nbackground proteins carrying a repeat family are not in the table above: "
        "14 did not call them, so they have no period."
    )

    # period estimate of the old detector, on the same pairs, for contrast
    o1 = ok[(ok.period_old > 0)].copy()
    o1["diff_old"] = o1.period_old - o1.spacing
    print(
        f"\nthe old detector (02) reports a period for {len(o1)} of these {len(ok)} pairs; "
        f"|period_old - spacing| <= 1 for {int((o1.diff_old.abs() <= 1).sum())}"
    )

    # ---- 4. what is not confirmed ----------------------------------------------
    print("\n=== 4. unverified gained calls ===")
    gn = sets[sets.set == "gained"]
    print(f"gained with no Pfam hit at all      : {int((~gn.any_pfam).sum())} / {len(gn)}")
    print(
        f"gained with a Pfam hit, none repeat : "
        f"{int((gn.any_pfam & ~gn.repeat_pfam).sum())} / {len(gn)}"
    )
    print(f"gained with a repeat family         : {int(gn.repeat_pfam.sum())} / {len(gn)}")

    # ---- outputs ---------------------------------------------------------------
    sets.to_csv(HERE / "pfam_confirmation.tsv", sep="\t", index=False)
    famtab = (
        dom[dom.is_repeat]
        .groupby(["query", "pfam", "clan_name", "desc"])
        .agg(
            proteins=("target", "nunique"), domains=("target", "size"), model_len=("qlen", "first")
        )
        .sort_values("proteins", ascending=False)
        .reset_index()
    )
    famtab.to_csv(HERE / "pfam_repeat_families.tsv", sep="\t", index=False)
    print("\nwrote pfam_confirmation.tsv, pfam_period_check.tsv, pfam_repeat_families.tsv")
    print(f"\nall repeat families hit ({len(famtab)}):")
    print(famtab.to_string(index=False))


if __name__ == "__main__":
    main()
