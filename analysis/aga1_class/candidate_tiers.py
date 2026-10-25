"""Tiered candidate list for the AGA1 class and the seeded hand-check sample (pandas, numpy).

Owner decisions of 2026-10-26 (memory note aga1-next-steps-decisions-2026-10-26):
  Tier 1 = narrow rule, OR (broad rule AND at least 2 of 3 models call the protein).
  Tier 2 = not tier 1, and (broad rule, OR (at least 2 of 3 models AND a signal peptide)).
  Tier 3 = not tier 1 or 2, a signal peptide, and among the 25 highest mean model ranks of its species.
  Tiers 1 and 2 together are capped at 30 novel proteins per species: tier 1 is never cut, tier 2 is cut by mean rank, and the cut
  proteins move to tier 3 with the reason "tier 2 overflow".
  Hand-check sample: ALL novel tier 1 proteins, plus a seeded sample of 25 novel tier 2 and 15 novel tier 3 proteins, spread over the
  species by taking them in turn in a sha256 order (seed + tier + accession), so the sample is reproducible and needs no random library.
  "Novel" = not among the 850 labelled proteins (those keep their tier for reference but are not reviewed again).

The mean rank is recomputed WITHIN each species from the three model scores (percentile of each score among that species' proteins,
averaged), so it does not depend on which proteomes were scored in the same run.

  python3 -I analysis/aga1_class/candidate_tiers.py assign --flags F1 F2 F3 --out tiers.tsv --accessions-out acc.txt --record-out rec.txt
  python3 -I analysis/aga1_class/candidate_tiers.py annotate --tiers tiers.tsv --annotations annotations.tsv --out candidates.tsv

Refuses an existing output. FATAL: and exit 1 on a missing column, a duplicate accession across the flag tables, or annotations that do not
cover the candidates.
"""

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

SEED = 20261026
CAP = 30
TOP3 = 25
N_SAMPLE = {"2": 25, "3": 15}
MODELS = ("logreg", "rf", "gboost")
NEED = [
    "accession",
    "species",
    "sp_yes",
    "gpi_yes",
    "broad_flag",
    "narrow_flag",
    "n_models_call",
    "in_850_role",
] + [f"score_{m}" for m in MODELS]
KEEP = [
    "accession", "species", "gene", "entry_name", "protein_name", "length", "sp_yes", "gpi_yes", "disorder_longest_run",
    "rep_coverage", "rep_n_copies", "broad_flag", "narrow_flag", "n_models_call", "score_logreg", "score_rf", "score_gboost",
    "mean_rank_species", "in_850_role", "in_850_tier",
]  # fmt: skip


class TierError(Exception):
    pass


def hash_key(seed, tier, accession):
    return hashlib.sha256(f"{seed}:{tier}:{accession}".encode()).hexdigest()


def spread_sample(df, tier, n, seed):
    """n accessions of df (one tier), species taken in turn, each species in sha256 order; all of them if fewer than n."""
    if len(df) <= n:
        return list(df["accession"])
    order = {}
    for sp, g in df.groupby("species", sort=True):
        order[sp] = sorted(g["accession"], key=lambda a: hash_key(seed, tier, a))
    species = sorted(order)
    out, i = [], 0
    while len(out) < n:
        for sp in species:
            if i < len(order[sp]) and len(out) < n:
                out.append(order[sp][i])
        i += 1
    return out


def assign(flags, seed=SEED):
    miss = [c for c in NEED if c not in flags.columns]
    if miss:
        raise TierError(f"flag tables lack columns: {miss}")
    if flags["accession"].duplicated().any():
        raise TierError("an accession occurs in more than one flag table or twice in one")
    df = flags.copy()
    for m in MODELS:
        df[f"_r_{m}"] = df.groupby("species")[f"score_{m}"].rank(pct=True)
    df["mean_rank_species"] = df[[f"_r_{m}" for m in MODELS]].mean(axis=1)
    broad, narrow, models = df.broad_flag == 1, df.narrow_flag == 1, df.n_models_call
    sp = df.sp_yes == 1
    novel = df.in_850_role.fillna("") == ""
    t1 = narrow | (broad & (models >= 2))
    t2 = ~t1 & (broad | ((models >= 2) & sp))
    df["tier"], df["tier_reason"] = "", ""
    df.loc[t1, ["tier", "tier_reason"]] = ["1", "narrow rule, or broad rule with 2 of 3 models"]
    df.loc[t2, ["tier", "tier_reason"]] = [
        "2",
        "broad rule alone, or 2 of 3 models with a signal peptide",
    ]
    for _species, g in df[novel].groupby("species"):
        over = len(g[g.tier == "1"]) + len(g[g.tier == "2"]) - CAP
        if over > 0:
            cut = (
                g[g.tier == "2"]
                .sort_values(["mean_rank_species", "accession"], ascending=[True, True])
                .head(over)
            )
            df.loc[cut.index, ["tier", "tier_reason"]] = [
                "3",
                "tier 2 overflow (cap 30 per species)",
            ]
    free = novel & sp & (df.tier == "")
    for _species, g in df[free].groupby("species"):
        top = g.sort_values(["mean_rank_species", "accession"], ascending=[False, True]).head(TOP3)
        df.loc[top.index, ["tier", "tier_reason"]] = [
            "3",
            f"top {TOP3} by mean model rank with a signal peptide, no flag",
        ]
    df["review_eligible"] = (novel & (df.tier != "")).astype(int)
    df["sample"] = ""
    t1n = df[(df.tier == "1") & novel]
    df.loc[t1n.index, "sample"] = "all tier 1"
    for t, n in N_SAMPLE.items():
        pool = df[(df.tier == t) & novel]
        for acc in spread_sample(pool, t, n, seed):
            df.loc[df.accession == acc, "sample"] = f"sampled tier {t} (seed {seed})"
    cols = [c for c in KEEP if c in df.columns] + [
        "tier",
        "tier_reason",
        "review_eligible",
        "sample",
    ]
    return (
        df[df.tier != ""][cols]
        .sort_values(["tier", "species", "mean_rank_species"], ascending=[True, True, False])
        .reset_index(drop=True)
    )


def record(tiers, seed):
    lines = [
        f"seed {seed}; cap {CAP} novel proteins per species for tiers 1 and 2; tier 3 = top {TOP3} per species; sample sizes {N_SAMPLE}"
    ]
    for t in ("1", "2", "3"):
        g = tiers[tiers.tier == t]
        n = g[g.review_eligible == 1]
        lines.append(
            f"tier {t}: {len(g)} proteins, {len(n)} novel, {int((g['sample'] != '').sum())} in the hand-check sample"
        )
        for sp, h in g.groupby("species"):
            lines.append(
                f"  {sp}: {len(h)} proteins, {int(h.review_eligible.sum())} novel, {int((h['sample'] != '').sum())} sampled"
            )
    return "\n".join(lines) + "\n"


def annotate(tiers, ann):
    miss = sorted(set(tiers.accession) - set(ann.accession))
    if miss:
        raise TierError(f"annotations lack {len(miss)} candidates, first: {miss[:5]}")
    if ann.accession.duplicated().any():
        raise TierError("duplicate accession in the annotations")
    cols = [
        "accession",
        "pfam",
        "interpro",
        "domains",
        "keywords",
        "subcellular_location",
        "reviewed",
        "uniprot_release",
    ]
    out = tiers.merge(ann[cols], on="accession", how="left", validate="1:1")
    return out


def _refuse(*paths):
    for p in paths:
        if p and Path(p).exists():
            raise TierError(f"{p} exists; nothing is overwritten")


def rd(path):
    return pd.read_csv(
        path,
        sep="\t",
        keep_default_na=False,
        na_values=[""],
        dtype={"tier": str, "in_850_role": str, "in_850_tier": str},
        low_memory=False,
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("assign")
    a.add_argument("--flags", nargs="+", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--accessions-out", required=True)
    a.add_argument("--record-out", required=True)
    a.add_argument("--seed", type=int, default=SEED)
    n = sub.add_parser("annotate")
    n.add_argument("--tiers", required=True)
    n.add_argument("--annotations", required=True)
    n.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "assign":
            _refuse(args.out, args.accessions_out, args.record_out)
            flags = pd.concat([rd(p) for p in args.flags], ignore_index=True)
            tiers = assign(flags, args.seed)
            tiers.to_csv(args.out, sep="\t", index=False, float_format="%.5g")
            Path(args.accessions_out).write_text("".join(f"{x}\n" for x in tiers.accession))
            Path(args.record_out).write_text(record(tiers, args.seed))
            print(Path(args.record_out).read_text(), end="")
        else:
            _refuse(args.out)
            tiers = rd(args.tiers)
            annotate(tiers, rd(args.annotations)).to_csv(args.out, sep="\t", index=False)
            print(f"wrote {args.out}")
    except (TierError, FileNotFoundError, KeyError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
