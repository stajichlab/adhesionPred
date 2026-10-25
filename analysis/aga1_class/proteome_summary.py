"""Per-species summary of the AGA1-class proteome test: hand rules, frozen classifier scores, candidates (pandas, numpy).

Inputs (all keyed by UniProt accession; none holds a sequence): the proteome feature table made by feature_classifier.py `features`
(it carries the SignalP, PredGPI, AIUPred and repeat-detector features), one score table per model made by feature_classifier.py
`apply` (fit on all tier A plus negatives of the 850), the 850-protein out-of-fold score table (to set thresholds), the 850-protein
feature table (roles of the known proteins), proteins.tsv of stage_proteomes.py (species, entry name, gene) and the UniProt FASTA
files (only the header text is read, for protein names).

Hand rules (same as feature_classifier.rule_flags): broad = signal peptide yes AND GPI yes AND longest AIUPred run >= 200;
narrow = broad AND the project repeat call (coverage >= 0.25 and copies >= 2.5). Classifier call = score at or above the score at
which the 850 NEGATIVES (hard + background, out-of-fold) have the chosen false-positive rate. The out-of-fold scores come from
models that did not see the protein, the proteome scores from the model fit on all 850, so the scales differ somewhat; the threshold
is therefore a rough calibration, and the report shows how many calls each species gets rather than reading them as precision.

  python3 -I analysis/aga1_class/proteome_summary.py --features F --scores logreg=S1 rf=S2 gboost=S3 --oof OOF --train-features TF \
      --proteins P --names FASTA [FASTA ...] --fpr 0.035 --out-dir NEW_DIR

Writes summary_by_species.tsv, name_based_recovery.tsv, known_proteins.tsv, candidates.tsv, thresholds.tsv, all_proteins_flags.tsv into
NEW_DIR (must not exist).
FATAL: and exit 1 on a missing column, an accession mismatch between tables or an existing --out-dir.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

NEG_ROLES = ("hard_negative", "background_assumed_negative")
TOP_N = 25
# UniProt protein names that mark an adhesin-like cell-wall protein. Names of unreviewed entries are often inferred from
# homology, so this is a rough, NOT independent, reference for how many of the annotated adhesin-like proteins a flag finds.
ADHESIN_NAME = re.compile(
    r"\bPA14\b|dhesin|Flocculation|gglutinin|Hyphally-regulated|Flo11|\bAls[0-9]|\bHwp|\bIff|Mucin|Cell wall adhesin",
    re.IGNORECASE,
)


class SummaryError(Exception):
    pass


def names_from_fasta(paths):
    out = {}
    for p in paths:
        for line in Path(p).read_text().splitlines():
            m = re.match(r">(?:sp|tr)\|([^|]+)\|\S+\s+(.*?)\s+OS=", line)
            if m:
                out[m.group(1)] = m.group(2)
    return out


def thresholds(oof, models, fpr):
    neg = oof[oof["role"].isin(NEG_ROLES)]
    if neg.empty:
        raise SummaryError("no negatives in the out-of-fold table")
    out = {}
    for m in models:
        col = f"score_{m}"
        if col not in oof.columns:
            raise SummaryError(f"out-of-fold table lacks {col}")
        out[m] = float(np.quantile(neg[col].to_numpy(float), 1 - fpr))
    return out


def build(feat, scores, oof, train, prot, names, fpr):
    need = [
        "accession",
        "sp_yes",
        "gpi_yes",
        "disorder_longest_run",
        "rep_coverage",
        "rep_n_copies",
        "length",
        "gene",
    ]
    miss = [c for c in need if c not in feat.columns]
    if miss:
        raise SummaryError(f"feature table lacks columns: {miss}")
    df = feat[["accession", "gene", "length", "sp_yes", "gpi_yes", "gpi_missing", "disorder_longest_run", "stalk_run_len",
               "rep_coverage", "rep_n_copies", "ctail_hydro"]].copy()  # fmt: skip
    p = prot.set_index("accession")
    if set(p.index) != set(df["accession"]):
        raise SummaryError("proteins.tsv and the feature table hold different accessions")
    df["species"] = df["accession"].map(p["species"])
    df["entry_name"] = df["accession"].map(p["entry_name"])
    df["protein_name"] = df["accession"].map(names).fillna("")
    for m, s in scores.items():
        if set(s["accession"]) != set(df["accession"]):
            raise SummaryError(f"score table {m} holds different accessions")
        df = df.merge(s, on="accession", how="left", validate="1:1")
    models = list(scores)
    thr = thresholds(oof, models, fpr)
    df["broad_flag"] = (
        (df.sp_yes == 1) & (df.gpi_yes == 1) & (df.disorder_longest_run >= 200)
    ).astype(int)
    df["repeat_call"] = ((df.rep_coverage >= 0.25) & (df.rep_n_copies >= 2.5)).astype(int)
    df["narrow_flag"] = df.broad_flag * df.repeat_call
    for m in models:
        df[f"call_{m}"] = (df[f"score_{m}"] >= thr[m]).astype(int)
    df["n_models_call"] = df[[f"call_{m}" for m in models]].sum(axis=1)
    df["mean_rank_pct"] = df[[f"score_{m}" for m in models]].rank(pct=True).mean(axis=1)
    known = train.set_index("accession")
    df["in_850_role"] = df["accession"].map(known["role"]).fillna("")
    df["in_850_tier"] = df["accession"].map(known["pos_tier"]).fillna("")
    return df, thr, models


def summarise(df, models):
    rows = []
    for sp, g in df.groupby("species", sort=True):
        row = {"species": sp, "n_proteins": len(g), "signal_peptide": int(g.sp_yes.sum()), "gpi": int(g.gpi_yes.sum()),
               "broad_flag": int(g.broad_flag.sum()), "narrow_flag": int(g.narrow_flag.sum())}  # fmt: skip
        for m in models:
            row[f"call_{m}"] = int(g[f"call_{m}"].sum())
        row["calls_2of3_models"] = int((g.n_models_call >= 2).sum())
        row["in_850"] = int((g.in_850_role != "").sum())
        a = g[(g.in_850_role == "class_positive") & (g.in_850_tier == "A")]
        row["tierA_in_proteome"] = len(a)
        row["tierA_broad"] = int(a.broad_flag.sum())
        row["tierA_narrow"] = int(a.narrow_flag.sum())
        row["tierA_2of3_models"] = int((a.n_models_call >= 2).sum())
        novel = g[g.in_850_role == ""]
        row["novel_broad"] = int(novel.broad_flag.sum())
        row["novel_narrow"] = int(novel.narrow_flag.sum())
        row["novel_2of3_models"] = int((novel.n_models_call >= 2).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def name_recovery(df, models):
    """Per species: proteins with an adhesin-like UniProt name, and how many each flag finds (flag = broad rule, narrow rule, or
    at least 2 of the models AND a signal peptide); plus how many proteins WITHOUT such a name the same flags call."""
    named = df["protein_name"].map(lambda t: bool(ADHESIN_NAME.search(t)) and "Zonadhesin" not in t)
    models_sp = (df.n_models_call >= 2) & (df.sp_yes == 1)
    flags = {"broad": df.broad_flag == 1, "narrow": df.narrow_flag == 1, "2of3_models_and_SP": models_sp,
             "broad_or_2of3_and_SP": (df.broad_flag == 1) | models_sp}  # fmt: skip
    rows = []
    for sp, g in df.groupby("species", sort=True):
        idx = g.index
        for scope, mask in (("adhesin_named", named[idx]), ("other", ~named[idx])):
            row = {"species": sp, "scope": scope, "n_proteins": int(mask.sum())}
            for k, f in flags.items():
                row[k] = int((f[idx] & mask).sum())
            rows.append(row)
    return pd.DataFrame(rows)


def candidates(df):
    novel = df[df.in_850_role == ""]
    keep = novel[(novel.broad_flag == 1) | (novel.narrow_flag == 1) | (novel.n_models_call >= 2)]
    top = novel.sort_values("mean_rank_pct", ascending=False).groupby("species").head(TOP_N)
    c = pd.concat([keep, top]).drop_duplicates("accession")
    c["why"] = np.select(
        [c.narrow_flag == 1, c.broad_flag == 1, c.n_models_call >= 2],
        ["narrow rule", "broad rule", "2 of 3 models"],
        "top 25 by mean rank",
    )
    return c.sort_values(["species", "mean_rank_pct"], ascending=[True, False])


def rd(path):
    return pd.read_csv(
        path, sep="\t", keep_default_na=False, na_values=[""], dtype={"pos_tier": str}
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--features", required=True)
    ap.add_argument("--scores", nargs="+", required=True, help="MODEL=FILE")
    ap.add_argument("--oof", required=True)
    ap.add_argument("--train-features", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--names", nargs="*", default=[])
    ap.add_argument("--fpr", type=float, default=0.035)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    try:
        out = Path(a.out_dir)
        if out.exists():
            raise SummaryError(f"{out} exists; nothing is overwritten")
        scores = {}
        for item in a.scores:
            m, _, f = item.partition("=")
            if not f:
                raise SummaryError(f"--scores expects MODEL=FILE, got {item!r}")
            scores[m] = rd(f)
        train = rd(a.train_features)
        train["pos_tier"] = train["pos_tier"].fillna("")
        df, thr, models = build(
            rd(a.features),
            scores,
            rd(a.oof),
            train,
            rd(a.proteins),
            names_from_fasta(a.names),
            a.fpr,
        )
    except (SummaryError, OSError, KeyError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    out.mkdir(parents=True)
    summ = summarise(df, models)
    summ.to_csv(out / "summary_by_species.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"model": m, "fpr_target": a.fpr, "score_threshold": t} for m, t in thr.items()]
    ).to_csv(out / "thresholds.tsv", sep="\t", index=False)
    cols = ["accession", "species", "gene", "entry_name", "protein_name", "length", "sp_yes", "gpi_yes", "disorder_longest_run",
            "rep_coverage", "rep_n_copies", "broad_flag", "narrow_flag"] + [f"score_{m}" for m in models] + [f"call_{m}" for m in models] + ["n_models_call", "mean_rank_pct"]  # fmt: skip
    df[df.in_850_role != ""][cols + ["in_850_role", "in_850_tier"]].to_csv(
        out / "known_proteins.tsv", sep="\t", index=False, float_format="%.5g"
    )
    name_recovery(df, models).to_csv(out / "name_based_recovery.tsv", sep="\t", index=False)
    c = candidates(df)
    c[cols + ["why"]].to_csv(out / "candidates.tsv", sep="\t", index=False, float_format="%.5g")
    df[cols + ["in_850_role", "in_850_tier"]].to_csv(
        out / "all_proteins_flags.tsv", sep="\t", index=False, float_format="%.5g"
    )
    print(summ.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
