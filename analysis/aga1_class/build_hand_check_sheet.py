"""Assemble the automatic evidence table of the hand-check sample (pandas; numbers and public annotation only, no sequence).

Joins `hand_check_sample.tsv` (candidate_tiers.py) with per-protein feature values (feature_classifier.py `features` tables of the same
proteomes: composition, tail hydrophobicity, stalk, repeat z-score) and writes one row per sample protein with the tier and its reason,
the predictor calls, the model scores and rank, the repeat and composition numbers, and the UniProt domain, keyword and location
annotation. Also writes `verdicts.tsv`: the table the OWNER fills in (empty `verdict`, `confidence`, `notes`, `reviewer`, `date`).
The verdict column takes: supported, plausible, not_wall_protein, enzyme, unknown (checked by `--check`).

  python3 -I analysis/aga1_class/build_hand_check_sheet.py build --sample S --features F1 F2 ... --out-dir NEW_DIR
  python3 -I analysis/aga1_class/build_hand_check_sheet.py leads --evidence evidence_auto.tsv --leads L1 L2 ... --out-dir NEW_DIR
  python3 -I analysis/aga1_class/build_hand_check_sheet.py check --verdicts verdicts.tsv --sample S

`leads` merges the agent-drafted lead files (public UniProt and PubMed evidence per protein; never verdicts) into `evidence_leads.tsv` and writes
`evidence_sheets.md`, one section per protein with the automatic numbers and the leads; it stops if a protein has no lead row, a lead row is
for an unknown or repeated accession, a column is missing, or a PMID is not a number. `build` refuses an existing --out-dir and stops (FATAL:, exit 1) if a sample protein has no feature row. `check` stops on an accession that
is not in the sample, a duplicate, a missing sample protein or a verdict outside the vocabulary; empty verdicts are reported, not an error.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

VERDICTS = ("supported", "plausible", "not_wall_protein", "enzyme", "unknown")
FEATURE_COLS = [
    "frac_S",
    "frac_T",
    "frac_ST",
    "frac_P",
    "frac_N",
    "frac_Q",
    "frac_G",
    "frac_C",
    "ctail_hydro",
    "stalk_run_len",
    "rep_z_seq",
]
SHEET_COLS = [
    "accession", "tier", "sample", "species", "gene", "entry_name", "protein_name", "length", "tier_reason", "sp_yes", "gpi_yes",
    "disorder_longest_run", "rep_coverage", "rep_n_copies", "rep_z_seq", "frac_ST", "frac_S", "frac_T", "frac_P", "frac_N", "frac_Q",
    "frac_G", "frac_C", "ctail_hydro", "stalk_run_len", "broad_flag", "narrow_flag", "n_models_call", "score_logreg", "score_rf",
    "score_gboost", "mean_rank_species", "pfam", "interpro", "domains", "keywords", "subcellular_location", "uniprot_release",
]  # fmt: skip
VERDICT_COLS = [
    "accession",
    "tier",
    "species",
    "gene",
    "protein_name",
    "verdict",
    "confidence",
    "notes",
    "reviewer",
    "date",
]


class SheetError(Exception):
    pass


def rd(path):
    return pd.read_csv(
        path, sep="\t", keep_default_na=False, na_values=[""], dtype=str, low_memory=False
    )


def build(sample, features):
    feat = pd.concat([f[["accession"] + FEATURE_COLS] for f in features], ignore_index=True)
    if feat["accession"].duplicated().any():
        raise SheetError("an accession occurs in more than one feature table")
    miss = sorted(set(sample["accession"]) - set(feat["accession"]))
    if miss:
        raise SheetError(f"{len(miss)} sample proteins have no feature row, first: {miss[:5]}")
    out = sample.drop(columns=[c for c in FEATURE_COLS if c in sample.columns]).merge(
        feat, on="accession", how="left", validate="1:1"
    )
    cols = [c for c in SHEET_COLS if c in out.columns]
    return out[cols]


def verdict_template(sheet):
    t = sheet[["accession", "tier", "species", "gene", "protein_name"]].copy()
    for c in ("verdict", "confidence", "notes", "reviewer", "date"):
        t[c] = ""
    return t[VERDICT_COLS]


def check(verdicts, sample):
    acc = set(sample["accession"])
    if verdicts["accession"].duplicated().any():
        raise SheetError("duplicate accession in the verdicts")
    extra = sorted(set(verdicts["accession"]) - acc)
    if extra:
        raise SheetError(f"verdicts for proteins that are not in the sample: {extra[:5]}")
    miss = sorted(acc - set(verdicts["accession"]))
    if miss:
        raise SheetError(
            f"{len(miss)} sample proteins are missing from the verdicts, first: {miss[:5]}"
        )
    v = verdicts["verdict"].fillna("")
    bad = sorted(set(v[(v != "") & ~v.isin(VERDICTS)]))
    if bad:
        raise SheetError(f"verdicts outside the vocabulary {VERDICTS}: {bad}")
    return int((v == "").sum())


LEAD_COLS = [
    "accession", "gene", "family_hint", "uniprot_function_or_similarity", "go_cellular_component", "pmids", "literature_lead",
    "wall_or_enzyme_flags", "lead_confidence", "problems",
]  # fmt: skip
CONFIDENCE = ("high", "medium", "low")


def merge_leads(evidence, leads_list):
    leads = pd.concat(leads_list, ignore_index=True)
    miss = [c for c in LEAD_COLS if c not in leads.columns]
    if miss:
        raise SheetError(f"lead files lack columns: {miss}")
    if leads["accession"].duplicated().any():
        raise SheetError("duplicate accession in the lead files")
    unknown = sorted(set(leads["accession"]) - set(evidence["accession"]))
    if unknown:
        raise SheetError(f"lead rows for proteins that are not in the sample: {unknown[:5]}")
    absent = sorted(set(evidence["accession"]) - set(leads["accession"]))
    if absent:
        raise SheetError(f"{len(absent)} sample proteins have no lead row, first: {absent[:5]}")
    for acc, text in zip(leads["accession"], leads["pmids"].fillna("none"), strict=True):
        if text != "none" and not all(x.strip().isdigit() for x in text.split(";")):
            raise SheetError(
                f"{acc}: pmids must be numbers separated by ';' or 'none', got {text!r}"
            )
    bad = sorted(set(leads["lead_confidence"].fillna("")) - set(CONFIDENCE))
    if bad:
        raise SheetError(f"lead_confidence outside {CONFIDENCE}: {bad}")
    merged = evidence.merge(
        leads[LEAD_COLS].drop(columns=["gene"]), on="accession", how="left", validate="1:1"
    )
    return merged


def render_sheets(merged):
    num = lambda x, d=2: "NA" if pd.isna(x) else f"{float(x):.{d}f}"  # noqa: E731
    out = ["# Evidence sheets for the hand-check sample", "",
           "Automatic numbers come from the pipeline; the leads are drafted from public UniProt and PubMed text and are NOT verdicts.", ""]  # fmt: skip
    for r in merged.sort_values(["tier", "species", "accession"]).itertuples(index=False):
        out += [
            f"## {r.accession}: {r.gene if isinstance(r.gene, str) else ''} ({r.species}, tier {r.tier})", "",
            f"**{r.protein_name}**, {r.length} aa. Sample: {r.sample}. Tier reason: {r.tier_reason}.", "",
            f"* Predictors: signal peptide {r.sp_yes}, GPI {r.gpi_yes}, longest disordered run {r.disorder_longest_run}, tier flags broad {r.broad_flag} / narrow {r.narrow_flag}, models calling {r.n_models_call} of 3, species rank {num(r.mean_rank_species, 3)}.",
            f"* Composition: Ser+Thr {num(r.frac_ST)}, Pro {num(r.frac_P)}, Cys {num(r.frac_C)}; repeat coverage {num(r.rep_coverage)} ({num(r.rep_n_copies, 1)} copies, z {num(r.rep_z_seq, 1)}); tail hydrophobicity {num(r.ctail_hydro, 1)}.",
            f"* UniProt: Pfam {r.pfam if isinstance(r.pfam, str) else 'none'}; domains {r.domains if isinstance(r.domains, str) else 'none'}; keywords {r.keywords if isinstance(r.keywords, str) else 'none'}; location {r.subcellular_location if isinstance(r.subcellular_location, str) else 'none'}.",
            f"* Family hint: {r.family_hint if isinstance(r.family_hint, str) else 'none'}. UniProt function or similarity: {r.uniprot_function_or_similarity}. GO component: {r.go_cellular_component}.",
            f"* Literature lead ({r.lead_confidence}): {r.literature_lead} PMIDs: {r.pmids}.",
            f"* Flags: {r.wall_or_enzyme_flags}. Problems: {r.problems}.", "",
            "Verdict: _(record it in `verdicts.tsv`)_", "",
        ]  # fmt: skip
    return "\n".join(out) + "\n"


def _write(df, path):
    df.to_csv(path, sep="\t", index=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--sample", required=True)
    b.add_argument("--features", nargs="+", required=True)
    b.add_argument("--out-dir", required=True)
    m = sub.add_parser("leads")
    m.add_argument("--evidence", required=True)
    m.add_argument("--leads", nargs="+", required=True)
    m.add_argument("--out-dir", required=True)
    c = sub.add_parser("check")
    c.add_argument("--verdicts", required=True)
    c.add_argument("--sample", required=True)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "build":
            out = Path(a.out_dir)
            if out.exists():
                raise SheetError(f"{out} exists; nothing is overwritten")
            sheet = build(rd(a.sample), [rd(f) for f in a.features])
            out.mkdir(parents=True)
            _write(sheet, out / "evidence_auto.tsv")
            _write(verdict_template(sheet), out / "verdicts.tsv")
            print(f"wrote {out}: {len(sheet)} proteins")
        elif a.cmd == "leads":
            out = Path(a.out_dir)
            if out.exists():
                raise SheetError(f"{out} exists; nothing is overwritten")
            merged = merge_leads(rd(a.evidence), [rd(f) for f in a.leads])
            out.mkdir(parents=True)
            _write(merged, out / "evidence_leads.tsv")
            (out / "evidence_sheets.md").write_text(render_sheets(merged))
            print(f"wrote {out}: {len(merged)} proteins")
        else:
            empty = check(rd(a.verdicts), rd(a.sample))
            print(f"verdicts ok; {empty} still empty")
    except (SheetError, OSError, KeyError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
