#!/usr/bin/env python3
"""PR 6: the SECOND specificity run, the three extra families of the owner decision of 2026-10-24.

DRAFT, not run. Stdlib only. Sibling of `pr6_family_specificity.py` (the first run, 30 families chosen by the rule
'coverage >= 0.5 OR E <= 1e-30'); that script is imported, not copied, and is not changed. Owner decision 2026-10-24: pursue
three families that the Pfam rule left out because their Pfam hits were weak partial-domain hits, each with a PANTHER family of
1,000+ fungal proteins: Big1, KNR4/SMI1, DFG10, in that order. The families come from the committed list
`pr6_family_specificity_extras_list.tsv`, NOT from the rule. The list is validated against the committed step 6.0 scan table
(accession, name, cluster, E-value and coverage of the hit).

The test is the one of the first run: hits over the 161 shared negatives, the four kept proteomes and the other Kiss positives,
with the family's OWN evidenced members as the denominator (D19: all members of the listed cluster). Each family has very few own
members (Big1 1, KNR4/SMI1 2, DFG10 1), the own-member recall is circular (`leakage=tuned_on_truth`), so the tier derived by
`cellsurface_sorting_hat.family_tier` is `smoke` at most. PANTHER columns are descriptive metadata, never a filter or truth.

THIS SCRIPT NEVER WRITES `data/sorting_hat/family_measure.tsv` OR `family_table.tsv`. `convert` here is a dry run only.

Commands:
  candidates   the extras list + scan table + clusters_missed.tsv -> the same four candidate files as the first run
               (candidate_clusters.tsv, candidate_families.tsv, family_members.tsv, families.txt) + extras_metadata.tsv.
  check-models the fetched `hmmfetch` model file must hold exactly the listed accessions and Pfam names (else it stops).
  summarize    as the first run (same hmmsearch inputs and tables); provenance.tsv names this list instead of the rule.
  convert      dry run only: prints the proposed rows' derived tiers; writes nothing.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from pathlib import Path

_BASE_PATH = Path(__file__).resolve().with_name("pr6_family_specificity.py")


def _load_base():
    mod = sys.modules.get("pr6_family_specificity")
    if mod is not None and Path(getattr(mod, "__file__", "")).resolve() == _BASE_PATH:
        return mod
    spec = importlib.util.spec_from_file_location("pr6_family_specificity", _BASE_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pr6_family_specificity"] = mod
    spec.loader.exec_module(mod)
    return mod


fs = _load_base()
ROOT = fs.ROOT
EXTRAS_LIST = Path(__file__).resolve().with_name("pr6_family_specificity_extras_list.tsv")
PANTHER_SUMMARY = ROOT / "docs/reports/data/panther_family_scoping/family_summary.tsv"
LIST_COLUMNS = (
    "order",
    "pfam_acc",
    "pfam_name",
    "cluster_id",
    "gene",
    "panther_family",
    "scan_ievalue",
    "scan_frac",
    "expected",
    "note",
    "approved_by",
)
# The owner-approved set and order (decision 2026-10-24): Big1, KNR4/SMI1, DFG10.
EXPECTED_FAMILIES = (
    ("PF27678", "BIG1_C", "KS0087"),
    ("PF09346", "SMI1_KNR4", "KS0006"),
    ("PF02544", "Steroid_dh", "KS0090"),
)
SELECTION = "owner list of 2026-10-24 (not the Pfam rule)"
EXTRAS_OUTPUTS = (*fs.CANDIDATE_OUTPUTS, "extras_metadata.tsv")


def read_extras(path: Path = EXTRAS_LIST) -> list[dict]:
    """The committed extras list, validated for shape: the three owner-approved rows in the approved order."""
    rows = fs.read_tsv(path)
    if not rows or tuple(rows[0]) != LIST_COLUMNS:
        raise ValueError(f"{path}: the columns must be {list(LIST_COLUMNS)}")
    got = tuple((r["pfam_acc"], r["pfam_name"], r["cluster_id"]) for r in rows)
    if got != EXPECTED_FAMILIES:
        raise ValueError(
            f"{path}: the rows must be exactly {list(EXPECTED_FAMILIES)} in this order, got {list(got)}"
        )
    if [r["order"] for r in rows] != ["1", "2", "3"]:
        raise ValueError(f"{path}: the order column must be 1, 2, 3")
    for r in rows:
        if r["expected"] not in {"specific", "broad_expected"} or not r["approved_by"]:
            raise ValueError(
                f"{path}: {r['pfam_acc']}: expected must be specific or broad_expected, approved_by is required"
            )
    return rows


def check_against_scan(extras: list[dict], scan_rows: list[dict], missed_rows: list[dict]) -> dict:
    """Each listed family must be a hit of its listed cluster in the committed scan table with the stated name, E-value and
    coverage. Returns {acc: {"hit": hit, "scan_row": row, "members": [accessions]}}."""
    scan = {r["cluster_id"]: r for r in scan_rows}
    missed = {r["cluster_id"]: [a for a in r["accessions"].split(",") if a] for r in missed_rows}
    found = {}
    for e in extras:
        cid, acc = e["cluster_id"], e["pfam_acc"]
        if cid not in scan or cid not in missed:
            raise ValueError(
                f"{acc}: cluster {cid} is not in the scan table and clusters_missed.tsv"
            )
        row = scan[cid]
        if row["class"] != "hit_other_pfam":
            raise ValueError(f"{acc}: cluster {cid} has class {row['class']}, not hit_other_pfam")
        hits = [h for h in fs.parse_other_pfam(row["other_pfam_families"]) if h["acc"] == acc]
        if len(hits) != 1:
            raise ValueError(
                f"{acc}: cluster {cid} has {len(hits)} scan hits of this family, not 1"
            )
        h = hits[0]
        if h["name"] != e["pfam_name"]:
            raise ValueError(
                f"{acc}: scan name {h['name']} differs from the list name {e['pfam_name']}"
            )
        if (
            h["ievalue"] != float(e["scan_ievalue"])
            or abs(h["frac"] - float(e["scan_frac"])) > 1e-9
        ):
            raise ValueError(
                f"{acc}: scan hit (E {h['ievalue']}, coverage {h['frac']}) differs from the list "
                f"(E {e['scan_ievalue']}, coverage {e['scan_frac']})"
            )
        if any(fs.rule_met(h)):
            raise ValueError(
                f"{acc}: the hit meets the Pfam rule, so it belongs to the first run, not the extras"
            )
        found[acc] = {"hit": h, "scan_row": row, "members": missed[cid]}
    return found


def panther_metadata(extras: list[dict], summary: Path = PANTHER_SUMMARY) -> list[dict]:
    """Descriptive PANTHER columns (family name, fungal size, scoping call). Never a filter, never truth."""
    table = {r["panther_family"]: r for r in fs.read_tsv(summary)}
    out = []
    for e in extras:
        row = table.get(e["panther_family"])
        if row is None:
            raise ValueError(
                f"{e['pfam_acc']}: {e['panther_family']} is not in {Path(summary).name}"
            )
        out.append(
            {
                "pfam_acc": e["pfam_acc"],
                "panther_family": e["panther_family"],
                "panther_name": row["name"],
                "panther_n_fungal": row["n_fungal"],
                "panther_broad_or_specific": row["broad_or_specific"],
            }
        )
    return out


def write_extras_candidates(
    out: Path,
    extras_list: Path = EXTRAS_LIST,
    scan: Path = fs.SCAN_RECORD,
    missed: Path = fs.MISSED,
    family_table: Path = fs.FAMILY_TABLE,
    summary: Path = PANTHER_SUMMARY,
) -> list[dict]:
    """The candidate files of the first run's format for the three listed families (same writers, same columns)."""
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"{out} exists and is not empty; nothing is overwritten")
    extras = read_extras(extras_list)
    found = check_against_scan(extras, fs.read_tsv(scan), fs.read_tsv(missed))
    meta = panther_metadata(extras, summary)
    in_table = {r["pfam_acc"] for r in fs.read_tsv(family_table) if r["pfam_acc"]}
    out.mkdir(parents=True, exist_ok=True)
    clusters, fams, members, names = [], [], [], []
    for e in extras:
        f = found[e["pfam_acc"]]
        r, h = f["scan_row"], f["hit"]
        clusters.append(
            [
                r["slice"], r["cluster_id"], r["group"], r["genes"], r["clades_fine"], r["strata_fine"],
                r["jgi_accept"], len(f["members"]), e["pfam_acc"], "owner_list_2026-10-24",
            ]
        )  # fmt: skip
        fams.append(
            [
                e["pfam_acc"], e["pfam_name"], 1, e["cluster_id"], r["slice"], 1, f"{h['frac']:.3f}",
                f"{h['ievalue']:.3g}", "yes" if fs.rule_met(h)[0] else "no",
                "yes" if fs.rule_met(h)[1] else "no", "yes" if e["pfam_acc"] in in_table else "no",
            ]
        )  # fmt: skip
        members.extend([e["pfam_acc"], e["cluster_id"], a] for a in f["members"])
        names.append(e["pfam_name"])
    fs.write_tsv(
        out / "candidate_clusters.tsv",
        ["slice", "cluster_id", "group", "genes", "clades_fine", "strata_fine", "jgi_accept",
         "n_members", "candidate_families", "rule"],
        clusters,
    )  # fmt: skip
    fs.write_tsv(
        out / "candidate_families.tsv",
        ["pfam_acc", "name", "n_clusters", "clusters", "slices", "n_hits", "max_frac", "min_ievalue",
         "rule_coverage", "rule_evalue", "in_family_table"],
        fams,
    )  # fmt: skip
    fs.write_tsv(out / "family_members.tsv", ["pfam_acc", "cluster_id", "accession"], members)
    (out / "families.txt").write_text("".join(f"{n}\n" for n in names))
    keys = [
        "pfam_acc",
        "panther_family",
        "panther_name",
        "panther_n_fungal",
        "panther_broad_or_specific",
    ]
    fs.write_tsv(
        out / "extras_metadata.tsv",
        ["order", *keys, "gene", "scan_ievalue", "scan_frac", "expected", "n_own_members"],
        [
            [e["order"], *[m[k] for k in keys], e["gene"], e["scan_ievalue"], e["scan_frac"], e["expected"],
             len(found[e["pfam_acc"]]["members"])]
            for e, m in zip(extras, meta, strict=False)
        ],
    )  # fmt: skip
    return [
        {**e, **m, "members": found[e["pfam_acc"]]["members"]}
        for e, m in zip(extras, meta, strict=False)
    ]


def check_models(models_hmm: Path, candidates: Path) -> int:
    """The fetched model file must hold exactly the candidate accessions and Pfam names (a missing name stops the job)."""
    want = {
        r["pfam_acc"]: r["name"] for r in fs.read_tsv(Path(candidates) / "candidate_families.tsv")
    }
    got, name = {}, None
    for ln in Path(models_hmm).read_text().splitlines():
        if ln.startswith("NAME"):
            name = ln.split()[1]
        elif ln.startswith("ACC"):
            if name is None:
                raise ValueError(f"{Path(models_hmm).name}: ACC line before a NAME line")
            got[ln.split()[1].split(".")[0]] = name
            name = None
    missing = sorted(want.items() - got.items())
    if missing or set(got) != set(want):
        raise ValueError(
            f"hmmfetch gave {sorted(got.items())}, expected {sorted(want.items())}; a name is missing from "
            f"Pfam-A or differs from the list: {missing}"
        )
    return len(got)


def summarize(*args, extras_list: Path = EXTRAS_LIST, **kw) -> list[dict]:
    """The first run's summarize, with provenance.tsv naming the list (and not the Pfam rule) as the selection."""
    cand = Path(args[0])
    rows_in = fs.read_tsv(cand / "candidate_families.tsv")
    if tuple((r["pfam_acc"], r["name"]) for r in rows_in) != tuple(
        (a, n) for a, n, _ in EXPECTED_FAMILIES
    ):
        raise ValueError(f"{cand}: the candidate tables are not the three extras families")
    rows = fs.summarize(*args, **kw)
    out = Path(args[6] if len(args) > 6 else kw["out"])
    prov = [
        r
        for r in fs.read_tsv(out / "provenance.tsv")
        if r["item"] not in {"script", "coverage_min", "evalue_max"}
    ]
    prov = [
        ["script", "analysis/wall_fusion_truth/pr6_family_specificity_extras.py"],
        ["selection", SELECTION],
        ["extras_list_sha256", hashlib.sha256(Path(extras_list).read_bytes()).hexdigest()],
        *([r["item"], r["value"]] for r in prov),
        ["tier_cap", "smoke (descriptive_shared_sets, own members only)"],
    ]
    fs.write_tsv(out / "provenance.tsv", ["item", "value"], prov)
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("candidates")
    c.add_argument("--list", type=Path, default=EXTRAS_LIST, dest="extras_list")
    c.add_argument("--scan", type=Path, default=fs.SCAN_RECORD)
    c.add_argument("--missed", type=Path, default=fs.MISSED)
    c.add_argument("--family-table", type=Path, default=fs.FAMILY_TABLE)
    c.add_argument("--panther-summary", type=Path, default=PANTHER_SUMMARY)
    c.add_argument("--out", type=Path, required=True)
    m = sub.add_parser("check-models")
    m.add_argument("--models", type=Path, required=True)
    m.add_argument("--candidates", type=Path, required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--candidates", type=Path, required=True)
    s.add_argument("--members-split", type=Path, default=fs.MEMBERS_SPLIT)
    s.add_argument("--neg-dir", type=Path, default=fs.NEG_DIR)
    s.add_argument("--domtbl-pos", type=Path, required=True)
    s.add_argument("--domtbl-neg", type=Path, required=True)
    s.add_argument("--proteome", type=fs.parse_proteome_arg, action="append", required=True)
    s.add_argument("--inputs", type=Path, default=fs.INPUTS_TSV)
    s.add_argument("--out", type=Path, required=True)
    v = sub.add_parser("convert")
    v.add_argument("--specificity", type=Path, required=True)
    v.add_argument("--source", required=True)
    v.add_argument("--date", required=True)
    v.add_argument("--pfam-release", default="38.2")
    v.add_argument("--hmmer-version", default="HMMER 3.4")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "candidates":
            rows = write_extras_candidates(
                a.out, a.extras_list, a.scan, a.missed, a.family_table, a.panther_summary
            )
            print(
                "3 extra families: "
                + ", ".join(f"{r['pfam_name']} ({len(r['members'])} own members)" for r in rows)
            )
        elif a.cmd == "check-models":
            print(f"{check_models(a.models, a.candidates)} models match the extras list")
        elif a.cmd == "summarize":
            given = {n: (d, f) for n, d, f in a.proteome}
            if sorted(given) != sorted(fs.KEPT_PROTEOMES):
                raise ValueError(
                    f"the proteomes must be exactly {list(fs.KEPT_PROTEOMES)} (one per species)"
                )
            rows = summarize(
                a.candidates, a.members_split, a.neg_dir, a.domtbl_pos, a.domtbl_neg,
                {n: given[n] for n in fs.KEPT_PROTEOMES}, a.out, fs.kept_proteome_counts(a.inputs),
            )  # fmt: skip
            print(f"{len(rows)} extra families summarized")
        else:
            for r in fs.convert(
                a.specificity, None, a.source, a.date, a.pfam_release, a.hmmer_version
            ):
                print(
                    f"{r['pfam_acc']}\t{r['name']}\ttier={r['_tier']}\talready_measured={r['_already_measured']}"
                )
            print("dry run: nothing written")
    except (ValueError, OSError) as err:
        print(f"pr6_family_specificity_extras: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
