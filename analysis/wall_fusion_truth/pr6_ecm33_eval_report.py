#!/usr/bin/env python3
"""Independent-set report of the frozen custom Ecm33 'all' model: recall per comparator, u, the rule, per clade, tier and subfamily.

DRAFT, not run on HPCC. Stdlib only; `python3 -I`. It reads the tables of pr6_ecm33_eval_members.py and of the HPCC jobs
(pr6_ecm33_eval_extract.sbatch, pr6_ecm33_eval_search.sbatch) and writes tables. It writes NO call, status, tier, family_measure or
family_table row, and no sequence: the custom model is report-only.

It REFUSES to score unless (1) analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json exists, says "frozen": true and carries an
owner_signature, (2) run.json of the search job records the sha256 of THAT file, (3) every result file listed in run.json has the recorded
sha256, (4) the cutoff of the freeze equals the committed cutoff.json. Scoring comes after the freeze; the proposal file is never accepted.

What n means. The unit is one UniRef50 cluster (its representative sequence). PANTHER is itself an HMM (a profile model built on a
phylogenetic tree), so recall here measures AGREEMENT between two profile models, not truth, and not adhesion. Every cell is `descriptive`
or `headline` (the headline column of the freeze); no interval, no `estimated`.

Columns of groups.tsv: `column` is headline | descriptive_1 | descriptive_2 (the freeze test_set.columns). `u` = clusters missed by EVERY
comparator (pfam_ga, pfam_relaxed, iteration_alignment_seeded, phmmer_seed_only, jackhmmer_seed_only). The group is testable when
u >= min_u; rule_ok when the custom model recovers at least ceil(u/2) of u.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

ITER_E = 1e-5
COMPARATORS = [
    "pfam_ga",
    "pfam_relaxed",
    "iteration_alignment_seeded",
    "phmmer_seed_only",
    "jackhmmer_seed_only",
]
COLUMNS = {
    "headline": "headline: PANTHER E <= 1e-10 and separated from every training member by cluster30",
    "descriptive_1": "descriptive: PANTHER E <= 1e-10, no separation",
    "descriptive_2": "descriptive: all new-cluster test clusters",
}
HIT_FILES = {
    "custom": "hits_custom.tsv",
    "pfam_ga": "hits_pfam_ga.tsv",
    "pfam_relaxed": "hits_pfam_relaxed.tsv",
    "iteration_alignment_seeded": "hits_iter_aln.tsv",
    "phmmer_seed_only": "hits_phmmer_seed.tsv",
    "jackhmmer_seed_only": "hits_jack_seed.tsv",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_frozen(path: Path) -> dict:
    if not Path(path).is_file():
        raise SystemExit(
            f"FATAL: freeze file {path} is missing; nothing is reported before the freeze"
        )
    fz = json.load(open(path))
    if fz.get("frozen") is not True:
        raise SystemExit(f"FATAL: freeze file {path} is not frozen")
    if not str(fz.get("owner_signature", "")).strip():
        raise SystemExit(f"FATAL: freeze file {path} has no owner_signature")
    try:
        fz["u"]["min_u"]
        fz["custom_model"]["cutoff_bits"]
        fz["searches"]["relaxed_pfam_cutoff_bits"]
        fz["test_set"]["headline_panther_evalue"]
    except KeyError as e:
        raise SystemExit(f"FATAL: freeze file lacks {e}") from e
    return fz


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r[h] for h in header])


def verify_run(results: Path, freeze_sha: str, bad: list[str]) -> dict:
    """run.json of the search job: the freeze sha256 and the sha256 of every result file."""
    p = results / "run.json"
    if not p.is_file():
        bad.append(f"{p} is missing")
        return {}
    run = json.load(open(p))
    if run.get("freeze_sha256") != freeze_sha:
        bad.append(
            "run.json was not made with this freeze file (sha256 differs): scoring before or without the freeze"
        )
    files = run.get("files", {})
    for name in list(HIT_FILES.values()):
        if name not in files:
            bad.append(f"run.json records no sha256 for {name}")
    for name, want in files.items():
        f = results / name
        if not f.is_file():
            bad.append(f"{name} is missing")
        elif sha256_file(f) != want:
            bad.append(f"{name} differs from the sha256 recorded in run.json")
    return run


def by_score(rows: list[dict], cutoff: float) -> dict[str, float]:
    return {r["id"]: float(r["score"]) for r in rows if float(r["score"]) >= cutoff}


def by_evalue(rows: list[dict], e: float = ITER_E) -> dict[str, float]:
    return {r["id"]: float(r["evalue"]) for r in rows if float(r["evalue"]) <= e}


def group_row(
    column: str,
    kind: str,
    name: str,
    ids: set[str],
    hits: dict[str, set[str]],
    min_u: int,
    n_proteins: int,
) -> dict:
    """One group of test clusters. `hits` maps `custom` and every comparator to the set of cluster ids it recovers (any id; restricted here)."""
    comp = set()
    for c in COMPARATORS:
        comp |= hits[c] & ids
    u_ids = ids - comp
    custom_u = hits["custom"] & u_ids
    u = len(u_ids)
    need = math.ceil(u / 2) if u else 0
    row = {
        "column": column,
        "group_type": kind,
        "group": name,
        "n_clusters": len(ids),
        "n_proteins": n_proteins,
        "custom": len(hits["custom"] & ids),
    }
    for c in COMPARATORS:
        row[c] = len(hits[c] & ids)
    row.update(
        {
            "any_comparator": len(comp),
            "u": u,
            "custom_recovers_of_u": len(custom_u),
            "min_u": min_u,
            "testable_u": "yes" if u >= min_u else "no",
            "ceil_u_half": need,
            "rule_ok": ("yes" if len(custom_u) >= need else "no") if u >= min_u else "not_testable",
            "recall_custom": f"{len(hits['custom'] & ids)}/{len(ids)}",
            "status": "headline" if column == "headline" else "descriptive",
        }
    )
    return row


GROUP_HEADER = (
    ["column", "group_type", "group", "n_clusters", "n_proteins", "custom"]
    + COMPARATORS
    + [
        "any_comparator",
        "u",
        "custom_recovers_of_u",
        "min_u",
        "testable_u",
        "ceil_u_half",
        "rule_ok",
        "recall_custom",
        "status",
    ]
)


def column_sets(clusters: list[dict], separated: dict[str, bool]) -> dict[str, list[dict]]:
    return {
        "headline": [
            c for c in clusters if c["headline_pe"] == "yes" and separated[c["cluster_id"]]
        ],
        "descriptive_1": [c for c in clusters if c["headline_pe"] == "yes"],
        "descriptive_2": list(clusters),
    }


def build_groups(
    clusters: list[dict], separated: dict[str, bool], hits: dict[str, set[str]], min_u: int
) -> list[dict]:
    rows = []
    for col, members in column_sets(clusters, separated).items():
        strata = [("all", "all", members)]
        for key, kind in (("clade_fine", "clade"), ("tier", "tier"), ("subfamily", "subfamily")):
            for v in sorted({c[key] for c in members}):
                strata.append((kind, v, [c for c in members if c[key] == v]))
        for kind, name, ms in strata:
            rows.append(
                group_row(
                    col,
                    kind,
                    name,
                    {c["cluster_id"] for c in ms},
                    hits,
                    min_u,
                    sum(int(c["n_proteins"]) for c in ms),
                )
            )
    return rows


CLUSTER_OUT = (
    [
        "cluster_id",
        "clade_fine",
        "tier",
        "subfamily",
        "n_proteins",
        "best_panther_evalue",
        "headline_pe",
        "separated",
        "custom_score",
        "custom_hit",
    ]
    + [f"{c}_hit" for c in COMPARATORS]
    + ["in_u"]
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--results",
        required=True,
        type=Path,
        help="output directory of pr6_ecm33_eval_search.sbatch",
    )
    ap.add_argument(
        "--clusters",
        required=True,
        type=Path,
        help="test_clusters.tsv of pr6_ecm33_eval_members.py select",
    )
    ap.add_argument(
        "--separation", required=True, type=Path, help="separation.tsv of the extract job"
    )
    ap.add_argument("--freeze", required=True, type=Path)
    ap.add_argument(
        "--cutoff-json", required=True, type=Path, help="the committed cutoff.json of the all model"
    )
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    fz = load_frozen(a.freeze)
    fz_sha = sha256_file(a.freeze)
    bad: list[str] = []
    cj = json.load(open(a.cutoff_json))
    if float(cj["cutoff"]) != float(fz["custom_model"]["cutoff_bits"]):
        bad.append("the cutoff of the freeze differs from the committed cutoff.json")
    if sha256_file(a.cutoff_json) != fz["custom_model"]["cutoff_json_sha256"]:
        bad.append("cutoff.json differs from the sha256 of the freeze")
    run = verify_run(a.results, fz_sha, bad)
    clusters = read_tsv(a.clusters)
    ids = {c["cluster_id"] for c in clusters}
    exp = fz["test_set"]
    if (
        len(clusters) != exp["expected_clusters"]
        or sum(int(c["n_proteins"]) for c in clusters) != exp["expected_proteins"]
    ):
        bad.append("the test set does not have the clusters and proteins of the freeze")
    sep = {r["cluster_id"]: r["separated"] == "yes" for r in read_tsv(a.separation)}
    if set(sep) != ids:
        bad.append("separation.tsv does not list exactly the test clusters")
    if bad:
        for b in bad:
            print("FATAL:", b, file=sys.stderr)
        return 1
    min_u = int(fz["u"]["min_u"])
    cut = float(fz["custom_model"]["cutoff_bits"])
    rel = float(fz["searches"]["relaxed_pfam_cutoff_bits"])
    scores = {
        "custom": by_score(read_tsv(a.results / HIT_FILES["custom"]), cut),
        "pfam_ga": {r["id"]: 0.0 for r in read_tsv(a.results / HIT_FILES["pfam_ga"])},
        "pfam_relaxed": by_score(read_tsv(a.results / HIT_FILES["pfam_relaxed"]), rel),
    }
    for c in ("iteration_alignment_seeded", "phmmer_seed_only", "jackhmmer_seed_only"):
        scores[c] = by_evalue(read_tsv(a.results / HIT_FILES[c]))
    hits = {k: set(v) for k, v in scores.items()}
    stray = {i for s in hits.values() for i in s} - ids
    if stray:
        print(
            f"FATAL: result ids that are not test clusters, e.g. {sorted(stray)[:3]}",
            file=sys.stderr,
        )
        return 1
    groups = build_groups(clusters, sep, hits, min_u)
    comp_all = set().union(*(hits[c] for c in COMPARATORS))
    cl_rows = []
    for c in clusters:
        i = c["cluster_id"]
        cl_rows.append(
            {
                **{
                    k: c[k]
                    for k in (
                        "cluster_id",
                        "clade_fine",
                        "tier",
                        "subfamily",
                        "n_proteins",
                        "best_panther_evalue",
                        "headline_pe",
                    )
                },
                "separated": "yes" if sep[i] else "no",
                "custom_score": scores["custom"].get(i, ""),
                "custom_hit": "yes" if i in hits["custom"] else "no",
                **{f"{k}_hit": "yes" if i in hits[k] else "no" for k in COMPARATORS},
                "in_u": "yes" if i not in comp_all else "no",
            }
        )
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "groups.tsv", GROUP_HEADER, groups)
    write_tsv(a.out / "clusters.tsv", CLUSTER_OUT, cl_rows)
    write_tsv(a.out / "u_clusters.tsv", CLUSTER_OUT, [r for r in cl_rows if r["in_u"] == "yes"])
    write_tsv(
        a.out / "not_separated.tsv",
        ["cluster_id", "clade_fine", "headline_pe"],
        [
            {k: c[k] for k in ("cluster_id", "clade_fine", "headline_pe")}
            for c in clusters
            if not sep[c["cluster_id"]]
        ],
    )
    top = {
        col: next(g for g in groups if g["column"] == col and g["group_type"] == "all")
        for col in COLUMNS
    }
    lines = [
        "PR 6 ecm33 independent-set report (DRAFT): the custom 'all' model is report-only; no call, status, tier, family_measure row or narrow_ok is written.",
        "PANTHER is itself an HMM: recall here is agreement between two profile models, not truth and not adhesion.",
        f"freeze sha256 {fz_sha}; custom cutoff {cut} bits (frozen), relaxed Pfam {rel} bits, min_u {min_u}, rule: custom recovers >= ceil(u/2) of u",
        f"search job commit: {run.get('commit', '')}",
    ]
    for col, g in top.items():
        lines.append(
            f"{col} ({COLUMNS[col]}): {g['n_clusters']} clusters ({g['n_proteins']} proteins); custom {g['recall_custom']}; "
            + "; ".join(f"{c} {g[c]}" for c in COMPARATORS)
            + f"; u {g['u']}; custom recovers {g['custom_recovers_of_u']} of u; testable {g['testable_u']}; rule {g['rule_ok']}"
        )
    (a.out / "summary.txt").write_text("\n".join(lines) + "\n")
    print((a.out / "summary.txt").read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
