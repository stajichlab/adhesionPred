#!/usr/bin/env python3
"""PR 6 step 6.0: summary of the whole-Pfam scan of the 44 clusters the 39-model panel misses.

DRAFT, not run. Stdlib only. Run with `/usr/bin/python3.12 -I` on HPCC or `python3 -I` locally.
It reads tables and writes tables. It writes no call, no status and no tier (nothing of
`cellsurface_sorting_hat` is imported). Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md,
decision 1 (the whole-Pfam scan comes first; only clusters with no Pfam hit become custom-HMM targets).

Input tables are `hmmscan --domtblout` files of the members of the 44 missed clusters against all of
Pfam-A 38.2 (orientation: target = model, query = sequence; `pr6_pfam_scan.sbatch` makes them):
  --domtbl-ga    made with `--cut_ga`: a row is a hit that passed the model's GA threshold
  --domtbl-weak  optional, made with `-E 1e-3 --domE 1e-3` and WITHOUT `--cut_ga`: descriptive only; a family
                 listed here for a member that has no GA hit of that family is a weak hit and changes no class

Class of a cluster (a cluster is hit when any of its members is):
  hit_panel_model  a GA hit to one of the 39 panel models (the 21 wall and fusion models and the 18 active
                   family rows). A cluster is "missed" by the models of its own slice, so a hit by a model of
                   another slice is possible (the baseline shows two wall rows hit by family models). This
                   class is reported LOUDLY (stderr banner, panel_hits.tsv, summary.txt).
  hit_other_pfam   a GA hit to a Pfam family outside the panel (and no panel hit)
  no_pfam_hit      no GA hit at all: these clusters (and only these) remain custom-HMM targets (decision 1)

Proposal column (a PROPOSAL for the owner, not a decision): a `hit_other_pfam` cluster whose best GA family
covers at least PROPOSAL_FRACTION of one member's sequence is a candidate to add that Pfam model to the panel
(`add_to_panel_proposal`); a smaller cover is `partial_domain_for_owner` (the hit may be one generic domain of a
larger protein; whether it counts as covering the family is the owner's call). The fraction is the union of the
alignment intervals of the family's GA domains divided by the sequence length.

Outputs (OUTPUTS): scan_per_cluster.tsv, scan_counts.tsv (per slice and function group), custom_targets.tsv,
panel_proposals.tsv (per Pfam family), panel_hits.tsv, summary.txt, provenance.tsv. Fusion clusters are
scanned like the others; the result is descriptive (no custom fusion model is built, decision 4).

Usage: python3 -I analysis/wall_fusion_truth/pr6_pfam_scan_summary.py --domtbl-ga GA --missed M --groups G --out OUT
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WALL_MODELS = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
FAMILY_TABLE = ROOT / "data/sorting_hat/family_table.tsv"
PROPOSAL_FRACTION = 0.5  # provisional; a proposal for the owner
SLICES = ("adhesion", "wall", "fusion")
CLASSES = ("hit_panel_model", "hit_other_pfam", "no_pfam_hit")
OUTPUTS = (
    "scan_per_cluster.tsv",
    "scan_counts.tsv",
    "custom_targets.tsv",
    "panel_proposals.tsv",
    "panel_hits.tsv",
    "summary.txt",
    "provenance.tsv",
)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel_models(wall: Path = WALL_MODELS, family: Path = FAMILY_TABLE) -> set[str]:
    """The 39 panel model accessions (without version): the wall and fusion models and the active family rows."""
    out = {r["id"] for r in read_tsv(wall)}
    out |= {r["pfam_acc"] for r in read_tsv(family) if r["active"] == "yes" and r["pfam_acc"]}
    return out


def parse_hmmscan_domtbl(path: Path, ga: bool) -> dict[str, dict[str, dict]]:
    """{sequence: {family accession without version: {name, intervals, best_score, best_ievalue}}}.

    The file must end with the `# [ok]` trailer. With ga=True its option line must show --cut_ga; with
    ga=False it must not (a table made with another threshold is refused). The orientation must be
    hmmscan's (column 2 is the Pfam accession of the model, column 4 the sequence name)."""
    lines = path.read_text().splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError(f"{path.name}: no '# [ok]' trailer (hmmscan did not finish)")
    opt = [ln for ln in lines if ln.startswith("# Option settings:")]
    has_ga = any("--cut_ga" in ln for ln in opt)
    if ga and not has_ga:
        raise ValueError(f"{path.name}: the option line does not show --cut_ga")
    if not ga and has_ga:
        raise ValueError(f"{path.name}: a weak-hit table must not be made with --cut_ga")
    hits: dict[str, dict[str, dict]] = defaultdict(dict)
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) < 22:
            raise ValueError(f"{path.name}: short domain line: {ln[:80]}")
        if not f[1].startswith("PF"):
            raise ValueError(
                f"{path.name}: column 2 is not a Pfam accession (hmmsearch orientation?): {ln[:80]}"
            )
        name, acc, seq = f[0], f[1].split(".")[0], f[3]
        ievalue, score = float(f[12]), float(f[13])
        a_from, a_to = int(f[17]), int(f[18])
        fam = hits[seq].setdefault(
            acc,
            {
                "name": name,
                "intervals": [],
                "best_score": score,
                "best_ievalue": ievalue,
                "qlen": int(f[5]),
            },
        )
        fam["intervals"].append((min(a_from, a_to), max(a_from, a_to)))
        fam["best_score"] = max(fam["best_score"], score)
        fam["best_ievalue"] = min(fam["best_ievalue"], ievalue)
    return {s: dict(v) for s, v in hits.items()}


def covered_fraction(intervals: list[tuple[int, int]], qlen: int) -> float:
    """The union of the closed alignment intervals divided by the sequence length."""
    if qlen <= 0:
        return 0.0
    total, cur_a, cur_b = 0, None, None
    for a, b in sorted(intervals):
        if cur_b is None or a > cur_b + 1:
            if cur_b is not None:
                total += cur_b - cur_a + 1
            cur_a, cur_b = a, b
        else:
            cur_b = max(cur_b, b)
    if cur_b is not None:
        total += cur_b - cur_a + 1
    return min(1.0, total / qlen)


def f3(x: float) -> str:
    return f"{x:.3f}"


def classify(cluster: dict, members: list[str], ga: dict, weak: dict, panel: set[str]) -> dict:
    """The scan result of one missed cluster (the class, the families and the proposal)."""
    panel_hits: dict[str, str] = {}
    other: dict[str, dict] = {}
    weak_fams: dict[str, str] = {}
    for acc in members:
        for fam, d in ga.get(acc, {}).items():
            frac = covered_fraction(d["intervals"], d["qlen"])
            if fam in panel:
                panel_hits[fam] = d["name"]
            else:
                best = other.get(fam)
                if best is None or frac > best["frac"]:
                    other[fam] = {
                        "name": d["name"],
                        "frac": frac,
                        "score": d["best_score"],
                        "ievalue": d["best_ievalue"],
                        "seq": acc,
                    }
        for fam, d in weak.get(acc, {}).items():
            if fam not in ga.get(acc, {}):
                weak_fams[fam] = d["name"]
    if panel_hits:
        cls, proposal = "hit_panel_model", "already_in_panel"
    elif other:
        cls = "hit_other_pfam"
        best_frac = max(o["frac"] for o in other.values())
        proposal = (
            "add_to_panel_proposal"
            if best_frac >= PROPOSAL_FRACTION
            else "partial_domain_for_owner"
        )
    else:
        cls, proposal = "no_pfam_hit", "custom_target"
    return {
        "class": cls,
        "proposal": proposal,
        "panel": panel_hits,
        "other": other,
        "weak": weak_fams,
    }


def fam_text(d: dict[str, str]) -> str:
    return ";".join(f"{a}:{n}" for a, n in sorted(d.items()))


def other_text(d: dict[str, dict]) -> str:
    return ";".join(
        f"{a}:{o['name']}:frac={f3(o['frac'])}:score={o['score']:g}:ievalue={o['ievalue']:g}:seq={o['seq']}"
        for a, o in sorted(d.items(), key=lambda kv: (-kv[1]["frac"], kv[0]))
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--domtbl-ga", type=Path, required=True)
    ap.add_argument("--domtbl-weak", type=Path)
    ap.add_argument(
        "--missed", type=Path, required=True, help="clusters_missed.tsv of the Task 5.4 record"
    )
    ap.add_argument(
        "--groups", type=Path, required=True, help="candidate_groups.tsv of the PR 6 spec"
    )
    ap.add_argument("--wall-models", type=Path, default=WALL_MODELS)
    ap.add_argument("--family-table", type=Path, default=FAMILY_TABLE)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.out.exists() and any(a.out.iterdir()):
        raise SystemExit(f"{a.out} is not empty; nothing is overwritten")

    panel = panel_models(a.wall_models, a.family_table)
    missed = read_tsv(a.missed)
    groups = {r["cluster_id"]: r["group"] for r in read_tsv(a.groups)}
    expected = {acc for r in missed for acc in r["accessions"].split(",")}
    if any(r["cluster_id"] not in groups for r in missed):
        raise SystemExit("a missed cluster has no function group in the groups table")
    ga = parse_hmmscan_domtbl(a.domtbl_ga, ga=True)
    weak = parse_hmmscan_domtbl(a.domtbl_weak, ga=False) if a.domtbl_weak else {}
    for name, table in (("GA", ga), ("weak", weak)):
        extra = sorted(set(table) - expected)
        if extra:
            raise SystemExit(
                f"{name} table: hits for sequences that are not members of a missed cluster: {extra[:5]}"
            )

    results = []
    for r in sorted(missed, key=lambda r: (SLICES.index(r["slice"]), r["cluster_id"])):
        members = r["accessions"].split(",")
        res = classify(r, members, ga, weak, panel)
        results.append((r, members, res))

    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        a.out / "scan_per_cluster.tsv",
        [
            "slice",
            "cluster_id",
            "group",
            "genes",
            "clades_fine",
            "strata_fine",
            "jgi_accept",
            "n_members",
            "class",
            "proposal",
            "panel_models",
            "other_pfam_families",
            "weak_families",
        ],
        [
            [
                r["slice"],
                r["cluster_id"],
                groups[r["cluster_id"]],
                r["genes"],
                r["clades_fine"],
                r["strata_fine"],
                r["jgi_accept"],
                len(m),
                res["class"],
                res["proposal"],
                fam_text(res["panel"]),
                other_text(res["other"]),
                fam_text(res["weak"]),
            ]
            for r, m, res in results
        ],
    )

    counts: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: dict.fromkeys(CLASSES, 0))
    for r, _, res in results:
        counts[(r["slice"], groups[r["cluster_id"]])][res["class"]] += 1
        counts[(r["slice"], "ALL")][res["class"]] += 1
    write_tsv(
        a.out / "scan_counts.tsv",
        ["slice", "group", "n_clusters", *CLASSES],
        [
            [s, g, sum(c.values()), *(c[k] for k in CLASSES)]
            for (s, g), c in sorted(
                counts.items(), key=lambda kv: (SLICES.index(kv[0][0]), kv[0][1] != "ALL", kv[0][1])
            )
        ],
    )

    write_tsv(
        a.out / "custom_targets.tsv",
        [
            "slice",
            "cluster_id",
            "group",
            "genes",
            "clades_fine",
            "strata_fine",
            "accessions",
            "weak_families",
        ],
        [
            [
                r["slice"],
                r["cluster_id"],
                groups[r["cluster_id"]],
                r["genes"],
                r["clades_fine"],
                r["strata_fine"],
                r["accessions"],
                fam_text(res["weak"]),
            ]
            for r, _, res in results
            if res["class"] == "no_pfam_hit"
        ],
    )

    fam_rows: dict[str, dict] = {}
    for r, _, res in results:
        if res["class"] != "hit_other_pfam" or res["proposal"] != "add_to_panel_proposal":
            continue
        for fam, o in res["other"].items():
            if o["frac"] < PROPOSAL_FRACTION:
                continue
            d = fam_rows.setdefault(
                fam, {"name": o["name"], "clusters": [], "slices": set(), "fracs": []}
            )
            d["clusters"].append(r["cluster_id"])
            d["slices"].add(r["slice"])
            d["fracs"].append(o["frac"])
    write_tsv(
        a.out / "panel_proposals.tsv",
        [
            "pfam_acc",
            "name",
            "n_clusters",
            "clusters",
            "slices",
            "min_fraction",
            "max_fraction",
            "status",
        ],
        [
            [
                fam,
                d["name"],
                len(d["clusters"]),
                ";".join(d["clusters"]),
                ";".join(sorted(d["slices"])),
                f3(min(d["fracs"])),
                f3(max(d["fracs"])),
                f"proposal_for_owner (cover >= {PROPOSAL_FRACTION})",
            ]
            for fam, d in sorted(fam_rows.items())
        ],
    )

    panel_rows = [
        [
            r["slice"],
            r["cluster_id"],
            r["genes"],
            r["accessions"],
            fam_text(res["panel"]),
            r["clades_fine"],
        ]
        for r, _, res in results
        if res["class"] == "hit_panel_model"
    ]
    write_tsv(
        a.out / "panel_hits.tsv",
        ["slice", "cluster_id", "genes", "accessions", "panel_models", "clades_fine"],
        panel_rows,
    )

    n = {k: sum(1 for _, _, res in results if res["class"] == k) for k in CLASSES}
    lines = [
        "PR 6 step 6.0: whole-Pfam scan of the clusters missed by the 39-model panel (DRAFT summary)",
        f"missed clusters: {len(results)} ({', '.join(f'{s} {sum(1 for r, _, _ in results if r['slice'] == s)}' for s in SLICES)})",
        *(f"{k}: {n[k]}" for k in CLASSES),
        f"custom-HMM targets (no_pfam_hit): {n['no_pfam_hit']}",
        f"panel proposals (families, cover >= {PROPOSAL_FRACTION}): {len(fam_rows)}",
        "descriptive only: no call, no status, no tier; the proposals are for the owner to decide",
    ]
    if panel_rows:
        lines.append(
            f"REPORT LOUD: {len(panel_rows)} missed cluster(s) are hit by a panel model: "
            + "; ".join(f"{p[1]} ({p[2]}: {p[4]})" for p in panel_rows)
        )
    (a.out / "summary.txt").write_text("\n".join(lines) + "\n")
    write_tsv(
        a.out / "provenance.tsv",
        ["item", "value"],
        [
            ["domtbl_ga_sha256", sha256(a.domtbl_ga)],
            ["domtbl_weak_sha256", sha256(a.domtbl_weak) if a.domtbl_weak else "none"],
            ["missed_sha256", sha256(a.missed)],
            ["groups_sha256", sha256(a.groups)],
            ["proposal_fraction", str(PROPOSAL_FRACTION)],
            ["panel_models", str(len(panel))],
        ],
    )
    print("\n".join(lines))
    if panel_rows:
        print("!" * 70, file=sys.stderr)
        print(
            f"REPORT LOUD: {len(panel_rows)} missed cluster(s) are hit by a panel model (see panel_hits.tsv)",
            file=sys.stderr,
        )
        print("!" * 70, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
