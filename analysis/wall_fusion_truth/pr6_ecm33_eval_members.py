#!/usr/bin/env python3
"""Independent-set evaluation of the frozen custom Ecm33 'all' model: the test set, its id lists, the cluster30 separation and the purity check.

DRAFT, not run on HPCC. Stdlib only; `python3 -I`. It reads committed tables (candidates.tsv of the scoping, members.tsv of the build) and
writes tables and id lists. It writes NO sequence, no call, no status, no tier, no family_measure row (the custom model is report-only).
Scoping: docs/reports/2026-10-23-ecm33-independent-set-scoping.md; candidates: docs/reports/data/ecm33_independent_set_candidates/.

Test set (every choice is a PROPOSAL until the owner signs analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json)
  - candidates.tsv rows with overlap_with_training == new_cluster: 725 proteins in 187 UniRef50 clusters (recounted here, not copied).
  - UNIT = one UniRef50 cluster. Its test sequence is the cluster REPRESENTATIVE (the sequence UniRef50 gives for the cluster id), as the
    training members are representatives. The representative need not be one of the PANTHER proteins of the cluster.
  - HEADLINE = clusters with at least one protein of PANTHER E <= 1e-10 (563 proteins, 109 clusters) AND no training member in the same
    MMseqs2 30 percent identity / 50 percent coverage cluster. DESCRIPTIVE columns: the same E cut without the separation; all 187.
  - Cluster attributes (proposal): clade = the most frequent clade_fine of the cluster's proteins (ties: alphabetical); tier = the best tier
    (T1 before T2 before T3); subfamily = SF12 if any protein is SF12, else the most frequent subfamily; best PANTHER E = the smallest.

Subcommands
  select    --candidates --out DIR     test_clusters.tsv, test_proteins.tsv (ids and metadata), test_cluster_ids.txt, headline_pe_cluster_ids.txt,
                                       counts.json
  separate  --clusters-tsv F --train-ids F --test-ids F --out F    separation.tsv from the `_cluster.tsv` of mmseqs easy-cluster
  purity    --candidates --members --out DIR                       the 481 training members versus the PANTHER-called UniRef50 clusters
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

HEADLINE_PANTHER_E = 1e-10  # PROPOSAL (freeze: headline_panther_evalue)
NEW = "new_cluster"
TIER_ORDER = ["T1_swissprot", "T2_named", "T3_panther_only"]


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r[h] for h in header])


def select_proteins(rows: list[dict]) -> list[dict]:
    """Candidate proteins whose UniRef50 cluster is not a training member."""
    return [r for r in rows if r["overlap_with_training"] == NEW]


def _modal(values: list[str]) -> str:
    c = collections.Counter(values)
    top = max(c.values())
    return sorted(v for v, n in c.items() if n == top)[0]


def cluster_table(proteins: list[dict], pe: float = HEADLINE_PANTHER_E) -> list[dict]:
    """One row per UniRef50 cluster of the selected proteins (the unit of the evaluation)."""
    by = collections.defaultdict(list)
    for p in proteins:
        by[p["uniref50"]].append(p)
    out = []
    for cid in sorted(by):
        ps = by[cid]
        clades = [p["clade_fine"] for p in ps]
        tiers = {p["tier"] for p in ps}
        sfs = [p["panther_subfamily"] for p in ps]
        best = min(float(p["panther_evalue"]) for p in ps)
        out.append(
            {
                "cluster_id": cid,
                "clade_fine": _modal(clades),
                "clade_multi": "yes" if len(set(clades)) > 1 else "no",
                "tier": next((t for t in TIER_ORDER if t in tiers), sorted(tiers)[0]),
                "tier_multi": "yes" if len(tiers) > 1 else "no",
                "subfamily": "SF12" if "SF12" in sfs else _modal(sfs),
                "n_proteins": len(ps),
                "best_panther_evalue": f"{best:.3g}",
                "headline_pe": "yes" if best <= pe else "no",
                "n_proteins_pe": sum(1 for p in ps if float(p["panther_evalue"]) <= pe),
            }
        )
    return out


def counts(proteins: list[dict], clusters: list[dict], pe: float = HEADLINE_PANTHER_E) -> dict:
    strong = [p for p in proteins if float(p["panther_evalue"]) <= pe]
    return {
        "n_proteins_all": len(proteins),
        "n_clusters_all": len(clusters),
        "n_proteins_pe": len(strong),
        "n_clusters_pe": sum(1 for c in clusters if c["headline_pe"] == "yes"),
        "n_clusters_with_two_clades": sum(1 for c in clusters if c["clade_multi"] == "yes"),
        "n_clusters_with_two_tiers": sum(1 for c in clusters if c["tier_multi"] == "yes"),
        "tiers_proteins": dict(collections.Counter(p["tier"] for p in proteins)),
        "tiers_clusters": dict(collections.Counter(c["tier"] for c in clusters)),
        "subfamily_proteins": dict(collections.Counter(p["panther_subfamily"] for p in proteins)),
        "subfamily_clusters": dict(collections.Counter(c["subfamily"] for c in clusters)),
        "clades_clusters": dict(collections.Counter(c["clade_fine"] for c in clusters)),
        "panther_evalue_threshold": pe,
    }


CLUSTER_HEADER = [
    "cluster_id",
    "clade_fine",
    "clade_multi",
    "tier",
    "tier_multi",
    "subfamily",
    "n_proteins",
    "best_panther_evalue",
    "headline_pe",
    "n_proteins_pe",
]
PROTEIN_HEADER = [
    "accession",
    "uniref50",
    "organism",
    "clade_fine",
    "tier",
    "panther_subfamily",
    "panther_evalue",
]


def separation(cluster_of: dict[str, str], train_ids: list[str], test_ids: list[str]) -> list[dict]:
    """The cluster30 separation. `cluster_of` is {member: representative} of an MMseqs2 clustering of the training members and the test
    sequences TOGETHER. A test sequence is NOT separated when its cluster holds at least one training member. A missing id is an error
    (a sequence that was not clustered would pass silently)."""
    missing = [i for i in list(train_ids) + list(test_ids) if i not in cluster_of]
    if missing:
        raise SystemExit(
            f"FATAL: {len(missing)} ids are absent from the clustering, e.g. {missing[:3]}"
        )
    if set(train_ids) & set(test_ids):
        raise SystemExit("FATAL: an id is both a training member and a test sequence")
    train_by_cluster = collections.defaultdict(list)
    for t in train_ids:
        train_by_cluster[cluster_of[t]].append(t)
    rows = []
    for t in test_ids:
        shared = sorted(train_by_cluster.get(cluster_of[t], []))
        rows.append(
            {
                "cluster_id": t,
                "mmseqs_cluster": cluster_of[t],
                "n_training_in_cluster": len(shared),
                "training_ids": ",".join(shared),
                "separated": "no" if shared else "yes",
            }
        )
    return rows


SEP_HEADER = ["cluster_id", "mmseqs_cluster", "n_training_in_cluster", "training_ids", "separated"]


def read_cluster_tsv(path: Path) -> dict[str, str]:
    out = {}
    for line in open(path):
        if line.strip():
            rep, mem = line.rstrip("\n").split("\t")[:2]
            out[mem] = rep
    return out


def purity(candidates: list[dict], members: list[dict]) -> tuple[list[dict], dict]:
    """The 481 training members versus the PANTHER-called UniRef50 clusters (descriptive; no job). A member is `in_panther_cluster` when its
    UniRef50 cluster id (hit_id of members.tsv) is the uniref50 of any PANTHER PTHR31018 candidate protein, whatever its overlap class."""
    called = {c["uniref50"] for c in candidates}
    rows = []
    for m in members:
        rows.append(
            {
                "hit_id": m["hit_id"],
                "clade_fine": m["clade_fine"],
                "is_query": m.get("is_query", ""),
                "in_panther_cluster": "yes" if m["hit_id"] in called else "no",
            }
        )
    by_clade = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        by_clade[r["clade_fine"]][0] += 1
        by_clade[r["clade_fine"]][1] += r["in_panther_cluster"] == "no"
    summ = {
        "n_members": len(rows),
        "n_not_in_panther_cluster": sum(r["in_panther_cluster"] == "no" for r in rows),
        "by_clade_members_and_not_in_panther": dict(sorted(by_clade.items())),
        "note": "descriptive; PANTHER is itself an HMM, so a member outside every PANTHER-called cluster is a disagreement between two profile models, not proof that it is not Ecm33",
    }
    return rows, summ


def cmd_select(a) -> int:
    prot = select_proteins(read_tsv(a.candidates))
    clu = cluster_table(prot)
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "test_clusters.tsv", CLUSTER_HEADER, clu)
    write_tsv(
        a.out / "test_proteins.tsv",
        PROTEIN_HEADER,
        [{k: p[k] for k in PROTEIN_HEADER} for p in sorted(prot, key=lambda p: p["accession"])],
    )
    (a.out / "test_cluster_ids.txt").write_text("".join(c["cluster_id"] + "\n" for c in clu))
    (a.out / "headline_pe_cluster_ids.txt").write_text(
        "".join(c["cluster_id"] + "\n" for c in clu if c["headline_pe"] == "yes")
    )
    cn = counts(prot, clu)
    (a.out / "counts.json").write_text(json.dumps(cn, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: cn[k] for k in cn if k.startswith("n_")}, sort_keys=True))
    return 0


def cmd_separate(a) -> int:
    rows = separation(
        read_cluster_tsv(a.clusters_tsv),
        [x.strip() for x in open(a.train_ids) if x.strip()],
        [x.strip() for x in open(a.test_ids) if x.strip()],
    )
    write_tsv(a.out, SEP_HEADER, rows)
    print(f"{sum(r['separated'] == 'yes' for r in rows)} of {len(rows)} test sequences separated")
    return 0


def cmd_purity(a) -> int:
    rows, summ = purity(read_tsv(a.candidates), read_tsv(a.members))
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        a.out / "purity_members.tsv",
        ["hit_id", "clade_fine", "is_query", "in_panther_cluster"],
        rows,
    )
    (a.out / "purity.json").write_text(json.dumps(summ, indent=1, sort_keys=True) + "\n")
    print(
        summ["n_not_in_panther_cluster"],
        "of",
        summ["n_members"],
        "training members are in no PANTHER-called cluster",
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--candidates", required=True, type=Path)
    s.add_argument("--out", required=True, type=Path)
    s.set_defaults(fn=cmd_select)
    p = sub.add_parser("separate")
    p.add_argument("--clusters-tsv", required=True, type=Path)
    p.add_argument("--train-ids", required=True, type=Path)
    p.add_argument("--test-ids", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.set_defaults(fn=cmd_separate)
    u = sub.add_parser("purity")
    u.add_argument("--candidates", required=True, type=Path)
    u.add_argument("--members", required=True, type=Path)
    u.add_argument("--out", required=True, type=Path)
    u.set_defaults(fn=cmd_purity)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
