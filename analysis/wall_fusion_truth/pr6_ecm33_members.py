#!/usr/bin/env python3
"""PR 6 steps 6.1 and 6.2 for ecm33 (KS0086, Q5AGC4): the family members and the held-out units.

DRAFT, not run on HPCC. Stdlib only; run with `python3 -I`. It reads a table and writes tables. It writes no call, no status, no tier
and NO SEQUENCE (ids, clades, scores and counts only). Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md
(decisions 2, 3, 4, 5, 7; sections 4.1 to 4.4).

Input
  --homologs   gap_homologs.tsv.gz of the round 2 record (docs/reports/data/pr6_gap_uniref50_2026-10-22/): one row per hit of
               the phmmer (one pass) and jackhmmer searches of the four gap queries against UniRef50 2025_03, with the
               fine clade of the hit's taxon (`clade_fine`, empty when the taxon has none) and `where`
               (self | inside_clade | outside_clade | non_fungal | no_taxid).

Members (owner decision of 2026-10-22: the CONSERVATIVE set)
  The family members are the phmmer hits of Q5AGC4 (one pass, E 1e-5) with `where` in self, inside_clade, outside_clade. The
  jackhmmer set is NOT used: it adds tens of thousands of non-fungal hits. Non-fungal hits are never members. A fungal hit
  whose taxon has no fine clade (the `(no clade)` bucket of the gap summary) is listed in `excluded.tsv` and is NOT a clade:
  it is neither a held-out unit nor counted toward the training rule. Each member is a UniRef50 REPRESENTATIVE (one sequence per
  50 percent identity cluster).

Held-out units (decision 4: wall = leave-clade-out at the fine level)
  A unit is a fine clade with at least MIN_UNIT_MEMBERS members. The training rule (decision 7, provisional) is: at least
  MIN_TRAIN_SEQ sequences in at least MIN_TRAIN_CLADES clades OUTSIDE the held-out unit. The independence rule of spec 4.4 (no
  30 percent neighbour of the held-out members in the training alignment) is applied later, by pr6_ecm33_build.py, after the
  MMseqs2 clustering of the members; this script only counts what is available before it.

What a member is NOT
  A member is a homolog of Q5AGC4 by sequence similarity (E 1e-5, UniRef50), not an independently curated Ecm33. The held-out
  "truth" of a unit is therefore defined by the same search that supplies the training members. Recall on a unit measures whether a
  model built from other clades recovers the family in an unseen clade; it is not evidence about adhesion and not call sensitivity.
  Members of one clade are related and are not independent clusters, so no cell is `estimated`.

Outputs (--out DIR)
  members.tsv      hit_id, query_accession (Q5AGC4 for the query itself), taxid, clade_fine, where, score, evalue, is_query
  excluded.tsv     the fungal hits without a fine clade, with the reason
  units.tsv        per fine clade: n_members, n_outside_seq, n_outside_clades, training_rule_met, testable, status
  member_ids.txt   the UniRef50 ids whose sequences the extract job fetches (the query itself is read from the Kiss control set)
  summary.txt      counts
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

QUERY = "Q5AGC4"
METHOD = "phmmer"
MEMBER_WHERE = ("self", "inside_clade", "outside_clade")
MIN_UNIT_MEMBERS = 5  # provisional, like the floor of 5 clusters (Kiss decision 4); no cluster notion exists for UniRef50 members
MIN_TRAIN_SEQ = 5  # decision 7, provisional
MIN_TRAIN_CLADES = 2
NO_CLADE = "(no clade)"


def read_hits(path: Path, query: str = QUERY, method: str = METHOD) -> list[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as f:
        return [
            r
            for r in csv.DictReader(f, delimiter="\t")
            if r["query"] == query and r["method"] == method
        ]


def select_members(hits: list[dict], query: str = QUERY) -> tuple[list[dict], list[dict]]:
    """(members with a fine clade, excluded fungal hits without one). Non-fungal and no-taxid hits are dropped silently
    (they are counted in the summary by the caller from the hits)."""
    members, excluded = [], []
    for h in hits:
        if h["where"] not in MEMBER_WHERE:
            continue
        clade = (h.get("clade_fine") or "").strip()
        row = {
            "hit_id": h["hit_accession"],
            "query_accession": query if h["where"] == "self" else "",
            "taxid": h["taxid"],
            "clade_fine": clade,
            "where": h["where"],
            "score": h["score"],
            "evalue": h["evalue"],
            "is_query": "yes" if h["where"] == "self" else "no",
        }
        if clade and clade != NO_CLADE:
            members.append(row)
        else:
            excluded.append(
                {**row, "reason": "fungal hit without a fine clade (not a clade; not counted)"}
            )
    members.sort(key=lambda r: (r["clade_fine"], r["hit_id"]))
    excluded.sort(key=lambda r: r["hit_id"])
    return members, excluded


def unit_table(
    members: list[dict],
    min_unit: int = MIN_UNIT_MEMBERS,
    min_seq: int = MIN_TRAIN_SEQ,
    min_clades: int = MIN_TRAIN_CLADES,
) -> list[dict]:
    """One row per fine clade. `training_rule_met`: at least min_seq members in at least min_clades OTHER clades.
    `testable`: the clade has at least min_unit members. status: unit | too_few_members | not_buildable."""
    by: dict[str, int] = {}
    for m in members:
        by[m["clade_fine"]] = by.get(m["clade_fine"], 0) + 1
    rows = []
    for clade in sorted(by):
        outside = {c: n for c, n in by.items() if c != clade}
        n_out = sum(outside.values())
        rule = n_out >= min_seq and len(outside) >= min_clades
        testable = by[clade] >= min_unit
        if not rule:
            status = "not_buildable"
        elif not testable:
            status = "too_few_members"
        else:
            status = "unit"
        rows.append(
            {
                "clade_fine": clade,
                "n_members": by[clade],
                "n_outside_seq": n_out,
                "n_outside_clades": len(outside),
                "training_rule_met": "yes" if rule else "no",
                "testable": "yes" if testable else "no",
                "status": status,
            }
        )
    return rows


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r[h] for h in header])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--homologs", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    hits = read_hits(a.homologs)
    if not hits:
        print(f"FATAL: no {METHOD} hits of {QUERY} in {a.homologs}", file=sys.stderr)
        return 1
    members, excluded = select_members(hits)
    units = unit_table(members)
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        a.out / "members.tsv",
        [
            "hit_id",
            "query_accession",
            "taxid",
            "clade_fine",
            "where",
            "score",
            "evalue",
            "is_query",
        ],
        members,
    )
    write_tsv(
        a.out / "excluded.tsv",
        [
            "hit_id",
            "query_accession",
            "taxid",
            "clade_fine",
            "where",
            "score",
            "evalue",
            "is_query",
            "reason",
        ],
        excluded,
    )
    write_tsv(
        a.out / "units.tsv",
        [
            "clade_fine",
            "n_members",
            "n_outside_seq",
            "n_outside_clades",
            "training_rule_met",
            "testable",
            "status",
        ],
        units,
    )
    ids = sorted({m["hit_id"] for m in members if m["is_query"] != "yes"})
    (a.out / "member_ids.txt").write_text("".join(i + "\n" for i in ids))
    n_non_fungal = sum(1 for h in hits if h["where"] == "non_fungal")
    n_units = sum(1 for u in units if u["status"] == "unit")
    summary = [
        f"query {QUERY}, method {METHOD}, UniRef50 representatives",
        f"hits {len(hits)}: non-fungal {n_non_fungal} (never members)",
        f"members with a fine clade {len(members)} in {len({m['clade_fine'] for m in members})} clades; fungal hits without a fine clade {len(excluded)} (listed, not counted)",
        f"held-out units (testable and training rule met) {n_units}; clades too small to test {sum(1 for u in units if u['status'] == 'too_few_members')}; not buildable {sum(1 for u in units if u['status'] == 'not_buildable')}",
        "every cell stays descriptive: members of a clade are related, not independent clusters",
    ]
    (a.out / "summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
