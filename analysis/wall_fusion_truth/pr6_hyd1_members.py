#!/usr/bin/env python3
"""PR 6 for hyd1_like (KS0017, A0ACE4KFX7, H. capsulatum, adhesion slice): the family members and the held-out units.

DRAFT, not run on HPCC. Stdlib only; run with `python3 -I`. It reads two committed tables and writes tables. It writes no call, no status, no
tier and NO SEQUENCE. Spec: docs/superpowers/specs/2026-10-23-pr6-hyd1-like-design.md. Owner decisions: spec 2026-10-22 (sections 4 to 7b) and
docs/reports/2026-10-22-pr6-drift-check.md section 4b (the rule is QUERY coverage >= 0.5; hyd1_like stays a target).

Inputs
  --drift-hits      docs/reports/data/pr6_drift_check_2026-10-22/drift_hits.tsv (one line per jackhmmer outside hit: class, query_cover, clade_fine)
  --hits-manifest   docs/reports/data/pr6_drift_check_2026-10-22/hits_manifest.tsv (taxid, evalue, score, in_phmmer of the same hits)

Members = the query itself + every drift hit of A0ACE4KFX7 with class `supported_homolog` AND query_cover >= 0.5 AND a named fine clade.
A supported hit with no fine clade (the `(no clade)` bucket) is listed in excluded.tsv: it is NOT a clade, is no held-out unit and does not count
toward the training rule. Inside-clade (Eurotiomycetes) hits are NOT members here: the drift check did not cover them (spec section 2, owner question 3).

Units. A unit is a named fine clade. `Eurotiomycetes` holds only the query, so that unit IS the Kiss unit (leave-cluster-out of KS0017: hold out
the one Kiss protein, train on the outside homologs). Status per unit:
  unit          >= MIN_UNIT_MEMBERS members and the training rule met
  kiss_single   the query's clade (1 member): a single-protein line, no recall fraction means anything, never a rule unit
  too_few_members / not_buildable
The training rule has two readings, both written to units.tsv (the freeze names the one in force; owner question 2):
  total       >= MIN_TRAIN_SEQ sequences in all outside clades together AND >= MIN_TRAIN_CLADES outside clades (what pr6_ecm33_members.py computes)
  per_clade   >= MIN_TRAIN_CLADES outside clades EACH holding >= MIN_TRAIN_SEQ sequences
Both are counted on the members, BEFORE the 30 percent independence rule (applied later by pr6_hyd1_build.py after the MMseqs2 clustering).

What a member is NOT: a curated hydrophobin or a proven hyd1_like family member. It is a jackhmmer-iteration hit of the Kiss protein that passed the drift
check; 27 of the 32 were reached only by iteration (in_phmmer = no). So the held-out truth of a clade unit is defined by the search that supplies
the training set (circular), and every cell stays `descriptive`.

Outputs (--out DIR): members.tsv, excluded.tsv, units.tsv, member_ids.txt (UniRef50 ids to fetch; the query is read from the Kiss control set), summary.txt
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

QUERY = "A0ACE4KFX7"
QUERY_CLADE = "Eurotiomycetes"  # clade_fine of the query in data/controls/wall-fusion-kiss2019/controls.tsv (a test checks)
MIN_QUERY_COVER = 0.5  # owner decision 2026-10-22 (drift report section 4b)
MIN_UNIT_MEMBERS = 5  # provisional, as for ecm33
MIN_TRAIN_SEQ = 5  # decision 7, provisional
MIN_TRAIN_CLADES = 2
NO_CLADE = "(no clade)"
MEMBER_HEADER = [
    "hit_id",
    "query_accession",
    "taxid",
    "clade_fine",
    "where",
    "score",
    "evalue",
    "is_query",
    "query_cover",
    "in_phmmer",
]


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def select_members(
    drift: list[dict], manifest: list[dict], query: str = QUERY, min_cover: float = MIN_QUERY_COVER
) -> tuple[list[dict], list[dict]]:
    """(members incl. the query, excluded supported hits without a named clade). Unsupported hits are dropped silently."""
    info = {m["hit"]: m for m in manifest if m["query"] == query}
    members = [
        {
            "hit_id": query,
            "query_accession": query,
            "taxid": "5037",
            "clade_fine": QUERY_CLADE,
            "where": "self",
            "score": "",
            "evalue": "",
            "is_query": "yes",
            "query_cover": "1.000",
            "in_phmmer": "",
        }
    ]
    excluded = []
    for h in drift:
        if h["query"] != query or h["class"] != "supported_homolog":
            continue
        if float(h["query_cover"]) < min_cover:
            continue
        m = info.get(h["hit"], {})
        clade = (h["clade_fine"] or "").strip()
        row = {
            "hit_id": h["hit"],
            "query_accession": "",
            "taxid": m.get("taxid", ""),
            "clade_fine": clade,
            "where": "outside_clade",
            "score": m.get("score", ""),
            "evalue": m.get("evalue", ""),
            "is_query": "no",
            "query_cover": h["query_cover"],
            "in_phmmer": m.get("in_phmmer", ""),
        }
        if clade and clade != NO_CLADE:
            members.append(row)
        else:
            excluded.append(row)
    members.sort(key=lambda r: (r["clade_fine"], r["hit_id"]))
    excluded.sort(key=lambda r: r["hit_id"])
    return members, excluded


def training_rules(outside: dict[str, int]) -> tuple[bool, bool]:
    """(total reading, per_clade reading) for {outside clade: n sequences}."""
    total = sum(outside.values()) >= MIN_TRAIN_SEQ and len(outside) >= MIN_TRAIN_CLADES
    per = sum(1 for n in outside.values() if n >= MIN_TRAIN_SEQ) >= MIN_TRAIN_CLADES
    return total, per


def unit_table(
    members: list[dict], reading: str = "total", min_unit: int = MIN_UNIT_MEMBERS
) -> list[dict]:
    if reading not in ("total", "per_clade"):
        raise ValueError(reading)
    by: dict[str, int] = {}
    for m in members:
        by[m["clade_fine"]] = by.get(m["clade_fine"], 0) + 1
    rows = []
    for clade in sorted(by):
        outside = {c: n for c, n in by.items() if c != clade}
        tot, per = training_rules(outside)
        met = tot if reading == "total" else per
        kiss = any(m["is_query"] == "yes" and m["clade_fine"] == clade for m in members)
        if kiss and by[clade] == 1:
            status = "kiss_single" if met else "not_buildable"
        elif not met:
            status = "not_buildable"
        elif by[clade] < min_unit:
            status = "too_few_members"
        else:
            status = "unit"
        rows.append(
            {
                "unit": f"clade:{clade}",
                "clade_fine": clade,
                "n_members": by[clade],
                "n_outside_seq": sum(outside.values()),
                "n_outside_clades": len(outside),
                "rule_total": "yes" if tot else "no",
                "rule_per_clade": "yes" if per else "no",
                "kiss_unit": "yes" if kiss else "no",
                "status": status,
            }
        )
    return rows


UNIT_HEADER = [
    "unit",
    "clade_fine",
    "n_members",
    "n_outside_seq",
    "n_outside_clades",
    "rule_total",
    "rule_per_clade",
    "kiss_unit",
    "status",
]


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
    ap.add_argument("--drift-hits", required=True, type=Path)
    ap.add_argument("--hits-manifest", required=True, type=Path)
    ap.add_argument("--reading", default="total", choices=["total", "per_clade"])
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    members, excluded = select_members(read_tsv(a.drift_hits), read_tsv(a.hits_manifest))
    if len(members) < 2:
        print(f"FATAL: no supported hit of {QUERY} in {a.drift_hits}", file=sys.stderr)
        return 1
    units = unit_table(members, a.reading)
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "members.tsv", MEMBER_HEADER, members)
    write_tsv(a.out / "excluded.tsv", MEMBER_HEADER, excluded)
    write_tsv(a.out / "units.tsv", UNIT_HEADER, units)
    ids = sorted(m["hit_id"] for m in members if m["is_query"] != "yes")
    (a.out / "member_ids.txt").write_text("".join(i + "\n" for i in ids))
    n_rule = {r["unit"]: r["status"] for r in units}
    summary = [
        f"query {QUERY} (KS0017), drift-supported hits with query cover >= {MIN_QUERY_COVER}, reading {a.reading}",
        f"members {len(members)} = query + {len(members) - 1} hits in {len({m['clade_fine'] for m in members})} named clades; excluded without a named clade {len(excluded)} (listed, not counted)",
        "units: " + "; ".join(f"{u} {s}" for u, s in sorted(n_rule.items())),
        "every cell stays descriptive: members are search hits, not independent curated truth",
    ]
    (a.out / "summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
