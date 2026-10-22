#!/usr/bin/env python3
"""Clusters per held-out clade and stratum, from the MMseqs2 cluster table of Task 5.3 (stdlib only).

Reads `clusters_split.tsv` (written by seed_split.py). Only clusters that hold a positive control
count; the 35 family representatives are not clustered and never appear here. For each slice and clade
(fine and coarse level) it counts the clusters present in the clade, split into

  clade_only  every positive member of the cluster is in this clade, and
  shared      the cluster also has positive members in another clade.

Status (owner decisions of 2026-10-22, report section 10): `estimated` needs at least FLOOR (5)
clusters in that count; fewer is `descriptive`; zero is `none`. The fusion slice is always
`descriptive`. The floor counts all clusters present (`status_any`). No performance number is made here.

Usage:  python3 -I analysis/wall_fusion_truth/clade_floors.py CLUSTERS_SPLIT_TSV [OUT_TSV]
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict

FLOOR = 5
LEVELS = (("fine", "clades_fine"), ("coarse", "clades_coarse"))
HEADER = [
    "level",
    "slice",
    "clade",
    "clusters_any",
    "clusters_clade_only",
    "clusters_shared",
    "status_any",
    "status_clade_only",
    "status_shared",
]


def status(slice_: str, n: int) -> str:
    if n == 0:
        return "none"
    return "estimated" if n >= FLOOR and slice_ != "fusion" else "descriptive"


def clade_floors(clusters: list[dict]) -> list[list]:
    pos = [r for r in clusters if "positive" in r["classes"].split(",")]
    out = []
    for level, col in LEVELS:
        t: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0, 0])
        for r in pos:
            clades = [x for x in r[col].split(",") if x]
            for x in clades:
                k = (r["slices"], x)
                t[k][0] += 1
                t[k][1 if len(clades) == 1 else 2] += 1
        for (sl, x), (a, o, s) in sorted(t.items()):
            out.append([level, sl, x, a, o, s, status(sl, a), status(sl, o), status(sl, s)])
    return out


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
    with open(argv[0], newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    text = (
        "\t".join(HEADER)
        + "\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in clade_floors(rows))
    )
    if len(argv) > 1:
        with open(argv[1], "w") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
