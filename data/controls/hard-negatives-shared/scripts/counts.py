# ruff: noqa
"""Fix gene names and write counts tables (markdown fragments) from controls.tsv. Usage: counts.py OUTDIR"""

import collections
import csv
import sys

OUT = sys.argv[1]
GENE = {"A0A0D1DYI3": "rsp3 (UMAG_03274)", "XP_045276491.1": "ENG2", "B6QUR1": "MP1 (PMAA_009820)"}
rows = list(csv.DictReader(open(f"{OUT}/controls.tsv"), delimiter="\t"))
cols = list(rows[0].keys())
for r in rows:
    r["gene"] = GENE.get(r["accession"], r["gene"])
with open(f"{OUT}/controls.tsv", "w") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
seqs = {l[1:].strip() for l in open(f"{OUT}/sequences.faa") if l.startswith(">")}


def table(key, title, rows_):
    g = collections.defaultdict(list)
    for r in rows_:
        g[key(r)].append(r)
    lines = [
        f"### {title}",
        "",
        "| Group | Rows | Sequences | Clusters | N1 rows | N2 rows | Clusters to 20 |",
        "|---|---|---|---|---|---|---|",
    ]
    for k in sorted(g):
        rs = g[k]
        ncl = len({r["cluster_id"] for r in rs})
        lines.append(
            f"| {k} | {len(rs)} | {sum(r['accession'] in seqs for r in rs)} | {ncl} | "
            f"{sum(r['evidence_level'] == 'N1' for r in rs)} | {sum(r['evidence_level'] == 'N2' for r in rs)} | {max(0, 20 - ncl)} |"
        )
    rs = rows_
    lines.append(
        f"| **all** | {len(rs)} | {sum(r['accession'] in seqs for r in rs)} | {len({r['cluster_id'] for r in rs})} | "
        f"{sum(r['evidence_level'] == 'N1' for r in rs)} | {sum(r['evidence_level'] == 'N2' for r in rs)} | |"
    )
    return "\n".join(lines) + "\n"


def cross(rows_):
    strata = sorted({r["stratum"] for r in rows_})
    clades = sorted({r["clade"] for r in rows_})
    lines = [
        "### Clusters per stratum and clade (rows in brackets)",
        "",
        "| Stratum | " + " | ".join(clades) + " | All |",
        "|---" * (len(clades) + 2) + "|",
    ]
    for s in strata:
        cells = []
        for c in clades:
            rs = [r for r in rows_ if r["stratum"] == s and r["clade"] == c]
            cells.append(f"{len({r['cluster_id'] for r in rs})} ({len(rs)})" if rs else "0")
        rs = [r for r in rows_ if r["stratum"] == s]
        cells.append(f"{len({r['cluster_id'] for r in rs})} ({len(rs)})")
        lines.append(f"| {s} | " + " | ".join(cells) + " |")
    cells = []
    for c in clades:
        rs = [r for r in rows_ if r["clade"] == c]
        cells.append(f"{len({r['cluster_id'] for r in rs})} ({len(rs)})")
    cells.append(f"{len({r['cluster_id'] for r in rows_})} ({len(rows_)})")
    lines.append("| **All** | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def serves(rows_):
    lines = [
        "### Rows and clusters per module served",
        "",
        "| serves | Rows | Clusters | Leakage none | Leakage other |",
        "|---|---|---|---|---|",
    ]
    for m in ("step1", "repeat", "family_domain", "adhesion_level"):
        rs = [r for r in rows_ if m in r["serves"].split(",")]
        lk = (
            "leakage_"
            + {
                "step1": "step1",
                "repeat": "repeat",
                "family_domain": "family_domain",
                "adhesion_level": "step1",
            }[m]
        )
        lines.append(
            f"| {m} | {len(rs)} | {len({r['cluster_id'] for r in rs})} | {sum(r[lk] == 'none' for r in rs)} | "
            f"{', '.join(f'{k}={v}' for k, v in collections.Counter(r[lk] for r in rs if r[lk] != 'none').items()) or '0'} |"
        )
    return "\n".join(lines) + "\n"


txt = [
    table(lambda r: r["stratum"], "Per stratum", rows),
    cross(rows),
    table(lambda r: r["clade"], "Per clade", rows),
    table(lambda r: f"{r['species']} ({r['taxon_id']})", "Per species", rows),
    serves(rows),
    table(
        lambda r: r["stratum"],
        "Per stratum, N1 rows only",
        [r for r in rows if r["evidence_level"] == "N1"],
    ),
    table(
        lambda r: r["stratum"],
        "Per stratum, without rows that share a cluster with a positive",
        [r for r in rows if r["mixed_cluster_with_positive"] == "no"],
    ),
]
open(f"{OUT}/_counts_fragment.md", "w").write("\n".join(txt))
print("\n".join(txt))
