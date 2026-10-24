"""Build curated_orthologs.tsv (ids only) from committed tables plus the read-only answers cached in CACHE.

Usage: python3 -I build_table.py REPO_ROOT OUT_TSV
The extra fields of the literature-named accessions (organism, clade, UniRef50 cluster) were read from UniProt on
2026-10-24 (README.md lists the calls). Nothing here contains a sequence.
"""

import csv
import sys
from pathlib import Path

root, out = Path(sys.argv[1]), Path(sys.argv[2])
EV = root / "docs/reports/data/pr6_ecm33_independent_eval_2026-10-24/extract"
cand = {
    r["accession"]: r
    for r in csv.DictReader(
        open(root / "docs/reports/data/ecm33_independent_set_candidates/candidates.tsv"),
        delimiter="\t",
    )
}
members = {
    r["hit_id"]
    for r in csv.DictReader(
        open(root / "docs/reports/data/pr6_ecm33_build_2026-10-23/members/members.tsv"),
        delimiter="\t",
    )
}
tests = {r["cluster_id"] for r in csv.DictReader(open(EV / "test_clusters.tsv"), delimiter="\t")}
sep = {
    r["cluster_id"]: r["separated"]
    for r in csv.DictReader(open(EV / "separation.tsv"), delimiter="\t")
}

swiss = "A6ZL22 A6ZY20 B5VE42 C7GQJ1 D4AZK9 O13960 O74411 P08459 P25380 P38248 Q12355 Q4WNS8 Q5AGC4".split()
# (accession, pmid, how the accession was resolved)
lit = [
    ("P38248", "15583168", "Swiss-Prot (S. cerevisiae Ecm33)"),
    ("Q12355", "15583168", "Swiss-Prot (S. cerevisiae Pst1)"),
    ("Q5AGC4", "", "Swiss-Prot (C. albicans Ecm33; the Kiss query)"),
    ("O13960", "20032302", "Swiss-Prot (S. pombe ecm33)"),
    ("Q4WNS8", "16804168", "Swiss-Prot (A. fumigatus ecm33)"),
    (
        "J4UGA7",
        "24549768",
        "UniProt name search, B. bassiana ARSEF 2860 (locus not checked against the paper)",
    ),
    (
        "G9MJ35",
        "30840111",
        "UniProt name search, T. virens Gv29-8 (locus not checked against the paper)",
    ),
    (
        "N1RSC8",
        "35183338",
        "UniProt name search, F. oxysporum f. sp. cubense race 4 (locus not checked against the paper)",
    ),
]
extra = {
    "N1RSC8": (
        "Fusarium oxysporum f. sp. cubense (strain race 4)",
        "Sordariomycetes",
        "UniRef50_A0A428P255",
    )
}

cols = [
    "accession",
    "set",
    "organism",
    "clade_fine",
    "reviewed",
    "pmid",
    "how_resolved",
    "uniref50",
    "uniref50_is_training_member",
    "uniref50_is_independent_test_cluster",
    "separated_30pct",
    "independent",
]
rows = []
for acc in swiss:
    c = cand[acc]
    rows.append(
        [
            acc,
            "swissprot_PTHR31018",
            c["organism"],
            c["clade_fine"],
            "yes",
            "",
            "UniProt reviewed, PTHR31018",
            c["uniref50"],
        ]
    )
for acc, pmid, how in lit:
    if acc in extra:
        org, clade, ur = extra[acc]
    else:
        org, clade, ur = cand[acc]["organism"], cand[acc]["clade_fine"], cand[acc]["uniref50"]
    rows.append(
        [
            acc,
            "literature_named",
            org,
            clade,
            "yes" if how.startswith("Swiss") else "no",
            pmid,
            how,
            ur,
        ]
    )
out_rows = []
for r in rows:
    ur = r[7]
    train = ur in members
    test = ur in tests
    s = sep.get(ur, "") if test else ""
    indep = "yes" if (not train and s == "yes") else "no"
    out_rows.append(
        r
        + [
            "yes" if train else "no",
            "yes" if test else "no",
            s if test else "not_applicable",
            indep,
        ]
    )
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    w.writerows(out_rows)
print(len(out_rows), "rows;", sum(r[-1] == "yes" for r in out_rows), "independent")
