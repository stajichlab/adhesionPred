#!/usr/bin/env python3.12
"""Compare an MMseqs2 cluster table with the cluster30 column of unilectin_fungal_sequences.tsv.

Usage: cluster_check.py MMSEQS_CLUSTER_TSV UNILECTIN_TABLE_TSV OUT_CLUSTERS_TSV [MODEL_SET ...]
Writes accession<TAB>cluster_id (representative accession) and prints whether the partitions are equal,
and the coverage count (sequences and clusters with a Pfam GA hit to a model set) after the three
exclusions of spec section 4.4 (PA14 class, yeast Emp L-type class, Ralstonia).
Default model sets: "five" (PF07938, PF21691, PF00337, PF07367, PF21595; revision-1 M1) and "PF07938"
(revision-2 M1). Extra sets: pass each as a comma-separated list of Pfam accessions, for example
PF07938,PF21691. The Ascomycota-only counts need the species list of review2/cov.py and are not printed here.
"""

import csv
import sys

clu_tsv, table_tsv, out = sys.argv[1:4]
SETS = {
    "five": {"PF07938", "PF21691", "PF00337", "PF07367", "PF21595"},
    "PF07938": {"PF07938"},
}
for extra in sys.argv[4:]:
    SETS[extra] = set(extra.split(","))
rep = {}
for line in open(clu_tsv):
    r, m = line.rstrip("\n").split("\t")
    rep[m] = r
rows = list(csv.DictReader(open(table_tsv), delimiter="\t"))
tab = {r["uniprot"]: r for r in rows if r["uniprot"] in rep}
with open(out, "w") as fh:
    for a in sorted(rep):
        fh.write(f"{a}\t{rep[a]}\n")
# partition equality
a_groups, b_groups = {}, {}
for a, r in tab.items():
    a_groups.setdefault(rep[a], set()).add(a)
    b_groups.setdefault(r["cluster30"], set()).add(a)
same = sorted(map(sorted, a_groups.values())) == sorted(map(sorted, b_groups.values()))
print(f"sequences {len(rep)}, clusters {len(a_groups)}, same partition as cluster30 column: {same}")


def excl(r):
    return (
        "PA14" in r["unilectin_class"]
        or "Emp" in r["unilectin_class"]
        or "Ralstonia" in r["species"]
    )


keep = {a: r for a, r in tab.items() if not excl(r)}
kc = {rep[a] for a in keep}
print(f"after exclusions: sequences {len(keep)}, clusters {len(kc)}")
for name, models in SETS.items():
    hit = {
        a
        for a, r in keep.items()
        if models & {h.split(":")[0] for h in r["pfam_ga_hits"].split(";") if h}
    }
    hc = {rep[a] for a in hit}
    print(
        f"model set {name} ({','.join(sorted(models))}): hit sequences {len(hit)}, clusters {len(hc)}"
    )
