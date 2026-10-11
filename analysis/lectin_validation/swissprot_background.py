#!/usr/bin/env python3.12
"""Describe the hits of the 15 M1/M2 models among all fungal Swiss-Prot reviewed proteins.

This is a description, not a measure of specificity: an unlabelled hit is not a false hit. The
labels come from the truth table; every other hit is `unlabelled`.
Usage: swissprot_background.py SP_HITS_TSV SWISSPROT_TSV_GZ TRUTH_TSV OUT_PREFIX
Writes OUT_PREFIX.summary.tsv (per model: hits, labelled positive, labelled negative, unlabelled)
and OUT_PREFIX.unlabelled.tsv (the unlabelled hits with name, species and score).
"""

import csv
import gzip
import sys
from collections import defaultdict

M12 = [
    "PF07938",
    "PF01476",
    "PF08881",
    "PF07367",
    "PF21691",
    "PF00337",
    "PF21595",
    "PF18647",
    "PF22803",
    "PF01419",
    "PF00652",
    "PF14200",
    "PF13385",
    "PF01607",
    "PF00187",
]
hits_path, sp_path, truth_path, prefix = sys.argv[1:5]
truth = {r["accession"]: r for r in csv.DictReader(open(truth_path), delimiter="\t")}
sp = {r["Entry"]: r for r in csv.DictReader(gzip.open(sp_path, "rt"), delimiter="\t")}
per = defaultdict(lambda: {"hits": 0, "pos": 0, "neg": 0, "unl": 0})
unl = []
for r in csv.DictReader(open(hits_path), delimiter="\t"):
    if r["model"] not in M12:
        continue
    d = per[r["model"]]
    d["hits"] += 1
    t = truth.get(r["accession"])
    if t:
        d["pos" if t["label"] == "positive" else "neg"] += 1
    else:
        d["unl"] += 1
        s = sp.get(r["accession"], {})
        unl.append(
            {
                "model": r["model"],
                "model_name": r["model_name"],
                "accession": r["accession"],
                "score": r["score"],
                "species": s.get("Organism", "")[:40],
                "protein_name": s.get("Protein names", "")[:70],
                "uniprot_keyword_lectin": "yes" if "Lectin" in s.get("Keywords", "") else "no",
            }
        )
with open(prefix + ".summary.tsv", "w") as fh:
    fh.write("model\thits_in_swissprot_fungi\tlabelled_positive\tlabelled_negative\tunlabelled\n")
    for m in M12:
        d = per[m]
        fh.write(f"{m}\t{d['hits']}\t{d['pos']}\t{d['neg']}\t{d['unl']}\n")
unl.sort(key=lambda x: (x["model"], -float(x["score"])))
with open(prefix + ".unlabelled.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(unl[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(unl)
distinct = set()
for r in csv.DictReader(open(hits_path), delimiter="\t"):
    if r["model"] in M12:
        distinct.add(r["accession"])
with open(prefix + ".distinct.tsv", "w") as fh:
    fh.write("item\tvalue\n")
    fh.write(f"distinct_proteins_with_a_hit_to_the_15_models\t{len(distinct)}\n")
    fh.write(f"swissprot_fungi_proteins\t{len(sp)}\n")
print(len(unl), "unlabelled hits")
