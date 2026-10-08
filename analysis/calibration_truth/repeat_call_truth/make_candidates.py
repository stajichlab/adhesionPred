#!/usr/bin/env python3
"""Label the curated proteins of two species for the repeat call (input of build_truth.py).

Inputs: curated_rows_mapped.tsv (curated rows joined to run protein IDs), uniprot_repeat_features.tsv
(analysis/calibration_truth/uniprot_repeat_features.py), and the adjudicated consensus table.
Output: repeat_truth_candidates.tsv. Rules are in build_truth.py.
"""

import csv
import sys

D = sys.argv[1] if len(sys.argv) > 1 else "."
CONS = (
    sys.argv[2]
    if len(sys.argv) > 2
    else "../../../data/controls/repeat-mechanism/adjudication/curation_table.consensus.tsv"
)
rows = list(csv.DictReader(open(f"{D}/curated_rows_mapped.tsv"), delimiter="\t"))
uni = {
    r["accession"]: r
    for r in csv.DictReader(open(f"{D}/uniprot_repeat_features.tsv"), delimiter="\t")
}
cons = {r["accession"]: r for r in csv.DictReader(open(CONS), delimiter="\t")}
out = []
for r in rows:
    if r["proteome"] not in ("Scer_S288C", "Calb_SC5314"):
        continue
    u, c = uni.get(r["accession"]), cons.get(r["accession"])
    n = int(u["n_repeat_features"]) if u else None
    paper = bool(c and c["final_label"] == "2a" and c["repeat_evidence_tier"] == "stated_in_paper")
    famonly = bool(
        c and c["final_label"] == "2a" and c["repeat_evidence_tier"] == "family_inference"
    )
    if paper or (n is not None and n >= 2):
        lab, why = 1, ("paper statement" if paper else "UniProt >=2 repeat features")
    elif famonly:
        lab, why = "", "family inference only (excluded)"
    elif n == 1:
        lab, why = "", "1 UniProt repeat feature (ambiguous, excluded)"
    elif n == 0:
        lab, why = 0, "no UniProt repeat feature (assumed negative)"
    else:
        lab, why = "", "no UniProt record"
    out.append(
        {
            "proteome": r["proteome"],
            "protein": r["protein"],
            "gene": r["gene"],
            "accession": r["accession"],
            "curated_class": r["cls"],
            "evidence": r["evidence"],
            "label": lab,
            "basis": why,
            "tuned_or_homolog": r["tuned_or_homolog"],
            "repeat_call": r["tandem_repeat_protein"],
        }
    )
w = csv.DictWriter(
    open(f"{D}/repeat_truth_candidates.tsv", "w"),
    fieldnames=list(out[0]),
    delimiter="\t",
    lineterminator="\n",
)
w.writeheader()
w.writerows(out)
print(len(out), "rows", file=sys.stderr)
