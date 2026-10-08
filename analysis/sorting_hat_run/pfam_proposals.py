#!/usr/bin/env python3
"""Propose a class for every Pfam hit in the review tables (owner decision 3, 2026-10-07).

These are PROPOSALS, not reviews. Every row has status "proposed". Rules use only facts in the
review tables (annotation text, length, Cys in the domain, TMHMM helices, the R0 signal peptide
call, the draft member flag). Anything that needs outside knowledge is "needs_expert". The owner
decides each family. Nothing here changes data/sorting_hat/family_table.tsv.

Classes: known_member, uncharacterised_true_member, false_domain_hit, receptor_like, needs_expert.

Rules, in order:
 1. Draft member -> known_member (the draft list is itself unreviewed).
 2. PA14, annotation names a glucosidase, no signal peptide -> false_domain_hit.
 3. CFEM, Hydrophobin, DewD (8-Cys families): annotation names another enzyme -> needs_expert;
    signal peptide not run (NA) -> needs_expert; fewer than 8 Cys in the domain -> needs_expert;
    8 Cys or more and 2+ helices, or 1 helix without a signal peptide -> receptor_like (low confidence);
    1 helix with a signal peptide -> needs_expert (the helix may be the signal anchor);
    8 Cys or more, no helix, signal peptide called -> uncharacterised_true_member (medium).
 4. Allergen_Asp_f_4, annotation says "allergen" -> uncharacterised_true_member (low: by similarity).
 5. All other hits -> needs_expert.

Usage: pfam_proposals.py --review-dir DIR --out-tsv FILE --out-md FILE
"""

import argparse
import collections
import csv
import re
import sys

CYS_FAMILIES = ("PF05730", "PF01185", "PF28987", "PF06766", "PF22354")
ENZYME = re.compile(
    r"chitinase|glucosidase|protease|peptidase|lipase|kinase|reductase|dehydrogenase", re.I
)
PROTEOMES = [
    "Afum_Af293_UniProt",
    "Afum_A1163",
    "Afum_W72310",
    "Scer_S288C",
    "Calb_SC5314",
    "Cimm_RS",
    "Cneo_H99",
]


def propose(r):
    acc = r["family"].split()[0]
    annot = r["annotation"]
    sp, cys, tm = r["signal_peptide"], int(r["domain_cys"]), r["n_tm"]
    tm_n = 0 if tm in ("NA", "") else int(tm)
    if r["draft_member"] == "yes":
        return "known_member", "medium", "in the draft member list (unreviewed)"
    if acc == "PF07691" and re.search("glucosidase", annot, re.I) and sp != "called":
        why = (
            "annotated beta-glucosidase; no signal peptide"
            if sp == "not_called"
            else "annotated beta-glucosidase; signal peptide not run"
        )
        return "false_domain_hit", "medium", why
    if acc in CYS_FAMILIES:
        if ENZYME.search(annot):
            return (
                "needs_expert",
                "low",
                f"annotation names an enzyme ({ENZYME.search(annot).group(0)}); {cys} Cys",
            )
        if sp == "NA":
            return "needs_expert", "low", f"signal peptide not run; {cys} Cys, {tm_n} helices"
        if cys < 8:
            return "needs_expert", "low", f"only {cys} Cys in the domain"
        if tm_n == 1 and sp == "called":
            return (
                "needs_expert",
                "low",
                f"{cys} Cys, 1 TMHMM helix with a signal peptide (the helix may be the signal anchor)",
            )
        if tm_n > 0:
            return "receptor_like", "low", f"{cys} Cys, {tm_n} TMHMM helices"
        if sp == "called":
            return "uncharacterised_true_member", "medium", f"{cys} Cys, signal peptide, no helix"
        return "needs_expert", "low", f"{cys} Cys, no signal peptide called"
    if acc == "PF25312" and re.search("allergen", annot, re.I):
        return (
            "uncharacterised_true_member",
            "low",
            "annotation says allergen-like (annotation by similarity, not IgE)",
        )
    return "needs_expert", "low", "no independent evidence in the table"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--review-dir", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-md", required=True)
    a = ap.parse_args()
    rows = []
    for p in PROTEOMES:
        with open(f"{a.review_dir}/{p}/review_table.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cls, conf, why = propose(r)
                rows.append(
                    {
                        "proteome": p,
                        "family": r["family"],
                        "protein": r["protein"],
                        "annotation": r["annotation"][:100],
                        "length": r["length"],
                        "domain_cys": r["domain_cys"],
                        "n_tm": r["n_tm"],
                        "signal_peptide": r["signal_peptide"],
                        "draft_member": r["draft_member"],
                        "proposed_class": cls,
                        "confidence": conf,
                        "basis": why,
                        "status": "proposed",
                    }
                )
    with open(a.out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        by[r["family"]][r["proposed_class"]] += 1
    with open(a.out_md, "w") as f:
        f.write(
            "| family | hits | known_member | uncharacterised_true_member | receptor_like | false_domain_hit | needs_expert |\n|---|---|---|---|---|---|---|\n"
        )
        for fam in sorted(by):
            c = by[fam]
            f.write(
                f"| {fam} | {sum(c.values())} | {c['known_member']} | {c['uncharacterised_true_member']} | {c['receptor_like']} | {c['false_domain_hit']} | {c['needs_expert']} |\n"
            )
    print(
        f"{len(rows)} hits; classes: {dict(collections.Counter(r['proposed_class'] for r in rows))}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
