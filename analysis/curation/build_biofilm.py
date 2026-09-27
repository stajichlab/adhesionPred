#!/usr/bin/env python
"""Build the draft biofilm label table (data/curated/biofilm/biofilm.tsv).

Biofilm involvement is kept separate from adhesion. Many biofilm genes are regulators or
matrix/cell-wall enzymes, not adhesins, and GO process annotations do not say whether a
protein promotes or restrains biofilm.

Source: QuickGO experimental (ECO:0000269 descendants) annotations in Fungi to
  GO:0042710 biofilm formation (is_a/part_of descendants), and
  GO:1900190 regulation of single-species biofilm formation (incl. positive/negative)

Putative classes (need expert confirmation):
  biofilm_surface    surface/secreted protein (candidate matrix, adhesin or cell-wall effector)
  biofilm_regulator  annotated to *regulation* of biofilm formation, or a transcription
                     factor / kinase / chromatin regulator by UniProt keywords
  biofilm_other      non-surface protein annotated to the biofilm process itself
Direction comes from regulation terms where available:
  promotes (positive regulation), restrains (negative regulation), unknown otherwise.

Usage:  python analysis/curation/build_biofilm.py
"""

import csv
import re

from curation_lib import (
    CURATED,
    REFERENCE_TAXA,
    apply_overrides,
    base_row,
    go_summary,
    is_surface,
    log,
    quickgo_experimental,
    uniprot_accessions,
    write_table,
)

CUR = CURATED / "biofilm"
PROCESS = {"GO:0042710": "biofilm formation"}
REGULATION = {"GO:1900190": "regulation of single-species biofilm formation"}
POSITIVE = "GO:1900192"
NEGATIVE = "GO:1900191"
# UniProt keywords marking transcription factors, signaling kinases and chromatin regulators
REGULATOR_KW = re.compile(
    r"Transcription regulation|DNA-binding|Kinase|Activator|Repressor|Chromatin regulator"
)


def main():
    log("QuickGO experimental biofilm annotations...")
    process = quickgo_experimental(PROCESS)
    regulation = quickgo_experimental(REGULATION)
    uni = uniprot_accessions(set(process) | set(regulation))

    adhesins = {}
    adh_path = CURATED / "adhesins" / "adhesins.tsv"
    if adh_path.exists():
        for r in csv.DictReader(open(adh_path), delimiter="\t"):
            if r["accession"]:
                adhesins[r["accession"]] = f"{r['cls']}:{r['evidence_level']}"

    out = {}
    for acc in sorted(set(process) | set(regulation)):
        u = uni.get(acc)
        if u is None:
            continue
        hits = process.get(acc, []) + regulation.get(acc, [])
        row = base_row(u, "go_exp")
        row.update(go_summary(hits))
        terms = {h[0] for h in hits}
        pos, neg = POSITIVE in terms, NEGATIVE in terms
        direction = (
            "mixed" if pos and neg else "promotes" if pos else "restrains" if neg else "unknown"
        )
        regulatory = bool(REGULATOR_KW.search(u["Keywords"]))
        if is_surface(u) and not regulatory:
            cls = "biofilm_surface"
        elif acc in regulation or regulatory:
            cls = "biofilm_regulator"
        else:
            cls = "biofilm_other"
        row.update(
            cls=cls,
            direction=direction,
            adhesin_table=adhesins.get(acc, ""),
            evidence_level="E1"
            if {"IDA", "IMP", "IGI", "EXP"} & set(row["evidence_codes"].split(","))
            else "E2",
            in_reference_genome="yes" if int(row["taxon_id"]) in REFERENCE_TAXA else "no",
            needs_review="yes",
        )
        out[acc] = row

    apply_overrides(out, CUR / "manual_overrides.tsv")
    cols = [
        "accession",
        "gene",
        "protein_name",
        "organism",
        "taxon_id",
        "in_reference_genome",
        "cls",
        "direction",
        "evidence_level",
        "adhesin_table",
        "evidence_codes",
        "pmids",
        "go_terms",
        "go_qualifiers",
        "source",
        "needs_review",
        "reviewed",
        "length",
        "signal_peptide",
        "gpi_anchor",
        "surface",
        "pfam",
        "evidence_summary",
    ]
    order = {"biofilm_surface": 0, "biofilm_regulator": 1, "biofilm_other": 2}
    write_table(
        out.values(),
        CUR / "biofilm.tsv",
        cols,
        lambda r: (order.get(r["cls"], 9), r["organism"], r["gene"]),
    )


if __name__ == "__main__":
    main()
