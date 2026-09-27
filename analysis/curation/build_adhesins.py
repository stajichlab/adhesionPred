#!/usr/bin/env python
"""Build the draft curated adhesin / hard-negative table (data/curated/adhesins/adhesins.tsv).

Sources, each recorded per row in `source`:
  go_exp      QuickGO experimental (ECO:0000269 descendants) annotations to fungal
              adhesion terms: cell adhesion, adhesion of symbiont to host, flocculation,
              biofilm cell adhesion, agglutination, cell-cell adhesion (with descendants)
  literature  data/curated/adhesins/literature_seeds.tsv (adhesins without GO experimental terms)
  hard_neg    data/curated/adhesins/hard_negative_seeds.tsv
  (then)      data/curated/adhesins/manual_overrides.tsv reclassifies specific accessions
  pfam_family reference-genome proteins that carry an adhesin-defining Pfam domain seen in
              an E1 adhesin (E2 candidates by family membership)

Putative classes (all rows need expert confirmation before training):
  adhesin                          E1 direct experimental adhesion/binding evidence,
                                   E2 family/domain membership or weaker evidence
  surface_other_adhesion_phenotype surface/secreted enzyme or structural protein whose
                                   mutant affects adhesion (likely indirect)
  indirect_regulator               non-surface protein with adhesion phenotype (excluded)
  hard_negative                    N1 characterized non-adhesive function,
                                   N2 no adhesion evidence but adhesin-like architecture

Usage (network access to ebi.ac.uk and rest.uniprot.org required):
    python analysis/curation/build_adhesins.py
"""

import re
from collections import defaultdict

from curation_lib import (
    CURATED,
    REFERENCE_TAXA,
    apply_overrides,
    base_row,
    go_summary,
    is_surface,
    log,
    quickgo_experimental,
    read_seeds,
    uniprot_accessions,
    uniprot_search,
    write_table,
)

CUR = CURATED / "adhesins"

ADHESION_GO = {
    "GO:0007155": "cell adhesion",
    "GO:0044406": "adhesion of symbiont to host",
    "GO:0000128": "flocculation",
    "GO:0043709": "cell adhesion involved in single-species biofilm formation",
    "GO:0000752": "agglutination involved in conjugation with cellular fusion",
    "GO:0098609": "cell-cell adhesion",
}
# Surface/secreted proteins with adhesion phenotypes whose primary role is enzymatic,
# structural or nutritional (reclassified from go_exp; still flagged for review).
INDIRECT_SURFACE = re.compile(
    r"^(SAP\d+|BGL2|XOG1|ENG1|GAM1|GCA2|CHT3|MNT\d|PHR1|MP65|SUN41|PGA7|PGA10|RBT5|CSA2|"
    r"CWH41|ROT2|CHS7|YPS\d+|TIR\d|CCW12|UTR2|ECM33|PRA1|DFI1|SUR7|PBR1|HXK1|DSE1|BIG1|"
    r"CHK1|IRE1|ACE2|ADA2|SSR1)$",
    re.I,
)
# Pfam domains that define adhesin families (checked against E1 members below).
# PA14 (PF07691) is also found in non-adhesins, so PA14-only hits stay E2 + review.
ADHESIN_PFAM = {
    "PF00624": "Flocculin repeat",
    "PF07691": "PA14",
    "PF10182": "Flo11",
    "PF05792": "Candida_ALS",
    "PF11766": "Candida_ALS_N",
    "PF11765": "Hyr1",
    "PF09770": "PAT1-like / S. pombe adhesin (check)",
}


def main():
    out = {}

    log("QuickGO experimental adhesion annotations...")
    ann = quickgo_experimental(ADHESION_GO)
    uni = uniprot_accessions(ann)
    for acc, hits in ann.items():
        u = uni.get(acc)
        if u is None:
            continue
        row = base_row(u, "go_exp")
        row.update(go_summary(hits))
        codes = row["evidence_codes"].split(",")
        gene = row["gene"]
        if not is_surface(u):
            row.update(
                cls="indirect_regulator",
                evidence_level="-",
                family="",
                evidence_summary="adhesion phenotype of a non-surface protein",
            )
        elif INDIRECT_SURFACE.match(gene or ""):
            row.update(
                cls="surface_other_adhesion_phenotype",
                evidence_level="-",
                family="",
                evidence_summary="surface enzyme/structural protein with adhesion phenotype",
            )
        else:
            strong = {"IDA", "IMP", "IGI", "EXP"} & set(codes)
            row.update(
                cls="adhesin",
                evidence_level="E1" if strong else "E2",
                family="",
                evidence_summary="experimental GO adhesion annotation (" + ",".join(codes) + ")",
            )
        out[acc] = row

    for path, source in [
        (CUR / "literature_seeds.tsv", "literature"),
        (CUR / "hard_negative_seeds.tsv", "hard_neg"),
        (CUR / "eurotiomycetes_seeds.tsv", "eurotiomycetes"),
    ]:
        log(f"Resolving {path.name}...")
        for s in read_seeds(path):
            hits = uniprot_search(s["uniprot_query"]) if s["uniprot_query"] else []
            if not hits:
                log(f"  unresolved (kept without accession): {s['gene']}")
                out[f"UNRESOLVED:{s['gene']}"] = {
                    "accession": "",
                    "gene": s["gene"],
                    "cls": s["class"],
                    "evidence_level": s["evidence_level"],
                    "family": s["family"],
                    "evidence_summary": s["evidence_summary"],
                    "pmids": s["pmids"],
                    "evidence_codes": "literature",
                    "source": source,
                    "reviewed": "",
                    "organism": "",
                    "taxon_id": "0",
                    "pfam": "",
                    "uniprot_query": s["uniprot_query"],
                }
                continue
            u = next((h for h in hits if h["Reviewed"] == "reviewed"), hits[0])
            row = out.get(u["Entry"]) or base_row(u, source)
            if u["Entry"] in out:
                row["source"] += f",{source}"
            row.update(
                cls=s["class"],
                evidence_level=s["evidence_level"],
                family=s["family"],
                evidence_summary=s["evidence_summary"],
                pmids=";".join(
                    sorted(set(filter(None, (row.get("pmids", "") + ";" + s["pmids"]).split(";"))))
                ),
                uniprot_query=s["uniprot_query"],
                n_query_hits=str(len(hits)),
                moonlighting=s.get("moonlighting", ""),
            )
            row.setdefault("evidence_codes", "literature")
            out[u["Entry"]] = row

    log("Pfam family members in reference genomes...")
    e1_pfams = defaultdict(set)
    for r in out.values():
        if r["cls"] == "adhesin" and r["evidence_level"] == "E1":
            for p in r["pfam"].split(";"):
                if p in ADHESIN_PFAM:
                    e1_pfams[p].add(r["gene"] or r["accession"])
    for pf, members in sorted(e1_pfams.items()):
        taxa = " OR ".join(f"organism_id:{t}" for t in REFERENCE_TAXA)
        for u in uniprot_search(f"xref:pfam-{pf} AND ({taxa})", size=500):
            if u["Entry"] in out:
                continue
            row = base_row(u, "pfam_family")
            row.update(
                cls="adhesin",
                evidence_level="E2",
                family=ADHESIN_PFAM[pf],
                evidence_codes="domain",
                pmids="",
                evidence_summary=f"shares {pf} ({ADHESIN_PFAM[pf]}) with E1 adhesins "
                f"{', '.join(sorted(members)[:5])}",
            )
            out[u["Entry"]] = row

    apply_overrides(out, CUR / "manual_overrides.tsv")

    for r in out.values():
        r["needs_review"] = (
            "no"
            if (
                r["cls"] == "adhesin"
                and r["evidence_level"] == "E1"
                and r["reviewed"] == "reviewed"
            )
            else "yes"
        )
        if (
            r["cls"] == "adhesin"
            and r["evidence_level"] == "E2"
            and "PF07691" in r["pfam"]
            and not ({"PF00624", "PF10182"} & set(r["pfam"].split(";")))
        ):
            r["evidence_summary"] += "; PA14-only (PA14 also occurs in non-adhesins)"
        r["in_reference_genome"] = "yes" if int(r["taxon_id"]) in REFERENCE_TAXA else "no"

    cols = [
        "accession",
        "gene",
        "protein_name",
        "organism",
        "taxon_id",
        "in_reference_genome",
        "cls",
        "evidence_level",
        "family",
        "evidence_summary",
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
        "moonlighting",
        "uniprot_query",
        "n_query_hits",
    ]
    order = {
        "adhesin": 0,
        "surface_other_adhesion_phenotype": 1,
        "hard_negative": 2,
        "indirect_regulator": 3,
    }
    write_table(
        out.values(),
        CUR / "adhesins.tsv",
        cols,
        lambda r: (order[r["cls"]], r["evidence_level"], r["organism"], r["gene"]),
    )


if __name__ == "__main__":
    main()
