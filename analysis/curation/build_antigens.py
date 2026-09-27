#!/usr/bin/env python
"""Build the draft fungal antigen table (data/curated/antigens/antigens.tsv).

Purpose: diagnostic-serology and immunodiagnostic marker discovery, with a focus on
Coccidioides. "Antigen" here means a protein with curated evidence of recognition by a
host immune response (T cell, B cell or MHC-elution assays), which is the label an
antigenicity predictor needs.

Source: IEDB Query API (https://query-api.iedb.org) antigen_search, filtered to fungal
source organisms. IEDB curates published T cell, B cell and MHC ligand assays; each row
records how many assays of each type support the protein. UniProt adds architecture
(signal peptide, GPI, Pfam) so antigens can be cross-referenced with the surface stage.

Putative classes (need expert confirmation):
  antigen_tcell     supported by T cell assays only
  antigen_bcell     supported by B cell (antibody) assays only
  antigen_both      supported by both
Evidence level is assay-count based, not a quality judgment:
  A1  >= 5 assays of the supporting type(s)
  A2  2-4 assays
  A3  a single assay

Usage:  python analysis/curation/build_antigens.py [--taxa 5500,5501,...]
"""

import argparse
import json
import urllib.parse

from curation_lib import CURATED, base_row, get, log, uniprot_accessions, write_table

CUR = CURATED / "antigens"
IEDB = "https://query-api.iedb.org/antigen_search"
# Focus taxa: Coccidioides first, then the other fungi with meaningful IEDB coverage.
DEFAULT_TAXA = {
    "Coccidioides posadasii": "focus",
    "Coccidioides immitis": "focus",
    "Aspergillus fumigatus": "comparison",
    "Candida albicans": "comparison",
    "Cryptococcus neoformans": "comparison",
    "Histoplasma capsulatum": "comparison",
    "Paracoccidioides brasiliensis": "comparison",
    "Blastomyces dermatitidis": "comparison",
}
SELECT = ",".join(
    [
        "parent_source_antigen_iri",
        "parent_source_antigen_names",
        "parent_source_antigen_source_org_name",
        "tcell_ids",
        "bcell_ids",
        "elution_ids",
        "pubmed_ids",
        "disease_names",
        "host_organism_names",
    ]
)


def iedb_antigens(organism):
    """IEDB antigen records whose source organism name matches exactly."""
    q = urllib.parse.urlencode(
        {
            "source_organism_names": "cs.{" + f'"{organism}"' + "}",
            "select": SELECT,
            "limit": 5000,
        }
    )
    d = json.loads(get(f"{IEDB}?{q}"))
    return d if isinstance(d, list) else []


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--taxa", default=None, help="comma-separated organism names (default: built-in list)"
    )
    args = ap.parse_args()
    taxa = dict.fromkeys(args.taxa.split(","), "focus") if args.taxa else DEFAULT_TAXA

    records = []
    for organism, role in taxa.items():
        hits = iedb_antigens(organism)
        log(f"IEDB {organism}: {len(hits)} antigen records")
        for h in hits:
            records.append((organism, role, h))

    # IEDB identifies proteins by IRI; only UNIPROT:<acc> rows can be resolved to sequences.
    accs = [
        h["parent_source_antigen_iri"].split(":", 1)[1]
        for _, _, h in records
        if (h["parent_source_antigen_iri"] or "").startswith("UNIPROT:")
    ]
    uni = uniprot_accessions(accs)
    log(f"resolved {len(uni)}/{len(set(accs))} UniProt accessions")

    out = {}
    unresolved = 0
    for organism, role, h in records:
        iri = h["parent_source_antigen_iri"] or ""
        acc = iri.split(":", 1)[1] if iri.startswith("UNIPROT:") else ""
        u = uni.get(acc)
        if u is None:
            unresolved += 1
            continue
        n_t = len(h["tcell_ids"] or [])
        n_b = len(h["bcell_ids"] or [])
        n_e = len(h["elution_ids"] or [])
        total = n_t + n_b
        cls = (
            "antigen_both"
            if n_t and n_b
            else "antigen_tcell"
            if n_t
            else "antigen_bcell"
            if n_b
            else "antigen_elution_only"
        )
        row = base_row(u, "iedb")
        row.update(
            cls=cls,
            evidence_level="A1" if total >= 5 else "A2" if total >= 2 else "A3",
            n_tcell_assays=str(n_t),
            n_bcell_assays=str(n_b),
            n_elution_assays=str(n_e),
            iedb_antigen_iri=iri,
            iedb_names="; ".join(h["parent_source_antigen_names"] or [])[:120],
            source_organism=organism,
            taxon_role=role,
            diseases="; ".join(sorted(set(h["disease_names"] or [])))[:120],
            hosts="; ".join(sorted(set(h["host_organism_names"] or [])))[:80],
            pmids=";".join(sorted(set(h["pubmed_ids"] or []))),
            needs_review="yes",
            evidence_summary=f"IEDB: {n_t} T cell, {n_b} B cell, {n_e} elution assays",
        )
        out[acc] = row

    if unresolved:
        log(f"{unresolved} IEDB records skipped (no UniProt accession in IEDB IRI)")

    cols = [
        "accession",
        "gene",
        "protein_name",
        "iedb_names",
        "source_organism",
        "taxon_role",
        "organism",
        "taxon_id",
        "cls",
        "evidence_level",
        "n_tcell_assays",
        "n_bcell_assays",
        "n_elution_assays",
        "diseases",
        "hosts",
        "pmids",
        "iedb_antigen_iri",
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
    order = {"focus": 0, "comparison": 1}
    write_table(
        out.values(),
        CUR / "antigens.tsv",
        cols,
        lambda r: (
            order[r["taxon_role"]],
            -int(r["n_tcell_assays"]) - int(r["n_bcell_assays"]),
            r["accession"],
        ),
    )


if __name__ == "__main__":
    main()
