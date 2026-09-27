#!/usr/bin/env python
"""Build the stage-1 surface glycoprotein table (data/curated/surface/surface.tsv).

This is the table that answers the central question of the review: **which fungal
cell-surface glycoproteins are adhesins and which are not?** The shipped classifier
detects the whole surface class at ~12% adhesin precision (see docs/model-review/),
so stage 2 must be trained *within* this population, using the non-adhesive members
as the discriminative negatives.

Population: every UniProt entry of a reference genome carrying a surface/secretion
keyword (Signal, GPI-anchor, Cell wall, Secreted). Each row is then cross-linked to
the other curated tables in data/curated/.

adhesion_status:
  adhesin             in adhesins.tsv as cls=adhesin (evidence level carried over)
  non_adhesin         in adhesins.tsv as hard_negative, or a surface protein with a
                      characterized non-adhesive function (surface_other_*)
  non_adhesin_putative  proposed automatically: an annotated catalytic activity (EC number
                      or enzyme keyword) accounts for the protein and no adhesion evidence
                      exists. Usable as training negatives, but flagged as proposed, since
                      a few adhesins are enzymatically active (e.g. Paracoccidioides gp43)
  unknown             no adhesion evidence either way (the majority; treat as unlabeled,
                      NOT as negative, in positive-unlabeled training)

Other columns mark biofilm involvement and antigen evidence so one protein can carry
several labels (a protein may be an adhesin *and* an antigen).

Usage:  python analysis/curation/build_surface.py
"""

import csv
import re

from curation_lib import (
    CURATED,
    REFERENCE_TAXA,
    UNIPROT_FIELDS,
    base_row,
    log,
    uniprot_accessions,
    uniprot_stream,
    write_table,
)

CUR = CURATED / "surface"
# Reference genomes plus the Coccidioides genomes for the antigen work.
GENOMES = dict(REFERENCE_TAXA)
GENOMES.update(
    {
        246410: "Coccidioides immitis RS",
        443226: "Coccidioides posadasii C735 delta SOWgp",
    }
)
SURFACE_KEYWORDS = [
    "KW-0732",
    "KW-0336",
    "KW-0134",
    "KW-0964",
]  # Signal, GPI-anchor, Cell wall, Secreted
FIELDS = UNIPROT_FIELDS + ",ec,protein_families"
# A surface protein with a characterized catalytic activity is a strong non-adhesin
# candidate: its function is accounted for, and enzymes are what the shipped model
# most often mistakes for adhesins (AA1 laccases, glucanases, proteases).
ENZYME_KW = re.compile(
    r"Hydrolase|Transferase|Oxidoreductase|Lyase|Isomerase|Ligase|Protease|Glycosidase"
)


def load_table(path, key="accession"):
    if not path.exists():
        log(f"  missing (skipped): {path.name}")
        return {}
    return {r[key]: r for r in csv.DictReader(open(path), delimiter="\t") if r.get(key)}


def main():
    adhesins = load_table(CURATED / "adhesins" / "adhesins.tsv")
    biofilm = load_table(CURATED / "biofilm" / "biofilm.tsv")
    antigens = load_table(CURATED / "antigens" / "antigens.tsv")

    out = {}
    kw = " OR ".join(f"keyword:{k}" for k in SURFACE_KEYWORDS)
    for taxon, name in GENOMES.items():
        hits = uniprot_stream(f"organism_id:{taxon} AND ({kw})", fields=FIELDS)
        log(f"{name}: {len(hits)} surface-keyword proteins")
        for u in hits:
            row = base_row(u, "uniprot_surface_kw")
            a, b, g = adhesins.get(u["Entry"]), biofilm.get(u["Entry"]), antigens.get(u["Entry"])

            if a and a["cls"] == "adhesin" and a.get("moonlighting", "").upper() == "YES":
                status, level = "moonlighting_excluded", a["evidence_level"]
            elif a and a["cls"] == "adhesin" and a["evidence_level"] == "E3":
                status, level = "unknown", ""
            elif a and a["cls"] == "adhesin":
                status, level = "adhesin", a["evidence_level"]
            elif a and a["cls"] == "hard_negative":
                status, level = "non_adhesin", a["evidence_level"]
            elif a and a["cls"] == "surface_other_adhesion_phenotype":
                status, level = "non_adhesin", "S1"
            elif u.get("EC number") or ENZYME_KW.search(u["Keywords"]):
                # Proposed, not curated: catalytic surface protein with no adhesion evidence.
                status, level = "non_adhesin_putative", "N3"
            else:
                status, level = "unknown", ""

            row.update(
                genome=name,
                adhesion_status=status,
                adhesion_evidence_level=level,
                adhesion_note=(
                    a["evidence_summary"][:100]
                    if a
                    else (
                        "catalytic activity annotated: " + (u.get("EC number") or "enzyme keyword")
                        if status == "non_adhesin_putative"
                        else ""
                    )
                ),
                ec=u.get("EC number", ""),
                protein_family=u.get("Protein families", "")[:60],
                biofilm_class=b["cls"] if b else "",
                biofilm_direction=b["direction"] if b else "",
                antigen_class=g["cls"] if g else "",
                antigen_assays=(f"T{g['n_tcell_assays']}/B{g['n_bcell_assays']}" if g else ""),
                labels=",".join(
                    filter(
                        None,
                        [
                            "adhesin" if status == "adhesin" else "",
                            "non_adhesin" if status == "non_adhesin" else "",
                            "non_adhesin_putative" if status == "non_adhesin_putative" else "",
                            "biofilm" if b else "",
                            "antigen" if g else "",
                        ],
                    )
                ),
            )
            out[u["Entry"]] = row

    # Curated labels whose organism is not one of GENOMES above -- SOWgp from a Coccidioides
    # isolate other than the reference strain, BAD1 from Blastomyces, Mp1p from Talaromyces --
    # would otherwise be dropped purely because this table is keyed on reference-proteome
    # taxon ids. They are the scarcest labels we have, so they are added explicitly.
    extra = [
        a for a, r in adhesins.items() if a not in out and r["cls"] in ("adhesin", "hard_negative")
    ]
    if extra:
        fetched = uniprot_accessions(extra, fields=FIELDS)
        for acc, u in fetched.items():
            a = adhesins[acc]
            b, gg = biofilm.get(acc), antigens.get(acc)
            if a.get("moonlighting", "").upper() == "YES":
                # Surface-localized cytoplasmic protein (e.g. Hsp60). Real host-binding
                # evidence, but no secretion signal, so it cannot be a stage-2 positive in a
                # signal-peptide-defined population. Held out explicitly.
                status = "moonlighting_excluded"
            elif a["evidence_level"] == "E3":
                # Domain-only guess; not strong enough to train on.
                status = "unknown"
            else:
                status = "adhesin" if a["cls"] == "adhesin" else "non_adhesin"
            row = base_row(u, "curated_literature")
            row.update(
                genome=f"{u['Organism'][:44]} [curated, non-reference]",
                adhesion_status=status,
                adhesion_evidence_level=a["evidence_level"],
                adhesion_note=a["evidence_summary"][:100],
                ec=u.get("EC number", ""),
                protein_family=a.get("family", "")[:60],
                biofilm_class=b["cls"] if b else "",
                biofilm_direction=b["direction"] if b else "",
                antigen_class=gg["cls"] if gg else "",
                antigen_assays=(f"T{gg['n_tcell_assays']}/B{gg['n_bcell_assays']}" if gg else ""),
                labels=",".join(
                    filter(None, [status, "biofilm" if b else "", "antigen" if gg else ""])
                ),
            )
            out[acc] = row
        log(f"added {len(fetched)} curated proteins from outside the reference genomes")

    cols = [
        "accession",
        "gene",
        "protein_name",
        "genome",
        "taxon_id",
        "adhesion_status",
        "adhesion_evidence_level",
        "biofilm_class",
        "biofilm_direction",
        "antigen_class",
        "antigen_assays",
        "labels",
        "length",
        "signal_peptide",
        "gpi_anchor",
        "pfam",
        "ec",
        "moonlighting",
        "protein_family",
        "reviewed",
        "adhesion_note",
        "source",
    ]
    order = {
        "adhesin": 0,
        "non_adhesin": 1,
        "non_adhesin_putative": 2,
        "moonlighting_excluded": 3,
        "unknown": 4,
    }
    write_table(
        out.values(),
        CUR / "surface.tsv",
        cols,
        lambda r: (r["genome"], order[r["adhesion_status"]], r["gene"] or "zzz", r["accession"]),
    )


if __name__ == "__main__":
    main()
