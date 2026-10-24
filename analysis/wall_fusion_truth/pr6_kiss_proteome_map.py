"""Map Kiss positives to protein ids of the four kept proteomes (pure functions, stdlib, deterministic).

Why: a proteome hit is keyed by the FASTA id of its proteome, and only Afum_Af293_UniProt headers carry a UniProt
accession (`sp|ACC|` or `tr|ACC|`). Calb_SC5314 uses CGD Assembly 22 ids (`C1_13450W_A`), Scer_S288C uses systematic
names (`YNL238W`) and Cimm_RS uses RefSeq ids (`XP_001239171.1`). A Kiss positive that sits inside a proteome must be
recognised by the id its proteome uses, otherwise it is counted as a non-member hit.

The committed map (`kiss_positive_proteome_map.tsv` in the PR 6 results record) holds one row per Kiss positive with
`status` mapped, not_in_kept_proteomes or unmapped. Nothing here touches the network or a sequence.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

KEPT_BY_TAXON = {
    "330879": "Afum_Af293_UniProt",  # Aspergillus fumigatus Af293
    "237561": "Calb_SC5314",  # Candida albicans SC5314
    "246410": "Cimm_RS",  # Coccidioides immitis RS
    "559292": "Scer_S288C",  # Saccharomyces cerevisiae S288C
}
# JGI project prefixes of the archived Kiss JGI positives; none of these species is a kept proteome.
JGI_PROJECT_SPECIES = {
    "Neucr2": "Neurospora crassa",
    "Ustma2_2": "Mycosarcoma maydis (Ustilago maydis)",
    "Ashgo1_1": "Eremothecium gossypii (Ashbya gossypii)",
}
STATUSES = ("mapped", "not_in_kept_proteomes", "unmapped")
MAP_COLUMNS = [
    "accession",
    "cluster_id",
    "organism_id",
    "organism",
    "status",
    "proteome",
    "proteome_id",
    "uniprot_ordered_locus",
    "basis",
]
_CALB = re.compile(r"^CAALFM_(C[1-7]|CR)(\d{5})([WC])A$")
_SCER = re.compile(r"^Y[A-P][LR]\d{3}[WC](?:-[A-Z])?$")
_XP_VERSION = re.compile(r"\.\d+$")


def calb_a22_id(ordered_locus: str) -> str | None:
    """`CAALFM_C206100WA` -> `C2_06100W_A`; None unless exactly one CAALFM allele-A name is present."""
    names = [t for t in ordered_locus.split() if t.startswith("CAALFM_")]
    if len(names) != 1:
        return None
    m = _CALB.match(names[0])
    return f"{m.group(1)}_{m.group(2)}{m.group(3)}_A" if m else None


def scer_systematic_id(ordered_locus: str) -> str | None:
    """The one SGD systematic name (`YNL238W`) among the ordered-locus names, else None."""
    names = [t for t in ordered_locus.split() if _SCER.match(t)]
    return names[0] if len(names) == 1 else None


def normalize_proteome_id(proteome: str, header_id: str) -> str:
    """The key a proteome hit is matched on: UniProt accession (Afum), ids as written (Calb, Scer), RefSeq without version (Cimm)."""
    if proteome == "Afum_Af293_UniProt":
        parts = header_id.split("|")
        return parts[1] if len(parts) >= 3 else header_id
    if proteome == "Cimm_RS":
        return _XP_VERSION.sub("", header_id)
    return header_id


def map_one(accession: str, cluster_id: str, uniprot: dict | None) -> dict:
    """One map row. `uniprot` is the UniProt TSV row (columns as downloaded) or None when there is none."""
    row = dict.fromkeys(MAP_COLUMNS, "")
    row.update(accession=accession, cluster_id=cluster_id)
    if accession.startswith("jgi:"):
        project = accession.split(":")[1] if accession.count(":") == 2 else ""
        species = JGI_PROJECT_SPECIES.get(project)
        if species:
            row.update(
                organism=species,
                status="not_in_kept_proteomes",
                basis=f"JGI project {project} is a species with no kept proteome",
            )
        else:
            row.update(status="unmapped", basis="JGI id with an unknown project prefix")
        return row
    if not uniprot:
        row.update(status="unmapped", basis="no UniProt entry returned")
        return row
    taxon = uniprot["Organism (ID)"]
    locus = uniprot["Gene Names (ordered locus)"]
    row.update(organism_id=taxon, organism=uniprot["Organism"], uniprot_ordered_locus=locus)
    proteome = KEPT_BY_TAXON.get(taxon)
    if proteome is None:
        row.update(
            status="not_in_kept_proteomes", basis=f"UniProt organism {taxon} has no kept proteome"
        )
    elif proteome == "Afum_Af293_UniProt":
        row.update(
            status="mapped",
            proteome=proteome,
            proteome_id=accession,
            basis="UniProt accession (sp|ACC| or tr|ACC| header)",
        )
    elif proteome == "Calb_SC5314":
        pid = calb_a22_id(locus)
        row.update(
            status="mapped" if pid else "unmapped",
            proteome=proteome,
            proteome_id=pid or "",
            basis="CAALFM ordered-locus name rewritten as CGD Assembly 22 id"
            if pid
            else "no single CAALFM ordered-locus name",
        )
    elif proteome == "Scer_S288C":
        pid = scer_systematic_id(locus)
        row.update(
            status="mapped" if pid else "unmapped",
            proteome=proteome,
            proteome_id=pid or "",
            basis="SGD systematic name from the UniProt ordered-locus field"
            if pid
            else "no single systematic name",
        )
    else:  # Cimm_RS: RefSeq ids are not in the UniProt export, so a positive here could not be mapped
        row.update(
            status="unmapped",
            proteome=proteome,
            basis="Cimm_RS needs a RefSeq id; none is mapped by this method",
        )
    return row


def build_map(members: list[dict], uniprot_rows: list[dict]) -> list[dict]:
    """One row per Kiss positive of `members_split.tsv`, in file order."""
    by_acc = {r["Entry"]: r for r in uniprot_rows}
    return [
        map_one(m["accession"], m["cluster_id"], by_acc.get(m["accession"]))
        for m in members
        if m["class"] == "positive"
    ]


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_map(path: Path, rows: list[dict]) -> None:
    with Path(path).open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MAP_COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in MAP_COLUMNS])


def lookup_from_map(rows: list[dict]) -> dict[tuple[str, str], str]:
    """(proteome, normalised proteome id) -> Kiss accession, for the mapped rows only."""
    out: dict[tuple[str, str], str] = {}
    for r in rows:
        if r["status"] == "mapped":
            key = (r["proteome"], normalize_proteome_id(r["proteome"], r["proteome_id"]))
            if key in out and out[key] != r["accession"]:
                raise ValueError(f"two Kiss positives map to {key}: {out[key]}, {r['accession']}")
            out[key] = r["accession"]
    return out


def classify_hit(
    proteome: str,
    header_id: str,
    own_accessions: set[str],
    lookup: dict[tuple[str, str], str],
) -> str:
    """`own` (an own member of the family), `other_kiss_positive` (a positive of another cluster) or `non_kiss`."""
    acc = lookup.get((proteome, normalize_proteome_id(proteome, header_id)))
    if acc is None:
        return "non_kiss"
    return "own" if acc in own_accessions else "other_kiss_positive"
