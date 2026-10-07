"""Audit: are there repeat-bearing adhesion or cell-wall proteins outside Saccharomycotina in UniProt?

Agent-C searched PubMed abstracts for protein-level adhesion evidence in non-Saccharomycotina fungi
(task 08) and found no new class-2a repeat adhesin. This script checks that result from another side.
It asks UniProt for fungal proteins outside Saccharomycotina (taxon 147537) that have a `Repeat`
feature and an adhesion or cell-wall annotation, and it records how strong each annotation is.

These are CANDIDATES, not controls. A UniProt record gives no proof that the protein is an adhesin
or that the repeats mediate adhesion. Use the output to see where a literature check is worth doing.

Output: one TSV row per UniProt entry. Columns say whether the entry is reviewed, how many repeat
features it has, whether it has a signal peptide and a GPI-anchor feature, the adhesion GO terms
with their evidence codes (experimental or not), and whether the accession is already in the
curation tables or in agent-C's files.

Usage: python3.12 uniprot_nonyeast_repeat_candidates.py --out FILE.tsv [--known FILE ...]
UniProt is queried live; a re-run can change the numbers. The fetch date is written to each row.
"""

import argparse
import csv
import datetime
import json
import subprocess
import urllib.parse
from collections import Counter

FUNGI = 4751
SACCHAROMYCOTINA = 147537
# KW-0130 Cell adhesion, KW-0134 Cell wall, GO:0007155 cell adhesion (UniProt expands to children).
QUERY = (
    f"(taxonomy_id:{FUNGI} AND NOT taxonomy_id:{SACCHAROMYCOTINA}) AND ft_repeat:* AND "
    "(keyword:KW-0130 OR keyword:KW-0134 OR go:0007155 OR protein_name:adhesin OR "
    "protein_name:agglutinin OR protein_name:flocculin OR protein_name:hydrophobin)"
)
EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}
ADHESION_GO = {"GO:0007155", "GO:0098609", "GO:0098610", "GO:0044406", "GO:0000128", "GO:0000501"}
COLUMNS = [
    "accession",
    "reviewed",
    "protein_name",
    "gene",
    "organism",
    "taxon_id",
    "phylum_or_class",
    "length",
    "n_repeat_features",
    "repeat_evidence_codes",
    "signal_peptide",
    "gpi_feature",
    "kw_cell_adhesion",
    "kw_cell_wall",
    "adhesion_go",
    "adhesion_go_experimental",
    "name_flags",
    "known_in_curation_tables",
    "fetched",
]


def fetch():
    url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=" + urllib.parse.quote(QUERY)
    text = subprocess.run(
        ["curl", "-s", "-m", "300", url], capture_output=True, text=True, check=True
    ).stdout
    return json.loads(text)["results"]


def lineage_group(entry):
    names = list(entry.get("lineage", []))
    for key in (
        "Basidiomycota",
        "Mucoromycota",
        "Chytridiomycota",
        "Microsporidia",
        "Taphrinomycotina",
        "Pezizomycotina",
        "Blastocladiomycota",
    ):
        if key in names:
            for sub in (
                "Eurotiomycetes",
                "Sordariomycetes",
                "Dothideomycetes",
                "Leotiomycetes",
                "Agaricomycetes",
                "Ustilaginomycotina",
                "Tremellomycetes",
                "Pucciniomycotina",
                "Pezizomycetes",
            ):
                if sub in names:
                    return f"{key}/{sub}"
            return key
    return "other"


def row_of(entry, known, today):
    feats = entry.get("features", [])
    repeats = [f for f in feats if f["type"] == "Repeat"]
    codes = Counter()
    for f in repeats:
        for ev in f.get("evidences") or [{}]:
            codes[ev.get("evidenceCode", "none")] += 1
    kws = {k["id"] for k in entry.get("keywords", [])}
    go_terms, go_exp = [], False
    for x in entry.get("uniProtKBCrossReferences", []):
        if x.get("database") != "GO":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties", [])}
        evidence = props.get("GoEvidenceType", "").split(":")[0]
        if x["id"] in ADHESION_GO:
            go_terms.append(f"{x['id']}:{evidence}")
            go_exp = go_exp or evidence in EXPERIMENTAL
    name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get(
        "value"
    ) or (entry.get("proteinDescription", {}).get("submissionNames") or [{}])[0].get(
        "fullName", {}
    ).get("value", "")
    genes = entry.get("genes") or [{}]
    gene = genes[0].get("geneName", {}).get("value", "") or genes[0].get("orderedLocusNames", [{}])[
        0
    ].get("value", "")
    low = name.lower()
    flags = [w for w in ("adhesin", "agglutinin", "flocculin", "hydrophobin") if w in low]
    acc = entry["primaryAccession"]
    return {
        "accession": acc,
        "reviewed": "yes" if entry["entryType"].startswith("UniProtKB reviewed") else "no",
        "protein_name": name,
        "gene": gene,
        "organism": entry["organism"]["scientificName"],
        "taxon_id": entry["organism"]["taxonId"],
        "phylum_or_class": lineage_group(entry["organism"]),
        "length": entry["sequence"]["length"],
        "n_repeat_features": len(repeats),
        "repeat_evidence_codes": ";".join(f"{k}:{v}" for k, v in sorted(codes.items())),
        "signal_peptide": "yes" if any(f["type"] == "Signal" for f in feats) else "no",
        "gpi_feature": "yes"
        if any(f["type"] == "Lipidation" and "GPI" in f.get("description", "") for f in feats)
        else "no",
        "kw_cell_adhesion": "yes" if "KW-0130" in kws else "no",
        "kw_cell_wall": "yes" if "KW-0134" in kws else "no",
        "adhesion_go": ";".join(go_terms),
        "adhesion_go_experimental": "yes" if go_exp else "no",
        "name_flags": ",".join(flags),
        "known_in_curation_tables": "yes" if acc in known else "no",
        "fetched": today,
    }


def known_accessions(paths):
    known = set()
    for p in paths:
        with open(p, encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("accession"):
                    known.add(r["accession"].strip())
    return known


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--known", nargs="*", default=[], help="tables with an `accession` column")
    args = ap.parse_args()
    known = known_accessions(args.known)
    today = datetime.date.today().isoformat()
    entries = fetch()
    rows = [row_of(e, known, today) for e in entries]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(
            sorted(rows, key=lambda r: (r["phylum_or_class"], r["organism"], r["accession"]))
        )
    print(f"{len(rows)} entries written to {args.out}; query: {QUERY}")


if __name__ == "__main__":
    main()
