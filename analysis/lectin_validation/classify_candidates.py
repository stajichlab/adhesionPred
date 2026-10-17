#!/usr/bin/env python3.12
"""Screen the fetched UniProt entries for experimental sugar-binding annotation.

Rule (written before the screen was read; spec 4.1 tier T2 for Swiss-Prot):
  reviewed entry, length >= 50, a FUNCTION text whose evidence includes ECO:0000269 and that
  matches both a binding word and a sugar word. The screen only proposes. A person reads each
  proposal. Nothing here uses a Pfam hit of the models under test: UniProt's own Pfam
  cross-references are printed for information.
Usage: classify_candidates.py UNIPROT_JSON_GZ UNILECTIN_TSV OUT_TSV
"""

import csv
import gzip
import json
import re
import sys

BIND = re.compile(
    r"\b(binds?|binding|specific(ally)? (for|to)|recogni[sz]es?|agglutinat\w+|hemagglutin\w+|sequester\w*)\b",
    re.I,
)
SUGAR = re.compile(
    r"(chitin|fucos|lactos|galactos|mannos|glcnac|acetylglucosamine|oligosaccharide|glycan|sialic|"
    r"carbohydrate|sugar|laminin|lactosamine|galnac|arabinos|inulin|xylos|glucos)",
    re.I,
)

ul = {r["uniprot"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
rows = []
for e in json.load(gzip.open(sys.argv[1], "rt")):
    acc = e["primaryAccession"]
    reviewed = e["entryType"].startswith("UniProtKB reviewed")
    length = e["sequence"]["length"]
    kws = [k["name"] for k in e.get("keywords", [])]
    pfam = sorted(
        {x["id"] for x in e.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"}
    )
    best = ""
    eco = False
    pmids = set()
    for c in e.get("comments", []):
        if c["commentType"] != "FUNCTION":
            continue
        for t in c.get("texts", []):
            ev = {x["evidenceCode"] for x in t.get("evidences", [])}
            if BIND.search(t["value"]) and SUGAR.search(t["value"]) and "ECO:0000269" in ev:
                eco = True
                best = t["value"]
                pmids |= {x["id"] for x in t["evidences"] if x["evidenceCode"] == "ECO:0000269"}
    rows.append(
        {
            "accession": acc,
            "reviewed": "yes" if reviewed else "no",
            "length": length,
            "organism": e["organism"]["scientificName"],
            "taxon_id": e["organism"]["taxonId"],
            "in_unilectin3d": "yes" if acc in ul else "no",
            "keyword_lectin": "yes" if "Lectin" in kws else "no",
            "uniprot_pfam": ",".join(pfam),
            "eco269_binding_text": "yes" if eco else "no",
            "pmids": ",".join(sorted(pmids)),
            "text": best[:300],
        }
    )
with open(sys.argv[3], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(
    len(rows), "entries;", sum(r["eco269_binding_text"] == "yes" for r in rows), "pass the screen"
)
