# ruff: noqa
"""Show UniProt TSV rows; --exp keeps rows whose function comment has ECO:0000269."""

import csv
import re
import sys

exp = "--exp" in sys.argv
width = 400
for a in sys.argv[1:]:
    if a.startswith("--w="):
        width = int(a[4:])
files = [a for a in sys.argv[1:] if not a.startswith("--")]
ADH = re.compile(
    r"adhe|attach|adher|agglutin|flocc|bind(s|ing)? to (host|epitheli|lamin|fibron|collagen|plastic)",
    re.I,
)
for fn in files:
    for r in csv.DictReader(open(fn), delimiter="\t"):
        f = r["Function [CC]"]
        if exp and "ECO:0000269" not in f:
            continue
        lin = r["Taxonomic lineage"]
        clade = (
            "Basidio"
            if "Basidiomycota" in lin
            else "Onyg"
            if "Onygenales" in lin
            else "Eurot"
            if "Eurotiales" in lin
            else "Pezizo"
            if "Pezizomycotina" in lin
            else "Taphrino"
            if "Taphrinomycotina" in lin
            else "Sacch"
            if "Saccharomycotina" in lin
            else "other"
        )
        go = r.get("Gene Ontology (biological process)", "")
        flag = "ADH!" if ADH.search(f) or ADH.search(go) else ""
        org = " ".join(r["Organism"].split()[:2])
        print(
            "\t".join(
                [
                    r["Entry"],
                    r["Reviewed"][:3],
                    r["Gene Names (primary)"],
                    org,
                    r["Organism (ID)"],
                    clade,
                    r["Length"],
                    "SP" if r["Signal peptide"] else "-",
                    r["Pfam"].rstrip(";"),
                    flag,
                    f[:width],
                ]
            )
        )
