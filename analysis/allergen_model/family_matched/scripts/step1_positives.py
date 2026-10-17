"""Step 1: Pfam domains for each positive. Writes positive_pfam.tsv and skip list."""

import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
from uni import get_all, pfams

ROOT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model"
OUT = ROOT + "/family_matched"
pos = list(csv.DictReader(open(ROOT + "/positives/positives.tsv"), delimiter="\t"))
accs = [r["accession"] for r in pos]
# UniProt accession pattern (6 or 10 char). NCBI protein IDs contain a dot.
UNI = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
uni_accs = [a for a in accs if UNI.match(a)]
ncbi_accs = [a for a in accs if a not in set(uni_accs)]
print(
    "positives",
    len(accs),
    "uniprot-pattern",
    len(uni_accs),
    "non-uniprot",
    len(ncbi_accs),
    file=sys.stderr,
)

found = {}
release = None
for i in range(0, len(uni_accs), 40):
    chunk = uni_accs[i : i + 40]
    q = " OR ".join(f"accession:{a}" for a in chunk)
    rows, release = get_all(
        f"({q})", "accession,xref_pfam,reviewed,organism_name,length,protein_name"
    )
    for r in rows:
        found[r["Entry"]] = r
    print("chunk", i, len(rows), file=sys.stderr)

with open(OUT + "/positive_pfam.tsv", "w") as f:
    f.write("accession\tin_uniprot\treviewed\tpfam\tlength\n")
    for a in uni_accs:
        r = found.get(a)
        if r is None:
            f.write(f"{a}\tfalse\t\t\t\n")
        else:
            f.write(f"{a}\ttrue\t{r['Reviewed']}\t{','.join(pfams(r['Pfam']))}\t{r['Length']}\n")
with open(OUT + "/positive_ncbi.txt", "w") as f:
    f.write("\n".join(ncbi_accs) + "\n")
with open(OUT + "/release.txt", "w") as f:
    f.write(str(release) + "\n")
print("release", release)
