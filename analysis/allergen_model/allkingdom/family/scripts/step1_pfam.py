"""Step 1: Pfam accessions for each positive (UniProt xref_pfam, batches of 40, OR query)."""

import sys
import time
import urllib.parse

sys.path.insert(
    0,
    "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/allkingdom/family/scripts",
)
from common import FAM, http_get, positives

pos = positives()
accs = sorted({p["accession"] for p in pos})
pfam = {}
found = set()
for i in range(0, len(accs), 40):
    batch = accs[i : i + 40]
    q = " OR ".join(f"accession:{a}" for a in batch)
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": q, "fields": "accession,xref_pfam,reviewed", "format": "tsv", "size": "500"}
    )
    body, _ = http_get(url)
    lines = body.splitlines()[1:]
    for line in lines:
        if not line:
            continue
        entry, pf, rev = (line.split("\t") + ["", ""])[:3]
        found.add(entry)
        pfam[entry] = [x for x in pf.split(";") if x]
    time.sleep(0.3)

with open(FAM + "/raw/positive_pfam.tsv", "w") as fh:
    fh.write("accession\tgroup\tpfam\n")
    group = {p["accession"]: p["group"] for p in pos}
    for a in accs:
        fh.write(f"{a}\t{group[a]}\t{','.join(pfam.get(a, []))}\n")

n_nopfam = sum(1 for a in accs if not pfam.get(a))
n_notfound = len(set(accs) - found)
print("positives", len(pos), "unique accessions", len(accs))
print("returned by UniProt", len(found), "not returned", n_notfound)
print("with Pfam", len(accs) - n_nopfam, "without Pfam", n_nopfam)
