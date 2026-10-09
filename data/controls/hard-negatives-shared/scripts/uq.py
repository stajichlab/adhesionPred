# ruff: noqa
"""Query UniProt REST, save TSV. Usage: uq.py NAME 'QUERY' OUTDIR"""

import subprocess
import sys
import urllib.parse

FIELDS = (
    "accession,reviewed,gene_primary,protein_name,organism_name,organism_id,lineage,length,"
    "cc_function,go_p,go_c,xref_pfam,ft_signal,keyword,lit_pubmed_id"
)
name, q, out = sys.argv[1], sys.argv[2], sys.argv[3]
url = (
    "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields="
    + FIELDS
    + "&query="
    + urllib.parse.quote(q)
)
txt = subprocess.run(
    ["curl", "-s", "-m", "120", url], capture_output=True, text=True, check=True
).stdout
open(f"{out}/{name}.tsv", "w").write(txt)
print(name, max(0, len(txt.splitlines()) - 1))
