import gzip
import sys
import urllib.parse
import urllib.request

D = "/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare/downloads/UP000000561.fasta.gz"
hdr = {}
for ln in gzip.open(D, "rt"):
    if ln[0] == ">":
        acc = ln.split("|")[1]
        hdr[acc] = ln.strip()


def lookup(q):
    u = (
        "https://rest.uniprot.org/uniprotkb/search?query="
        + urllib.parse.quote(f"(organism_id:5270 OR organism_id:237631) AND ({q})")
        + "&fields=accession,gene_orf,gene_names,protein_name&format=tsv&size=10"
    )
    r = (
        urllib.request.urlopen(
            urllib.request.Request(u, headers={"User-Agent": "curation/1.0"}), timeout=60
        )
        .read()
        .decode()
    )
    return r.strip().split("\n")[1:]


for q in sys.argv[1:]:
    for r in lookup(q):
        a = r.split("\t")[0]
        print(q, "|", r, "| IN_PROTEOME" if a in hdr else "| NOT_IN_PROTEOME")
