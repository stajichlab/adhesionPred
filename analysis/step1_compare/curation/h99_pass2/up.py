import sys
import urllib.parse
import urllib.request


def q(query):
    u = (
        "https://rest.uniprot.org/uniprotkb/search?query="
        + urllib.parse.quote(query)
        + "&fields=accession,gene_names,gene_orf,protein_name,length,reviewed&format=tsv&size=25"
    )
    print("Q:", query)
    print(urllib.request.urlopen(u, timeout=60).read().decode())


for a in sys.argv[1:]:
    q(a)
