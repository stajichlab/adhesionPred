import json
import sys
import urllib.parse

from fetch import get

q = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 25
u = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
    + urllib.parse.quote(q)
    + f"&format=json&pageSize={n}&resultType=lite"
)
j = json.loads(get(u))
print("HITS", j["hitCount"])
for r in j["resultList"]["result"]:
    print(
        r.get("pmid"),
        r.get("pmcid"),
        r.get("pubYear"),
        r.get("title")[:95],
        "|",
        r.get("journalTitle", "")[:20],
    )
