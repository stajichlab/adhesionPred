import json
import sys
import urllib.parse
import urllib.request


def s(q, n=25):
    u = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
        + urllib.parse.quote(q)
        + f"&format=json&resultType=lite&pageSize={n}"
    )
    j = json.loads(urllib.request.urlopen(u, timeout=60).read())
    print("Q:", q, "| hitCount", j["hitCount"])
    for r in j["resultList"]["result"]:
        print(" ", r.get("pmid"), r.get("pmcid"), r.get("pubYear"), r.get("title", "")[:110])


for q in sys.argv[1:]:
    s(q)
