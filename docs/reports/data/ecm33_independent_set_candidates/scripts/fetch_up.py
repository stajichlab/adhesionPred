# ruff: noqa
# Scoping script kept exactly as run on 2026-10-23 (stdlib only).
import urllib.request, urllib.parse, re, sys

q = sys.argv[1]
out = sys.argv[2]
fields = sys.argv[3]
u = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
    {"query": q, "fields": fields, "format": "tsv", "size": 500}
)
first = True
n = 0
with open(out, "w") as f:
    while u:
        r = urllib.request.urlopen(
            urllib.request.Request(u, headers={"User-Agent": "adhesionPred-scoping/1.0"}),
            timeout=100,
        )
        t = r.read().decode().split("\n")
        if first:
            f.write(t[0] + "\n")
            first = False
        body = [x for x in t[1:] if x]
        n += len(body)
        f.write("\n".join(body) + "\n")
        l = r.headers.get("Link")
        m = re.search(r'<([^>]+)>; rel="next"', l or "")
        u = m.group(1) if m else None
print(n, r.headers.get("X-Total-Results"))
