# ruff: noqa
# Scoping script kept exactly as run on 2026-10-23 (stdlib only).
import csv, json, sys, urllib.request, urllib.parse, time

ids = sorted({l.strip() for l in open(sys.argv[1]) if l.strip()}, key=int)
out = {}


def get(u):
    for k in range(6):
        try:
            return json.load(
                urllib.request.urlopen(
                    urllib.request.Request(u, headers={"User-Agent": "adhesionPred-scoping/1.0"}),
                    timeout=40,
                )
            )
        except Exception as e:
            print("retry", e, file=sys.stderr)
            time.sleep(2)
    raise SystemExit("fail")


for i in range(0, len(ids), 40):
    ch = ids[i : i + 40]
    q = " OR ".join("tax_id:" + x for x in ch)
    d = get(
        "https://rest.uniprot.org/taxonomy/search?"
        + urllib.parse.urlencode(
            {"query": q, "fields": "id,rank,scientific_name,lineage", "format": "json", "size": 100}
        )
    )
    for r in d["results"]:
        lin = {}
        for l in r.get("lineage", []):
            lin[l["rank"]] = l["scientificName"]
        if r.get("rank"):
            lin[r["rank"]] = r["scientificName"]
        phylum = lin.get("phylum", "")
        sub = lin.get("subphylum", "")
        cls = lin.get("class", "")
        if phylum == "Ascomycota":
            fine = cls if sub == "Pezizomycotina" else sub
        elif phylum == "Basidiomycota":
            fine = cls
        else:
            fine = phylum
        out[str(r["taxonId"])] = (fine, phylum)
json.dump(out, open(sys.argv[2], "w"))
print(len(ids), len(out))
