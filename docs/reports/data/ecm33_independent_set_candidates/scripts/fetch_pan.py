# ruff: noqa
# Scoping script kept exactly as run on 2026-10-23 (stdlib only).
import json, urllib.request, csv, time

u = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/panther/PTHR31018/taxonomy/uniprot/4751?page_size=100"
rows = []


def get(u):
    for k in range(8):
        try:
            return json.load(
                urllib.request.urlopen(
                    urllib.request.Request(u, headers={"User-Agent": "adhesionPred-scoping/1.0"}),
                    timeout=40,
                )
            )
        except Exception as e:
            print("retry", k, e, flush=True)
            time.sleep(3)
    raise SystemExit("fail " + u)


while u:
    d = get(u)
    for r in d["results"]:
        m = r["metadata"]
        e = r["entries"][0]["entry_protein_locations"][0]
        rows.append(
            [
                m["accession"],
                m["source_database"],
                m["source_organism"]["taxId"],
                m["source_organism"]["scientificName"],
                m.get("gene") or "",
                m["name"],
                m["length"],
                e.get("model", ""),
                e.get("score", ""),
                e["fragments"][0]["start"],
                e["fragments"][-1]["end"],
            ]
        )
    u = d.get("next")
    print(len(rows), flush=True)
w = csv.writer(open("pan.tsv", "w"), delimiter="\t")
w.writerow("accession reviewed taxid organism gene name length model score start end".split())
w.writerows(rows)
print("done", len(rows))
