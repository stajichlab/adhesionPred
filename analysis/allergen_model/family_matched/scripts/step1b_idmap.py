"""Step 1b: map NCBI protein IDs to UniProtKB via the UniProt ID mapping API (stdlib)."""

import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(__file__))
OUT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/family_matched"
ids = [x.strip() for x in open(OUT + "/positive_ncbi.txt") if x.strip()]


def post(url, data):
    body = urllib.parse.urlencode(data).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=body), timeout=120) as r:
        return json.loads(r.read())


def get_json(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.loads(r.read())


result = {}
for db in ["EMBL-GenBank-DDBJ_CDS", "RefSeq_Protein"]:
    todo = [i for i in ids if i not in result]
    if not todo:
        break
    job = post(
        "https://rest.uniprot.org/idmapping/run",
        {"from": db, "to": "UniProtKB", "ids": ",".join(todo)},
    )
    jid = job["jobId"]
    for _ in range(120):
        st = get_json(f"https://rest.uniprot.org/idmapping/status/{jid}")
        if "results" in st or "failedIds" in st and "jobStatus" not in st:
            break
        if st.get("jobStatus") == "FINISHED":
            break
        time.sleep(2)
    res = get_json(f"https://rest.uniprot.org/idmapping/results/{jid}?size=500")
    for m in res.get("results", []):
        frm, to = m["from"], m["to"]
        acc = to["primaryAccession"] if isinstance(to, dict) else to
        result.setdefault(frm, (acc, db))
    print(db, "mapped so far", len(result), file=sys.stderr)
with open(OUT + "/positive_ncbi_map.tsv", "w") as f:
    f.write("ncbi_id\tuniprot\tmap_db\n")
    for i in ids:
        if i in result:
            f.write(f"{i}\t{result[i][0]}\t{result[i][1]}\n")
        else:
            f.write(f"{i}\t\t\n")
print("mapped", len(result), "of", len(ids))
