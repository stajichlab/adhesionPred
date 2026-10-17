#!/usr/bin/env python3.12
"""Fetch UniLectin3D fungal structures (documents) and RCSB ligand data.

Step 1: query the UniLectin search API (POST, module unilectin3d, rows 100) for the phrase
"Fungal and yeast lectins" and for every class, family and species name that comes back
(the API returns at most 100 documents per call and ignores `start`). Keep documents with
origin "Fungal and yeast lectins".
Step 2: for every PDB id, query the RCSB GraphQL API for the non-polymer and branched
(carbohydrate) entities and the primary citation.

Output (--out-dir): unilectin_docs.json.gz, rcsb_ligands.json.gz
Usage: fetch_unilectin.py --out-dir DIR
"""

import argparse
import gzip
import json
import sys
import time
import urllib.request

UL = "https://unilectin.unige.ch/api/search-module"
RCSB = "https://data.rcsb.org/graphql"

GQL = """query($id:String!){ entry(entry_id:$id){
  rcsb_primary_citation{pdbx_database_id_PubMed title}
  nonpolymer_entities{pdbx_entity_nonpoly{comp_id name}}
  branched_entities{pdbx_entity_branch{type} rcsb_branched_entity_name_com{name}
    pdbx_entity_branch_descriptor{type descriptor}}
}}"""


def post(url, payload, tries=3):
    data = json.dumps(payload).encode()
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, data=data, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as exc:  # network errors: retry
            print("retry", i, exc, file=sys.stderr)
            time.sleep(3)
    raise RuntimeError("failed: " + str(payload)[:100])


def ul_query(q):
    return post(UL, {"module": "unilectin3d", "rows": 100, "query": q})["response"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    docs = {}
    seen = set()
    queue = ["Fungal and yeast lectins"]
    while queue:
        q = queue.pop(0)
        if q in seen:
            continue
        seen.add(q)
        resp = ul_query(q)
        for d in resp["docs"]:
            if d.get("origin") != "Fungal and yeast lectins":
                continue
            docs[d["pdb"]] = d
            for k in ("class", "family", "species"):
                v = d.get(k)
                if v and v not in seen:
                    queue.append(v)
    print("structures", len(docs), "queries", len(seen))
    with gzip.open(a.out_dir + "/unilectin_docs.json.gz", "wt") as fh:
        json.dump(list(docs.values()), fh)
    out = {}
    for pdb in sorted(docs):
        r = post(RCSB, {"query": GQL, "variables": {"id": pdb}})
        out[pdb] = (r.get("data") or {}).get("entry")
    with gzip.open(a.out_dir + "/rcsb_ligands.json.gz", "wt") as fh:
        json.dump(out, fh)
    print("rcsb entries", sum(1 for v in out.values() if v))


main()
