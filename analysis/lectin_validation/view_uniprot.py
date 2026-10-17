#!/usr/bin/env python3.12
"""Print a compact view of UniProt JSON entries. Usage: view_uniprot.py JSON_GZ [ACC ...]"""

import gzip
import json
import sys

d = json.load(gzip.open(sys.argv[1], "rt"))
want = set(sys.argv[2:])
for e in d:
    a = e["primaryAccession"]
    if want and a not in want:
        continue
    fn = [c for c in e.get("comments", []) if c["commentType"] == "FUNCTION"]
    txt = " ".join(t["value"] for c in fn for t in c.get("texts", []))
    codes = sorted(
        {
            ev["evidenceCode"] + ":" + ev.get("id", "")
            for c in fn
            for t in c.get("texts", [])
            for ev in t.get("evidences", [])
        }
    )
    print(
        a,
        e["entryType"][:9],
        e["organism"]["taxonId"],
        e["organism"]["scientificName"][:35],
        e["sequence"]["length"],
    )
    print(
        "  ",
        e.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")[:60],
    )
    print("   FN:", txt[:500])
    print("   EV:", ",".join(codes))
