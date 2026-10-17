#!/usr/bin/env python3.12
"""Fetch UniProt JSON records for a list of accessions.

Usage: fetch_uniprot.py ACC_LIST OUT.json.gz
ACC_LIST: one accession per line (# comments allowed). Output: a JSON list of entries.
The release is printed (X-UniProt-Release headers). UniProt data is CC BY 4.0.
"""

import gzip
import json
import sys
import time
import urllib.parse
import urllib.request

accs = []
for line in open(sys.argv[1]):
    line = line.split("#")[0].strip()
    if line and line not in accs:
        accs.append(line)
out = []
for i in range(0, len(accs), 40):
    chunk = accs[i : i + 40]
    q = "accession:(" + " OR ".join(chunk) + ")"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "size": 100}
    )
    for _attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                rel = r.headers.get("X-UniProt-Release"), r.headers.get("X-UniProt-Release-Date")
                out.extend(json.load(r)["results"])
            break
        except Exception as exc:
            print("retry", exc, file=sys.stderr)
            time.sleep(3)
    else:
        raise SystemExit("failed")
got = {e["primaryAccession"] for e in out}
secondary = {}
for e in out:
    for s in e.get("secondaryAccessions", []):
        secondary[s] = e["primaryAccession"]
missing = [a for a in accs if a not in got]
print("requested", len(accs), "got", len(out), "release", rel)
print("missing (not primary):", missing)
with gzip.open(sys.argv[2], "wt") as fh:
    json.dump(out, fh)
