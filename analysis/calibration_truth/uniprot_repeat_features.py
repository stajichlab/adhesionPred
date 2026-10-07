"""Fetch UniProt repeat annotations for the accessions of the repeat-mechanism curation table.

For each accession it records how many `Repeat` features the UniProt entry has, the evidence codes of
those features, and how many Region, Compositional bias or Domain features mention "repeat".
This is an evidence tier of its own ("database annotation"). It is independent of the agents that
filled the table and of our repeat detectors. It is not an experimental result: most repeat features
carry no evidence code or a rule-based or curator-inferred code. It under-counts: well-known repeat
proteins can lack a repeat feature.

Usage: python3.12 uniprot_repeat_features.py --table curation_table.tsv --out uniprot_repeat_features.tsv
UniProt is queried live, so a re-run can change the numbers. The fetch date is written to the table.
"""

import argparse
import csv
import datetime
import json
import re
import subprocess
import time
import urllib.parse
from collections import Counter

COLUMNS = [
    "accession",
    "reviewed",
    "length",
    "n_repeat_features",
    "repeat_evidence_codes",
    "n_region_mentions_repeat",
    "region_descriptions",
    "fetched",
]


def fetch_entries(accessions, batch=40):
    for i in range(0, len(accessions), batch):
        query = " OR ".join(f"accession:{a}" for a in accessions[i : i + batch])
        url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=" + urllib.parse.quote(
            query
        )
        text = subprocess.run(
            ["curl", "-s", "-m", "90", url], capture_output=True, text=True, check=True
        ).stdout
        yield from json.loads(text)["results"]
        time.sleep(0.3)


def summarise(entry, fetched):
    feats = entry.get("features", [])
    repeats = [f for f in feats if f["type"] == "Repeat"]
    codes = Counter()
    for f in repeats:
        evs = f.get("evidences") or [{}]
        for ev in evs:
            codes[ev.get("evidenceCode", "none")] += 1
    regions = [
        f.get("description", "")
        for f in feats
        if f["type"] in ("Region", "Compositional bias", "Domain")
        and re.search("repeat", f.get("description", ""), re.IGNORECASE)
    ]
    return {
        "accession": entry["primaryAccession"],
        "reviewed": "yes" if entry["entryType"].startswith("UniProtKB reviewed") else "no",
        "length": entry["sequence"]["length"],
        "n_repeat_features": len(repeats),
        "repeat_evidence_codes": ";".join(f"{k}:{v}" for k, v in sorted(codes.items())),
        "n_region_mentions_repeat": len(regions),
        "region_descriptions": " | ".join(regions[:3]),
        "fetched": fetched,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--table", required=True, help="a curation table with an `accession` column")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    with open(args.table, encoding="utf-8-sig", newline="") as fh:
        accessions = sorted(
            {r["accession"].strip() for r in csv.DictReader(fh, delimiter="\t")} - {""}
        )
    today = datetime.date.today().isoformat()
    found = {}
    for entry in fetch_entries(accessions):
        found[entry["primaryAccession"]] = summarise(entry, today)
    missing = [a for a in accessions if a not in found]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for a in accessions:
            if a in found:
                w.writerow(found[a])
    print(f"{len(found)} of {len(accessions)} accessions fetched; missing: {missing[:10]}")


if __name__ == "__main__":
    main()
