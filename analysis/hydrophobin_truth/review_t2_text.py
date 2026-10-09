#!/usr/bin/python3.12
"""Task L1b: for every T2 entry, take the experimentally supported FUNCTION, SUBCELLULAR LOCATION and SUBUNIT text and
classify it as hydrophobin_related or other. Nothing is relabelled. Usage: review_t2_text.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

RELATED = re.compile(
    r"hydrophobin|rodlet|hydrophobic|surface tension|spore wall|aerial|amphipath|self-?assembl|wettab|surface-active|conidia|hyphae|fruiting",
    re.I,
)
TYPES = {"FUNCTION", "SUBCELLULAR LOCATION", "SUBUNIT"}


def classify(text):
    return "hydrophobin_related" if RELATED.search(text or "") else "other"


def experimental_text(entry):
    out = []
    for c in entry.get("comments", []):
        if c.get("commentType") not in TYPES:
            continue
        for t in c.get("texts", []):
            if any(
                e.get("evidenceCode") == "ECO:0000269" and e.get("source") == "PubMed"
                for e in t.get("evidences", [])
            ):
                out.append(t["value"])
    return " ".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    entries = {
        e["primaryAccession"]: e for e in json.load(gzip.open(d / "uniprot_entries.json.gz", "rt"))
    }
    tiers = {r["accession"]: r for r in csv.DictReader(open(d / "truth_all.tsv"), delimiter="\t")}
    n_other = 0
    with open(d / "t2_text_review.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "entry_name", "class", "experimental_text"])
        for acc, r in tiers.items():
            if r["tier"] != "T2":
                continue
            text = experimental_text(entries[acc])
            cl = classify(text)
            n_other += cl == "other"
            w.writerow([acc, r["entry_name"], cl, text[:300]])
            if cl == "other":
                print(acc, r["entry_name"], "|", text[:160])
    print("T2 entries with a text that is not hydrophobin_related:", n_other)
    return 0


if __name__ == "__main__":
    sys.exit(main())
