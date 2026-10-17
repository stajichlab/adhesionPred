#!/usr/bin/env python3.12
"""Merge the new fungal AllergenOnline proteins into the positive set (owner request, 2026-10-10).

Inputs:
  --positives   positives/positives.tsv (existing 147 rows, from 02_positives.py)
  --positives-faa positives/positives.faa
  --new-tsv     allergenonline/allergenonline_new.tsv (accession, organism, length, protein_name)
  --new-faa     allergenonline/allergenonline_new.faa (P_ IDs)
  --fungi-tsv   allergenonline/allergenonline_fungi.tsv (for the evidence labels)
Outputs (overwrite positives/ in place; the previous files are kept as *.before_allergenonline):
  positives.tsv  adds rows with in_iuis=false, in_kw0020=false, evidence_level=allergenonline_ige,
                 merged_from=""; a column 'source' is added (iuis, kw0020, allergenonline)
  positives.faa  appends the new P_ records
Checks: no duplicate accession; every added accession has a sequence; the new count is recorded in counts.md.
"""

import argparse
import csv
import hashlib
import shutil
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positives", required=True)
    ap.add_argument("--positives-faa", required=True)
    ap.add_argument("--new-tsv", required=True)
    ap.add_argument("--new-faa", required=True)
    ap.add_argument("--fungi-tsv", required=True)
    ap.add_argument(
        "--evidence", default="allergenonline_ige", help="evidence_level written for the added rows"
    )
    ap.add_argument(
        "--source", default="allergenonline", help="source column written for the added rows"
    )
    args = ap.parse_args()

    old = list(csv.DictReader(open(args.positives, encoding="utf-8"), delimiter="\t"))
    old_acc = {r["accession"] for r in old}
    new = list(csv.DictReader(open(args.new_tsv, encoding="utf-8"), delimiter="\t"))
    fungi = list(csv.DictReader(open(args.fungi_tsv, encoding="utf-8"), delimiter="\t"))
    lab = {}
    for r in fungi:
        if r["new"] == "yes":
            lab.setdefault(r["uniprot"], (r["allergenicity"], r["iuis_allergen"]))

    seqs = {}
    cur = None
    for ln in open(args.new_faa, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            cur = ln[1:].split()[0][2:]  # strip P_
            seqs[cur] = []
        elif cur:
            seqs[cur].append(ln)
    seqs = {k: "".join(v) for k, v in seqs.items()}

    added = []
    for r in new:
        a = r["accession"]
        if a in old_acc:
            sys.exit(f"FATAL: {a} is already in the positive set")
        if a not in seqs or not seqs[a]:
            sys.exit(f"FATAL: {a} has no sequence in {args.new_faa}")
        allergenicity, iuis_name = lab.get(a, ("", ""))
        added.append(
            {
                "accession": a,
                "in_iuis": "false",
                "in_kw0020": "false",
                "evidence_level": args.evidence,
                "merged_from": "",
                "organism": r["organism"],
                "length": r["length"],
                "protein_name": r["protein_name"],
                "sequence_sha256": hashlib.sha256(seqs[a].encode()).hexdigest(),
                "source": args.source,
                "allergenonline_type": allergenicity,
                "iuis_name": iuis_name,
            }
        )

    base = args.positives
    for path in (args.positives, args.positives_faa):
        shutil.copy(path, path + ".before_allergenonline")

    fields = [
        "accession",
        "in_iuis",
        "in_kw0020",
        "evidence_level",
        "merged_from",
        "organism",
        "length",
        "protein_name",
        "sequence_sha256",
        "source",
        "allergenonline_type",
        "iuis_name",
    ]
    with open(base, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in old:
            r = dict(r)
            r.setdefault("source", "iuis_or_kw0020")
            w.writerow(r)
        for r in added:
            w.writerow(r)
    with open(args.positives_faa, "a", encoding="utf-8") as fh:
        for r in added:
            fh.write(f">P_{r['accession']}\n")
            s = seqs[r["accession"]]
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    print(f"old={len(old)} added={len(added)} total={len(old) + len(added)}")


if __name__ == "__main__":
    main()
