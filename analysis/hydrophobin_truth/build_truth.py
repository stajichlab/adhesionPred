#!/usr/bin/python3.12
"""Build the hydrophobin truth table (task H2 of the hydrophobin validation plan).

Steps (subcommands):
  fetch    fetch UniProt JSON (with evidence codes) for the accessions in sp_hydrophobin_query.tsv
  tiers    assign an evidence tier to each entry and write truth_all.tsv
  split    cluster-wise development/test split from clusters.tsv, written once

Tiers (spec section 4.1):
  T2  reviewed entry, name has hydrophobin or rodlet, and a FUNCTION, SUBCELLULAR LOCATION or SUBUNIT
      comment has experimental evidence (ECO:0000269) from PubMed
  T3  name matches, no such evidence (rule, similarity, or other comment types only)
  none  name does not match: never a positive

T1 (a paper's protein-level evidence added by hand) comes from manual_entries.tsv, with a PMID per row.
Pfam is not used for any tier.

Usage: build_truth.py fetch|tiers|split [--dir DIR]
Reads plain or .gz input.
"""

import argparse
import csv
import gzip
import json
import random
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

NAME_RE = re.compile(r"hydrophobin|rodlet", re.I)
EXPERIMENTAL = "ECO:0000269"
COMMENT_TYPES = {"FUNCTION", "SUBCELLULAR LOCATION", "SUBUNIT"}


def _open(path, mode="rt"):
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def protein_names(entry):
    d = entry.get("proteinDescription", {})
    names = []
    for key in ("recommendedName", "submissionNames", "alternativeNames"):
        v = d.get(key)
        if not v:
            continue
        for n in v if isinstance(v, list) else [v]:
            names.append(n.get("fullName", {}).get("value", ""))
            for s in n.get("shortNames", []):
                names.append(s.get("value", ""))
    return names


def is_named(entry):
    return any(NAME_RE.search(n) for n in protein_names(entry))


def assign_tier(entry):
    """Return (tier, pmids). tier is None when the name does not match."""
    if not is_named(entry):
        return None, []
    reviewed = entry.get("entryType", "").startswith("UniProtKB reviewed")
    pmids = []
    for c in entry.get("comments", []):
        if c.get("commentType") not in COMMENT_TYPES:
            continue
        for t in c.get("texts", []):
            for e in t.get("evidences", []):
                if e.get("evidenceCode") == EXPERIMENTAL and e.get("source") == "PubMed":
                    pmids.append(e["id"])
    if reviewed and pmids:
        return "T2", sorted(set(pmids))
    return "T3", []


def split_clusters(clusters, seed, dev_fraction=0.3):
    """clusters: {protein: cluster}. Whole clusters go to dev or test. Same seed, same split."""
    ids = sorted(set(clusters.values()))
    rng = random.Random(seed)
    rng.shuffle(ids)
    n_dev = max(1, round(len(ids) * dev_fraction))
    dev = set(ids[:n_dev])
    return {p: ("dev" if c in dev else "test") for p, c in clusters.items()}


def fetch(d):
    rows = list(csv.DictReader(open(d / "sp_hydrophobin_query.tsv"), delimiter="\t"))
    accs = [r["Entry"] for r in rows]
    out = []
    for i in range(0, len(accs), 40):
        q = "(" + " OR ".join(f"accession:{a}" for a in accs[i : i + 40]) + ")"
        url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
            {"query": q, "format": "json"}
        )
        with urllib.request.urlopen(url, timeout=120) as r:
            out.extend(json.load(r)["results"])
    got = {e["primaryAccession"] for e in out}
    missing = sorted(set(accs) - got)
    with gzip.open(d / "uniprot_entries.json.gz", "wt") as fh:
        json.dump(out, fh)
    print(f"fetched {len(out)} of {len(accs)} entries; missing {missing}")


def tiers(d):
    entries = json.load(_open(d / "uniprot_entries.json.gz"))
    rows = []
    for e in entries:
        tier, pmids = assign_tier(e)
        seq = e["sequence"]["value"]
        rows.append(
            {
                "accession": e["primaryAccession"],
                "entry_name": e["uniProtkbId"],
                "species": e["organism"]["scientificName"],
                "taxon_id": e["organism"].get("taxonId", ""),
                "reviewed": int(e.get("entryType", "").startswith("UniProtKB reviewed")),
                "name_match": int(is_named(e)),
                "tier": tier or "none",
                "pmids": ",".join(pmids),
                "length": len(seq),
                "sequence": seq,
            }
        )
    manual = d / "manual_entries.tsv"
    if manual.exists():
        for r in csv.DictReader(open(manual), delimiter="\t"):
            if not r["pmid"].strip():
                raise SystemExit(f"{manual}: every T1 row needs a PMID ({r['accession']})")
            rows.append({**dict.fromkeys(rows[0], ""), **r, "tier": "T1"})
    with open(d / "truth_all.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter

    print("tier counts:", dict(Counter(r["tier"] for r in rows)))


def split(d, seed=20261008):
    clusters = {}
    with _open(d / "clusters.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            clusters[r["id"]] = r["cluster"]
    sp = split_clusters(clusters, seed)
    path = d / "split.tsv"
    if path.exists():
        raise SystemExit(f"{path} exists: the split is made once")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "cluster", "part", "seed"])
        for p, part in sorted(sp.items()):
            w.writerow([p, clusters[p], part, seed])
    from collections import Counter

    print("split:", dict(Counter(sp.values())))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=["fetch", "tiers", "split"])
    ap.add_argument("--dir", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    d = Path(a.dir)
    {"fetch": fetch, "tiers": tiers, "split": split}[a.cmd](d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
