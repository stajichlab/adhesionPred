#!/usr/bin/env python3
"""Widen the repeat-call truth tables with reviewed, secreted UniProt proteins of the two species.

Population: UniProtKB reviewed entries with the keyword Signal (KW-0732) for S. cerevisiae S288C
(559292) and C. albicans SC5314 (237561). Each entry is mapped to the run proteome by an exact
sequence match (an entry whose sequence matches no protein, or several, is skipped). Label: 1 when the
entry has two or more `Repeat` features, 0 when it has none and no Region, Compositional bias or Domain
feature that mentions "repeat". Entries with one feature or a text mention are left out.
NOTE: a negative is an absent annotation, not a proof; `Repeat` features also mark repeat DOMAINS (WD40,
ankyrin) that the repeat detectors are not built to find.
The curated labels of make_candidates.py are kept; evidence of a repeat beats an assumed negative.

Usage: expand_with_uniprot_secreted.py --candidates repeat_truth_candidates.tsv --fasta PROT=PATH ... --out DIR
Needs `mmseqs` on PATH. UniProt is queried live; the release is recorded in the output.
"""

import argparse
import csv
import datetime
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_truth import cluster, read_fasta  # noqa: E402

TAXA = {"Scer_S288C": 559292, "Calb_SC5314": 237561}


def fetch(taxon):
    query = f"organism_id:{taxon} AND reviewed:true AND keyword:KW-0732"
    url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=" + urllib.parse.quote(query)
    with urllib.request.urlopen(url, timeout=300) as resp:
        release = resp.headers.get("X-UniProt-Release", "")
        data = json.load(resp)
    return data["results"], release


def summarize(entry):
    feats = entry.get("features", [])
    n_rep = sum(1 for f in feats if f["type"] == "Repeat")
    mention = sum(
        1
        for f in feats
        if f["type"] in ("Region", "Compositional bias", "Domain", "Motif")
        and re.search("repeat", f.get("description", ""), re.I)
    )
    return n_rep, mention


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--fasta", action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cur = list(csv.DictReader(open(a.candidates), delimiter="\t"))
    out = Path(a.out)
    for spec in a.fasta:
        prot, path = spec.split("=", 1)
        seqs = read_fasta(path)
        by_seq = {}
        for k, s in seqs.items():
            by_seq.setdefault(s, []).append(k)
        entries, release = fetch(TAXA[prot])
        rows, skipped = {}, {"no exact match": 0, "several matches": 0, "ambiguous label": 0}
        for e in entries:
            hits = by_seq.get(e["sequence"]["value"], [])
            if not hits:
                skipped["no exact match"] += 1
                continue
            if len(hits) > 1:
                skipped["several matches"] += 1
                continue
            n_rep, mention = summarize(e)
            if n_rep >= 2:
                lab, why = "1", "UniProt reviewed, secreted, >=2 repeat features"
            elif n_rep == 0 and mention == 0:
                lab, why = "0", "UniProt reviewed, secreted, no repeat feature (assumed negative)"
            else:
                skipped["ambiguous label"] += 1
                continue
            rows[hits[0]] = {
                "protein": hits[0],
                "gene": e.get("genes", [{}])[0].get("geneName", {}).get("value", ""),
                "accession": e["primaryAccession"],
                "label": lab,
                "basis": why,
                "curated_class": "uniprot_secreted",
                "evidence": "-",
                "tuned_or_homolog": "unknown",
            }
        n_uni = len(rows)
        for r in cur:  # curated labels are kept; a repeat claim beats an assumed negative
            if r["proteome"] != prot or r["label"] not in ("0", "1"):
                continue
            old = rows.get(r["protein"])
            if old and old["label"] == "1" and r["label"] == "0":
                continue
            rows[r["protein"]] = {
                k: r[k]
                for k in (
                    "protein",
                    "gene",
                    "accession",
                    "label",
                    "basis",
                    "curated_class",
                    "evidence",
                    "tuned_or_homolog",
                )
            }
        keep = list(rows.values())
        rep = cluster({r["protein"]: seqs[r["protein"]] for r in keep})
        with open(out / f"truth_v2.{prot}.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["id", "label", "cluster"])
            for r in keep:
                w.writerow([r["protein"], r["label"], rep[r["protein"]]])
        with open(out / f"truth_v2.{prot}.annotated.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(
                [
                    "id",
                    "gene",
                    "accession",
                    "label",
                    "cluster",
                    "basis",
                    "curated_class",
                    "evidence",
                    "tuned_or_homolog",
                ]
            )
            for r in keep:
                w.writerow(
                    [
                        r["protein"],
                        r["gene"],
                        r["accession"],
                        r["label"],
                        rep[r["protein"]],
                        r["basis"],
                        r["curated_class"],
                        r["evidence"],
                        r["tuned_or_homolog"],
                    ]
                )
        pos = [r for r in keep if r["label"] == "1"]
        neg = [r for r in keep if r["label"] == "0"]
        print(
            f"{prot}: UniProt release {release} fetched {datetime.date.today()}: {len(entries)} entries, "
            f"{n_uni} usable, skipped {skipped}; truth_v2 {len(pos)} positives in "
            f"{len({rep[r['protein']] for r in pos})} clusters, {len(neg)} negatives in "
            f"{len({rep[r['protein']] for r in neg})} clusters",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
