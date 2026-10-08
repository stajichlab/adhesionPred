#!/usr/bin/env python3
"""Truth set for the repeat call under the owner's definition of 2026-10-08.

`tandem_repeat_protein` means a repeating motif or array as in FLO11. A repeating DOMAIN counts only if it
is a known domain associated with adhesion (owner's examples: Ser/Thr-rich tandem repeats, the Thr-rich
functional amyloid core, Hwp1 repeats, Iff/Hyr repeats).

Labels from the UniProt `Repeat` features of each entry (a feature description is a number such as "1" or
"2-5", a family name such as "ALS 3", or a domain name such as "WD 4"):
  array feature   : unnamed/numbered, or a name on ADHESION_FAMILIES (ALS, PIR, HYR, IFF, HWP, FLO, ...)
  domain feature  : a name on GLOBULAR_DOMAINS (WD, LRR, BNR, Sel1, PbH, CXXCXGXG, ANK, TPR, ...)
  other named     : any other name (kept out of both classes)
Positive = at least 2 array features. Negative = no Repeat feature at all and no Region, Compositional bias
or Domain feature that mentions "repeat" (an ASSUMED negative: an absent annotation). A protein whose only
repeats are globular domains, or whose names are unknown, is left out: it is neither positive nor negative.
Paper statements from the adjudicated curation table are positive. Evidence of a repeat beats an assumed
negative when one protein has several rows.

Population: reviewed UniProt entries with a signal peptide (KW-0732) plus the curated adhesin-table proteins,
mapped to the run proteome by exact sequence. Clusters: MMseqs2, 30% identity, coverage 0.5.

Usage: build_truth_v3.py --candidates repeat_truth_candidates.tsv --fasta PROT=PATH ... --out DIR
Needs `mmseqs` on PATH. UniProt is queried live.
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
ADHESION_FAMILIES = re.compile(
    r"^(ALS|PIR|HYR|IFF|HWP|HPF|FLO|EPA|AWP|SRP|GLEYA|FLOCCULIN|CWP|PGA)", re.I
)
GLOBULAR_DOMAINS = re.compile(
    r"^(WD|LRR|BNR|SEL1?|PBH\d?|CXXCXGXG|ANK|TPR|HEAT|ARM|KELCH|RCC\d?|EF-HAND|PPR|NHL|PQQ|RLD|TIM|CBS|ZF|ZINC)\b",
    re.I,
)


def kind(description):
    """'array' | 'domain' | 'other' for one Repeat feature description."""
    d = re.sub(r";.*$", "", description.strip())  # drop '; approximate' and '; truncated'
    name = re.sub(r"[\d\s\-]+$", "", d).strip()
    if not name:
        return "array"  # numbered only, as in "1", "2-5"
    if GLOBULAR_DOMAINS.match(name):
        return "domain"
    if ADHESION_FAMILIES.match(name):
        return "array"
    return "other"


def fetch_by_query(query):
    url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=" + urllib.parse.quote(query)
    with urllib.request.urlopen(url, timeout=300) as resp:
        release = resp.headers.get("X-UniProt-Release", "")
        return json.load(resp)["results"], release


def fetch_accessions(accs, batch=40):
    out = []
    for i in range(0, len(accs), batch):
        q = " OR ".join(f"accession:{a}" for a in accs[i : i + batch])
        out += fetch_by_query(q)[0]
    return out


def summarize(entry):
    feats = entry.get("features", [])
    kinds = [kind(f.get("description", "")) for f in feats if f["type"] == "Repeat"]
    mention = sum(
        1
        for f in feats
        if f["type"] in ("Region", "Compositional bias", "Domain", "Motif")
        and re.search("repeat", f.get("description", ""), re.I)
    )
    return kinds.count("array"), kinds.count("domain"), kinds.count("other"), mention


def label_of(n_arr, n_dom, n_oth, mention):
    if n_arr >= 2:
        return "1", "UniProt >=2 array-type repeat features"
    if n_arr == 0 and n_dom == 0 and n_oth == 0 and mention == 0:
        return "0", "no UniProt repeat feature (assumed negative)"
    if n_arr == 0 and n_dom > 0 and n_oth == 0:
        return "", "only globular repeat domains (excluded)"
    return "", "ambiguous repeat annotation (excluded)"


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
        entries, release = fetch_by_query(
            f"organism_id:{TAXA[prot]} AND reviewed:true AND keyword:KW-0732"
        )
        accs = sorted({r["accession"] for r in cur if r["proteome"] == prot and r["accession"]})
        have = {e["primaryAccession"] for e in entries}
        entries += fetch_accessions([x for x in accs if x not in have])
        rows, skipped = {}, {"no exact match": 0, "several matches": 0, "left out": 0}
        why_out = {}
        for e in entries:
            hits = by_seq.get(e["sequence"]["value"], [])
            if not hits:
                skipped["no exact match"] += 1
                continue
            if len(hits) > 1:
                skipped["several matches"] += 1
                continue
            lab, why = label_of(*summarize(e))
            if lab == "":
                skipped["left out"] += 1
                why_out[why] = why_out.get(why, 0) + 1
                continue
            gene = (e.get("genes") or [{}])[0].get("geneName", {}).get("value", "")
            row = {
                "protein": hits[0],
                "gene": gene,
                "accession": e["primaryAccession"],
                "label": lab,
                "basis": why,
                "curated_class": "uniprot",
                "evidence": "-",
                "tuned_or_homolog": "unknown",
            }
            old = rows.get(hits[0])
            if old is None or (lab == "1" and old["label"] == "0"):
                rows[hits[0]] = row
        n_uni = len(rows)
        by_acc = {e["primaryAccession"]: e for e in entries}
        for r in cur:  # curated proteins keep the protein ID of the curated table (their sequence can differ slightly from UniProt)
            if r["proteome"] != prot or not r["accession"]:
                continue
            e = by_acc.get(r["accession"])
            base = {
                k: r[k]
                for k in (
                    "protein",
                    "gene",
                    "accession",
                    "curated_class",
                    "evidence",
                    "tuned_or_homolog",
                )
            }
            if r["basis"] == "paper statement" and r["label"] == "1":
                new = {**base, "label": "1", "basis": "paper statement"}
            elif e is not None:
                lab, why = label_of(*summarize(e))
                if lab == "":
                    rows.pop(r["protein"], None)
                    continue
                new = {**base, "label": lab, "basis": why}
            else:
                continue
            old = rows.get(r["protein"])
            if old is None or new["label"] == "1" or old["label"] != "1":
                rows[r["protein"]] = new
        keep = list(rows.values())
        rep = cluster({r["protein"]: seqs[r["protein"]] for r in keep})
        with open(out / f"truth_v3.{prot}.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["id", "label", "cluster"])
            for r in keep:
                w.writerow([r["protein"], r["label"], rep[r["protein"]]])
        with open(out / f"truth_v3.{prot}.annotated.tsv", "w", newline="") as f:
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
                        r["id"] if "id" in r else r["protein"],
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
            f"{prot}: UniProt release {release} fetched {datetime.date.today()}: {len(entries)} entries, {n_uni} labelled, skipped {skipped} {why_out}; "
            f"truth_v3 {len(pos)} positives in {len({rep[r['protein']] for r in pos})} clusters, "
            f"{len(neg)} negatives in {len({rep[r['protein']] for r in neg})} clusters",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
