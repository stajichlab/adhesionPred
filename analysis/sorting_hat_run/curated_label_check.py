#!/usr/bin/env python3
"""Exploratory check: what do the sorting hat calls say about proteins with curated labels?

Descriptive only. This is not a measure of accuracy. The curated table is small, built from
literature, and some of its proteins tuned modules (flagged `tuned_or_homolog`). Output: one TSV with
a row per curated protein found in a run, and a summary of counts per proteome, curated class and
call: all rows; without the flagged proteins; and without them and with strong evidence only
(E1 or E2 for adhesins, N1 or N2 for hard negatives; E3 rows are mostly labels inherited from a shared domain).

Id mapping, cached in --idmap:
  Af293 UniProt   accession inside the protein ID (sp|ACC|NAME)
  S288C           gene name in the FASTA header
  C. albicans     CGD API (gene name -> systematic name), cached
  C. immitis RS   exact sequence match to the UniProt entry of the curated accession
Usage: curated_label_check.py --work WORKDIR --out-prefix PREFIX [--idmap FILE]
"""

import argparse
import collections
import csv
import gzip
import json
import re
import sys
import time
import urllib.request

TUNED_ACC = {"Q4WXJ1", "P41746", "Q8NK60", "Q8NK61", "Q96V71", "A0A0E1RVD3", "Q1E3R8"}
TUNED_GENE = re.compile(
    r"^(ALS\d|FLO\d|SOWGP|PRA|ROD[A-G]$|CAL[ABC]$|CTS1|EPA\d|HWP1|IFF|HYR)", re.I
)
CALLS = [
    "signal_peptide_protein[R0]",
    "tandem_repeat_protein",
    "wall_family_domain",
    "cell_wall_adhesion_candidate[R0]",
]
PROTEOMES = {
    "Afum_Af293_UniProt": "Aspergillus fumigatus",
    "Scer_S288C": "Saccharomyces cerevisiae",
    "Calb_SC5314": "Candida albicans",
    "Cimm_RS": "Coccidioides immitis",
}


def calls_of(path):
    d = collections.defaultdict(dict)
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[r["protein"]][r["call"] + ("[" + r["variant"] + "]" if r["variant"] else "")] = r[
                "value"
            ]
    return d


def fasta_headers(path):
    h = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                h[line[1:].split()[0]] = line[1:].strip()
    return h


def fasta_seqs(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).rstrip("*")
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).rstrip("*")
    return seqs


def cgd_systematic(gene, cache):
    if gene in cache:
        return cache[gene]
    try:
        d = json.load(
            urllib.request.urlopen(f"https://www.candidagenome.org/api/locus/{gene}", timeout=40)
        )
        cache[gene] = d["results"]["Candida albicans SC5314"]["feature_name"]
    except Exception:
        cache[gene] = ""
    time.sleep(0.2)
    return cache[gene]


def uniprot_seq(acc):
    try:
        t = (
            urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=40)
            .read()
            .decode()
        )
        return "".join(t.splitlines()[1:])
    except Exception:
        return ""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--work", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--curated", default="data/curated/adhesins/adhesins.tsv")
    ap.add_argument("--idmap", default=None)
    a = ap.parse_args()
    cur = list(csv.DictReader(open(a.curated, encoding="utf-8-sig"), delimiter="\t"))
    cache = {}
    if a.idmap:
        try:
            cache = json.load(open(a.idmap))
        except Exception:
            cache = {}
    cgd = cache.setdefault("cgd", {})
    found = []
    for prot, org in PROTEOMES.items():
        calls = calls_of(f"{a.work}/{prot}/out/calls.long.tsv.gz")
        rows = [r for r in cur if org in r["organism"]]
        if prot == "Afum_Af293_UniProt":
            ids = {
                re.match(r"^(?:sp|tr)\|([A-Z0-9]+)\|", p).group(1): p
                for p in calls
                if re.match(r"^(?:sp|tr)\|", p)
            }
            pairs = [(r, ids.get(r["accession"])) for r in rows]
        elif prot == "Scer_S288C":
            hd = fasta_headers(f"{a.work}/Scer_S288C.faa")
            g2id = {h.split()[1].upper(): p for p, h in hd.items() if len(h.split()) > 1}
            pairs = [(r, g2id.get(r["gene"].upper()) if r["gene"] else None) for r in rows]
        elif prot == "Calb_SC5314":
            pairs = []
            for r in rows:
                sysname = cgd_systematic(r["gene"], cgd) if r["gene"] else ""
                pairs.append((r, sysname if sysname in calls else None))
        else:
            seqs = fasta_seqs(f"{a.work}/Cimm_RS.faa")
            by_seq = collections.defaultdict(list)
            for p, s in seqs.items():
                by_seq[s].append(p)
            pairs = []
            for r in rows:
                hit = by_seq.get(uniprot_seq(r["accession"]), [])
                pairs.append((r, hit[0] if len(hit) == 1 else None))
        for r, pid in pairs:
            if pid is None or pid not in calls:
                continue
            tuned = r["accession"] in TUNED_ACC or bool(TUNED_GENE.match((r["gene"] or "").strip()))
            found.append(
                {
                    "proteome": prot,
                    "protein": pid,
                    "gene": r["gene"],
                    "accession": r["accession"],
                    "cls": r["cls"],
                    "evidence": r["evidence_level"],
                    "tuned_or_homolog": "yes" if tuned else "no",
                    **{c: calls[pid].get(c, "NA") for c in CALLS},
                }
            )
    if a.idmap:
        json.dump(cache, open(a.idmap, "w"), indent=1)
    with open(a.out_prefix + ".rows.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(found[0]), delimiter="\t")
        w.writeheader()
        w.writerows(found)
    lines = ["proteome\tcurated_class\tset\tn\t" + "\t".join(c + "=called" for c in CALLS)]
    for prot in PROTEOMES:
        for cls in (
            "adhesin",
            "hard_negative",
            "surface_other_adhesion_phenotype",
            "indirect_regulator",
        ):
            for label, keep in (
                ("all", lambda r: True),
                ("excluding_tuned_or_homolog", lambda r: r["tuned_or_homolog"] == "no"),
                (
                    "E1_E2_or_N_and_not_tuned",
                    lambda r: r["tuned_or_homolog"] == "no"
                    and r["evidence"] in ("E1", "E2", "N1", "N2"),
                ),
            ):
                sel = [r for r in found if r["proteome"] == prot and r["cls"] == cls and keep(r)]
                if sel:
                    lines.append(
                        f"{prot}\t{cls}\t{label}\t{len(sel)}\t"
                        + "\t".join(str(sum(r[c] == "called" for r in sel)) for c in CALLS)
                    )
    open(a.out_prefix + ".summary.tsv", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
