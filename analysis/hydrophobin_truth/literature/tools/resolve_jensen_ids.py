#!/usr/bin/python3.12
"""Resolve the Jensen 2010 gene IDs (protein_lists.tsv) to protein sequences from local FungiDB-31 files,
the local A1163 UniProt proteome, and UniProt for the Af293 AFUA_ IDs. Writes jensen_sequences.faa and
jensen_resolution.tsv (id, species, source file, matched header, length, status). Nothing is guessed:
an ID that is not found is listed as unresolved.
"""

import csv
import json
import re
import urllib.parse
import urllib.request
from collections import Counter

FDB = "/bigdata/operations/pkgadmin/srv/projects/db/FungiDB/release-31/FungiDB-31_{}_AnnotatedProteins.fasta"
SPECIES = {
    "A. oryzae RIB40": "AoryzaeRIB40",
    "A. niger CBS 513.88": "AnigerCBS513-88",
    "A. niger ATCC 1015": "AnigerATCC1015",
    "E. nidulans FGSC A4": "AnidulansFGSCA4",
    "A. flavus NRRL 3357": "AflavusNRRL3357",
    "A. terreus NIH 2624": "AterreusNIH2624",
    "A. clavatus NRRL 1": "AclavatusNRRL1",
}
A1163 = "/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat/Afum_A1163.faa"


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name:
        out.append((name, "".join(buf)))
    return out


def find_fdb(entries, key):
    pat = re.compile(r"(^|[ |=])" + re.escape(key) + r"(-T|-t26_1|\.\d)?(-p\d+)?( |$)")
    hits = [(h, s) for h, s in entries if pat.search(h.split(" | location")[0])]
    return hits


def uniprot_gene(gene, taxon):
    q = f"(gene:{gene}) AND (taxonomy_id:{taxon})"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "fields": "accession,gene_names,organism_name,sequence"}
    )
    with urllib.request.urlopen(url, timeout=60) as r:
        res = json.load(r)["results"]
    return [
        (f"{x['primaryAccession']}|{x['organism']['scientificName']}", x["sequence"]["value"])
        for x in res
    ]


rows = [
    r for r in csv.DictReader(open("protein_lists.tsv"), delimiter="\t") if r["pmid"] == "21182770"
]
cache = {}
res_rows, seqs = [], []
for r in rows:
    pid, sp = r["protein_id"], r["species"]
    hits, src = [], ""
    if sp in SPECIES:
        f = FDB.format(SPECIES[sp])
        cache.setdefault(f, read_fasta(f))
        key = pid
        if sp == "A. niger ATCC 1015" and pid.startswith("JGI"):
            key = "ASPNIDRAFT_" + pid[3:]
        if sp == "E. nidulans FGSC A4":
            key = re.sub(r"\.\d+$", "", pid)
        hits, src = find_fdb(cache[f], key), f.split("/")[-1]
    elif sp == "A. fumigatus A1163":
        cache.setdefault(A1163, read_fasta(A1163))
        hits = [(h, s) for h, s in cache[A1163] if f"GN={pid} " in h]
        src = "Afum_A1163.faa"
    elif sp == "A. fumigatus AF293":
        hits = uniprot_gene(pid, 330879)
        src = "UniProt REST gene:" + pid
    status = "resolved" if len(hits) == 1 else ("ambiguous" if hits else "unresolved")
    if len(hits) == 1:
        h, s = hits[0]
        seqs.append((f"jensen|{pid}", s))
    res_rows.append(
        [
            pid,
            sp,
            src,
            hits[0][0][:90] if hits else "",
            len(hits[0][1]) if len(hits) == 1 else "",
            status,
            len(hits),
        ]
    )
with open("jensen_sequences.faa", "w") as f:
    for k, s in seqs:
        f.write(f">{k}\n{s}\n")
with open("jensen_resolution.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["protein_id", "species", "source", "matched_header", "length", "status", "n_hits"])
    w.writerows(res_rows)
print(Counter(r[5] for r in res_rows))
for r in res_rows:
    if r[5] != "resolved":
        print(r[:3], r[5], r[6])
