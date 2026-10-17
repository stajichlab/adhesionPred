#!/usr/bin/env python3.12
"""Non-allergen negatives for the all-kingdom positives (owner decision 2026-10-10, option 1).

Positives: allergenonline/allergenonline_all (AllergenOnline non-fungal rows mapped to UniProt, 1,679 proteins),
  grouped as plant (AllergenOnline types Aero Plant, Food Plant, Gliadin) or animal (all other types).
Negatives: reviewed Swiss-Prot entries, same filters as the fungal weak negatives (no KW-0020, no allergen comment,
  no antigen word, no Cancer, length 10-850, protein-level evidence, annotation score 4-5), excluding human (9606).
  Plant pool: taxonomy_id 33090 (Viridiplantae). Animal pool: taxonomy_id 33208 (Metazoa), excluding 9606.
Sampling: within each group, length-matched 1:1 per 50-residue bin (seed 20261010), excluding every positive accession.

Outputs (in --outdir): kingdom_negatives.faa (IDs K_<accession>), kingdom_negatives.tsv (accession, group, length, bin,
  organism, sha256), counts.md (per group: positives, pool size, sampled, shortfall bins).
"""

import argparse
import collections
import csv
import datetime
import hashlib
import os
import random
import urllib.parse
import urllib.request

BIN = 50
SEED = 20261010
BASE = (
    "reviewed:true AND NOT keyword:KW-0020 AND NOT cc_allergen:allergic AND NOT antigen AND NOT Cancer "
    "AND length:[10 TO 850] AND existence:1 AND (annotation_score:4 OR annotation_score:5)"
)
POOLS = {"plant": "taxonomy_id:33090", "animal": "taxonomy_id:33208 AND NOT taxonomy_id:9606"}
PLANT_TYPES = {"Aero Plant", "Food Plant", "Gliadin"}


def stream(query):
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {
            "query": f"{query} AND {BASE}",
            "format": "tsv",
            "fields": "accession,organism_name,length,sequence",
        }
    )
    req = urllib.request.Request(url, headers={"User-Agent": "adhesionPred-kingdom-negatives/1.0"})
    with urllib.request.urlopen(req, timeout=900) as resp:
        return resp.read().decode("utf-8").rstrip("\n").split("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--allergenonline-fungi-tsv",
        required=True,
        help="allergenonline_all/allergenonline_fungi.tsv",
    )
    ap.add_argument(
        "--positives-tsv",
        required=True,
        help="positives/positives.tsv (fungal set, excluded from negatives)",
    )
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(SEED)

    group_of = {}
    for r in csv.DictReader(open(args.allergenonline_fungi_tsv), delimiter="\t"):
        if r["new"] != "yes" or not r["uniprot"]:
            continue
        group_of[r["uniprot"]] = "plant" if r["type"] in PLANT_TYPES else "animal"
    pos_len = collections.defaultdict(list)
    pos_acc = set(group_of)
    for r in csv.DictReader(open(args.positives_tsv), delimiter="\t"):
        pos_acc.add(r["accession"])
    # positive lengths per group, from the all-kingdom FASTA's TSV
    lenmap = {}
    for r in csv.DictReader(
        open(
            os.path.join(os.path.dirname(args.allergenonline_fungi_tsv), "allergenonline_new.tsv")
        ),
        delimiter="\t",
    ):
        lenmap[r["accession"]] = int(r["length"])
    for a, g in group_of.items():
        if a in lenmap:
            pos_len[g].append(lenmap[a])
    print({g: len(v) for g, v in pos_len.items()})

    summary = {}
    chosen = []
    for g, query in POOLS.items():
        lines = stream(query)
        pool = collections.defaultdict(list)
        for ln in lines[1:]:
            f = ln.split("\t")
            acc, org, L, seq = f[0], f[1], int(f[2]), f[3].strip()
            if acc in pos_acc or not seq:
                continue
            pool[L // BIN].append((acc, org, L, seq))
        want = collections.Counter(L // BIN for L in pos_len[g])
        shortfall = {}
        for b in sorted(want):
            cands = sorted(pool.get(b, []))
            k = want[b]
            take = min(k, len(cands))
            if take < k:
                shortfall[b] = (k, len(cands))
            chosen += [(g, *c) for c in rng.sample(cands, take)] if take else []
        summary[g] = {
            "positives": len(pos_len[g]),
            "pool_after_exclusion": sum(len(v) for v in pool.values()),
            "sampled": sum(1 for c in chosen if c[0] == g),
            "shortfall_bins": shortfall,
        }

    with (
        open(os.path.join(args.outdir, "kingdom_negatives.faa"), "w") as fa,
        open(os.path.join(args.outdir, "kingdom_negatives.tsv"), "w") as tsv,
    ):
        tsv.write("accession\tgroup\tlength\tlength_bin\torganism\tsequence_sha256\n")
        for g, acc, org, L, seq in chosen:
            fa.write(f">K_{acc}\n")
            for i in range(0, len(seq), 60):
                fa.write(seq[i : i + 60] + "\n")
            tsv.write(
                f"{acc}\t{g}\t{L}\t{L // BIN}\t{org}\t{hashlib.sha256(seq.encode()).hexdigest()}\n"
            )

    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "counts.md"), "w") as md:
        md.write("# Non-allergen negatives for the all-kingdom positives\n\n")
        md.write(
            f"- Run date: {today}; UniProt release 2026_03; seed {SEED}; length bin {BIN} aa; 1:1 per bin per group\n\n"
        )
        md.write(
            "| Group | Positives | Pool after exclusion | Sampled | Bins short of candidates |\n|---|---|---|---|---|\n"
        )
        for g, s in summary.items():
            md.write(
                f"| {g} | {s['positives']} | {s['pool_after_exclusion']} | {s['sampled']} | {s['shortfall_bins'] or 'none'} |\n"
            )
    print(summary)


if __name__ == "__main__":
    main()
