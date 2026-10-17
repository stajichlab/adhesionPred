#!/usr/bin/env python3.12
"""Select the `assumed` (weak) negatives: fungal Swiss-Prot reviewed proteins with no lectin signal.

Filters (all from UniProt fields, none from the hmmsearch scan of this work):
  reviewed (Swiss-Prot), not a fragment, length >= 50, no keyword Lectin, no UniProt Pfam
  cross-reference to any of the 20 models (15 M1/M2 and 5 M3), no word of the lectin screen in
  the protein name or function text, and not one of the accessions in --exclude.
Then a random sample of N rows, with a fixed seed, with at most --max-per-species rows per species.
These rows have NO reason that the protein is not a lectin. They are `assumed`.
Usage: select_assumed.py SWISSPROT_TSV_GZ EXCLUDE_ACC_LIST OUT_ACC_LIST [--n 500] [--seed 20261011]
"""

import argparse
import csv
import gzip
import random
import re

PF = {
    "PF07938",
    "PF01476",
    "PF08881",
    "PF07367",
    "PF21691",
    "PF00337",
    "PF21595",
    "PF18647",
    "PF22803",
    "PF01419",
    "PF00652",
    "PF14200",
    "PF13385",
    "PF01607",
    "PF00187",
    "PF03388",
    "PF00262",
    "PF07915",
    "PF02157",
    "PF11721",
}
WORDS = re.compile(
    r"(lectin|agglutinin|galectin|lysm|chitin[- ]binding|carbohydrate[- ]binding|glycan[- ]binding|"
    r"hemagglutin|sugar[- ]binding|mannose[- ]binding|fucose|calnexin|calreticulin|ricin|jacalin|"
    r"cyanovirin|fruit ?body|malectin|OS-9|EMP4[67])",
    re.I,
)
ap = argparse.ArgumentParser()
ap.add_argument("sp")
ap.add_argument("exclude")
ap.add_argument("out")
ap.add_argument("--n", type=int, default=500)
ap.add_argument("--seed", type=int, default=20261011)
ap.add_argument("--max-per-species", type=int, default=15)
a = ap.parse_args()
excl = {x.strip() for x in open(a.exclude) if x.strip()}
pool = []
for r in csv.DictReader(gzip.open(a.sp, "rt"), delimiter="\t"):
    if r["Entry"] in excl:
        continue
    if r["Fragment"].strip():
        continue
    if int(r["Length"]) < 50:
        continue
    if "Lectin" in r["Keywords"]:
        continue
    if {x for x in r["Pfam"].split(";") if x} & PF:
        continue
    if WORDS.search(r["Protein names"]) or WORDS.search(r["Function [CC]"]):
        continue
    pool.append(r)
print("eligible", len(pool))
rnd = random.Random(a.seed)
rnd.shuffle(pool)
per = {}
out = []
for r in pool:
    sp = r["Organism (ID)"]
    if per.get(sp, 0) >= a.max_per_species:
        continue
    per[sp] = per.get(sp, 0) + 1
    out.append(r["Entry"])
    if len(out) >= a.n:
        break
open(a.out, "w").write("\n".join(out) + "\n")
print("selected", len(out), "species", len(per))
