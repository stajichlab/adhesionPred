"""Steps 2-3: query reviewed fungal Swiss-Prot negatives per positive Pfam family; cap per family."""

import collections
import csv
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
from uni import get_all

ROOT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model"
OUT = ROOT + "/family_matched"
CAP = int(os.environ.get("FAM_CAP", "10"))
SEED = 20261010

pos_rows = list(csv.DictReader(open(OUT + "/positive_pfam.tsv"), delimiter="\t"))
pos_acc = {r["accession"] for r in pos_rows}
fam_pos = collections.defaultdict(set)  # pfam -> positives
for r in pos_rows:
    for pf in [x for x in r["pfam"].split(",") if x]:
        fam_pos[pf].add(r["accession"])
cand_acc = {
    r["accession"]
    for r in csv.DictReader(open(ROOT + "/candidates/candidates.tsv"), delimiter="\t")
}
# family names from UniProt are not in the positive table; use the Pfam accession from xref and name lookup below
fams = sorted(fam_pos)
print("families", len(fams), file=sys.stderr)

rng = random.Random(SEED)
selected, per_family, pool_sizes, pool_overlap_cand, pool_overlap_pos = {}, {}, {}, {}, {}
query_log = []
release = open(OUT + "/release.txt").read().strip()
for fam in fams:
    q = (
        f"(taxonomy_id:4751) AND (xref:pfam-{fam}) AND (reviewed:true) "
        f"AND (NOT keyword:KW-0020) AND (NOT cc_allergen:allergic) AND (length:[10 TO 850])"
    )
    rows, rel = get_all(q, "accession,reviewed,organism_name,length,protein_name,xref_pfam")
    query_log.append((fam, q, len(rows)))
    pool = [
        r
        for r in rows
        if r["Entry"] not in pos_acc and r["Entry"] not in cand_acc and r["Entry"] not in selected
    ]
    pool_sizes[fam] = len(rows)
    pool_overlap_cand[fam] = sum(1 for r in rows if r["Entry"] in cand_acc)
    pool_overlap_pos[fam] = sum(1 for r in rows if r["Entry"] in pos_acc)
    pool.sort(key=lambda r: r["Entry"])
    pick = rng.sample(pool, min(CAP, len(pool))) if pool else []
    per_family[fam] = []
    for r in pick:
        selected[r["Entry"]] = {"family": fam, "row": r}
        per_family[fam].append(r["Entry"])
    print(fam, len(rows), len(pool), len(pick), file=sys.stderr)

json.dump(
    {
        "queries": query_log,
        "release": release,
        "cap": CAP,
        "seed": SEED,
        "pool_sizes": pool_sizes,
        "pool_overlap_cand": pool_overlap_cand,
        "pool_overlap_pos": pool_overlap_pos,
        "per_family": per_family,
        "selected": {k: {"family": v["family"], "row": v["row"]} for k, v in selected.items()},
    },
    open(OUT + "/scripts/step2_state.json", "w"),
    indent=1,
)
print("selected", len(selected))
