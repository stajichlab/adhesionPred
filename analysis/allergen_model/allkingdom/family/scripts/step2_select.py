"""Step 2-3: per (family, group) UniProt candidates, exclusions, seeded cap of 10, first-family-wins."""

import csv
import hashlib
import json
import random
import sys
import time
import urllib.parse

sys.path.insert(
    0,
    "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/allkingdom/family/scripts",
)
from common import BASE, FAM, http_get, positives, read_tsv, uniprot_search_all

SEED = 20261010
CAP = 10
TAX = {"plant": "taxonomy_id:33090", "animal": "(taxonomy_id:33208 AND NOT taxonomy_id:9606)"}
BASEQ = (
    "reviewed:true AND xref:pfam-{pf} AND NOT keyword:KW-0020 AND NOT cc_allergen:allergic "
    "AND NOT antigen AND NOT cancer AND NOT taxonomy_id:9606 AND NOT taxonomy_id:10239 "
    "AND length:[10 TO 850] AND {tax}"
)
FIELDS = "accession,length,organism_name,protein_name,reviewed"

pos = positives()
pos_acc = {p["accession"] for p in pos}
pfam_of = {}
for r in read_tsv(FAM + "/raw/positive_pfam.tsv"):
    pfam_of[r["accession"]] = [x for x in r["pfam"].split(",") if x]
kin = read_tsv(BASE + "/allkingdom/kingdom_negatives.tsv")
kin_acc = {r["accession"] for r in kin}

pairs = set()
for p in pos:
    for pf in pfam_of.get(p["accession"], []):
        pairs.add((pf, p["group"]))
pairs = sorted(pairs)
print("family-group pairs", len(pairs), "families", len({f for f, _ in pairs}))

# Step 2: fetch candidates per pair
cand = {}
meta_rows = []
for pf, grp in pairs:
    q = BASEQ.format(pf=pf, tax=TAX[grp])
    rows, total = uniprot_search_all(q, FIELDS)
    cand[(pf, grp)] = rows
    with open(FAM + "/raw/candidates.tsv", "a") as fh:
        pass
    meta_rows.append((pf, grp, total, len(rows)))
    print(pf, grp, "uniprot total", total, "rows", len(rows), flush=True)

# Save raw candidates (before exclusion)
with open(FAM + "/raw/candidates_raw.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["family", "group", "accession", "length", "organism", "protein_name", "reviewed"])
    for (pf, grp), rows in sorted(cand.items()):
        for r in rows:
            w.writerow(
                [pf, grp, r["Entry"], r["Length"], r["Organism"], r["Protein names"], r["Reviewed"]]
            )

# Family names from InterPro
fname = {}
for pf in sorted({f for f, _ in pairs}):
    try:
        body, _ = http_get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}")
        fname[pf] = json.loads(body)["metadata"]["name"]["name"]
    except Exception as exc:  # noqa: BLE001
        fname[pf] = ""
        print("name fetch failed", pf, exc)
    time.sleep(0.2)

# Steps 3: exclusion, sampling
rng = random.Random(SEED)
picked_all = set()
picks = []
stats = []
for pf, grp in pairs:
    rows = cand[(pf, grp)]
    n_raw = len(rows)
    n_kin = sum(1 for r in rows if r["Entry"] in kin_acc)
    n_pos = sum(1 for r in rows if r["Entry"] in pos_acc)
    pool = sorted(
        {
            r["Entry"]: r
            for r in rows
            if r["Entry"] not in pos_acc
            and r["Entry"] not in kin_acc
            and r["Entry"] not in picked_all
        }.values(),
        key=lambda r: r["Entry"],
    )
    n_pool = len(pool)
    chosen = rng.sample(pool, CAP) if n_pool > CAP else pool
    for r in sorted(chosen, key=lambda r: r["Entry"]):
        picked_all.add(r["Entry"])
        picks.append((r, pf, grp))
    stats.append((pf, grp, n_raw, n_pos, n_kin, n_pool, len(chosen)))

# Fetch sequences for picked accessions
acc_list = [r["Entry"] for r, _, _ in picks]
seqs = {}
for i in range(0, len(acc_list), 100):
    b = acc_list[i : i + 100]
    url = "https://rest.uniprot.org/uniprotkb/accessions?" + urllib.parse.urlencode(
        {"accessions": ",".join(b), "format": "fasta"}
    )
    body, _ = http_get(url)
    cur = None
    for line in body.splitlines():
        if line.startswith(">"):
            cur = line.split("|")[1]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    time.sleep(0.2)
seq = {k: "".join(v) for k, v in seqs.items()}
missing = [a for a in acc_list if a not in seq]
print("picked", len(acc_list), "sequences", len(seq), "missing", len(missing))

with open(FAM + "/family_negatives.faa", "w") as fh:
    for r, _pf, _grp in picks:
        s = seq.get(r["Entry"], "")
        fh.write(f">K_{r['Entry']}\n")
        for j in range(0, len(s), 60):
            fh.write(s[j : j + 60] + "\n")

with open(FAM + "/family_negatives.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(
        [
            "accession",
            "group",
            "family",
            "family_name",
            "organism",
            "length",
            "sha256",
            "reviewed",
            "overlap_with_kingdom_negatives",
        ]
    )
    for r, pf, grp in picks:
        s = seq.get(r["Entry"], "")
        w.writerow(
            [
                r["Entry"],
                grp,
                pf,
                fname.get(pf, ""),
                r["Organism"],
                len(s),
                hashlib.sha256(s.encode()).hexdigest(),
                "yes" if r["Reviewed"] == "reviewed" else "no",
                "yes" if r["Entry"] in kin_acc else "no",
            ]
        )

with open(FAM + "/raw/pair_stats.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(
        [
            "family",
            "group",
            "uniprot_rows",
            "in_positives",
            "in_kingdom_negatives",
            "pool_after_exclusion",
            "sampled",
        ]
    )
    for s in stats:
        w.writerow(s)
with open(FAM + "/raw/pairs_meta.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["family", "group", "uniprot_total_results", "rows_fetched"])
    for m in meta_rows:
        w.writerow(m)
print("done; picks", len(picks))
