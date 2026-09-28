import csv
import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis._common.paths import FUNGI5K_DUCKDB  # noqa: E402

DB = str(FUNGI5K_DUCKDB)
RS = "FA2214EC"
con = duckdb.connect(DB, read_only=True)
con.execute("PRAGMA memory_limit='3GB'")
best = {}
for line in open("cocci_antigens/idmap.m8"):
    f = line.rstrip("\n").split("\t")
    q, t, pid = f[0], f[1], float(f[2])
    if q not in best or pid > best[q][1]:
        best[q] = (t, pid)
print(f"id map: {len(best)} of ~9139 reference proteins have a Fungi_5k match")
sp = set(
    con.execute(
        f"SELECT DISTINCT replace(protein_id,'.protein','') p FROM signalp WHERE protein_id LIKE '{RS}%'"
    )
    .fetchdf()
    .p
)
tm = set(
    con.execute(
        f"SELECT DISTINCT replace(protein_id,'.protein','') p FROM tmhmm WHERE protein_id LIKE '{RS}%'"
    )
    .fetchdf()
    .p
)
print(f"Fungi_5k RS: {len(sp)} proteins with signal peptide, {len(tm)} with TM helix")
# sanity: the anchors
for name, prot in [
    ("SOWgp", "CIMG_04613-t26_1-p1"),
    ("Ag2/PRA", "CIMG_09696-t26_1-p1"),
    ("PRA3", "CIMG_02492-t26_1-p1"),
    ("CF antigen", "CIMG_02795-t26_1-p1"),
]:
    fid = best.get(prot, ("", 0))[0]
    print(f"   {name:<12}{prot:<22}fungi5k={fid or 'NO MATCH':<22}SP={fid in sp}  TM={fid in tm}")
# rebuild shortlist WITH secretion required
rank = [
    x
    for x in csv.DictReader(open("cocci_antigens/cocci_antigen_ranking.tsv"), delimiter="\t")
    if x["is_orthogroup_representative"] == "yes"
]


def f(x, k):
    return float(x[k])


mapped = [x for x in rank if x["protein"] in best]
secreted = [x for x in mapped if best[x["protein"]][0] in sp]
print(
    f"\nof {len(rank)} orthogroups: {len(mapped)} have annotation, {len(secreted)} have a signal peptide"
)
sl = [
    x
    for x in secreted
    if f(x, "prevalence") >= 0.95
    and f(x, "max_fungal_crossreact_pid") == 0
    and f(x, "human_homolog_pid") == 0
]
print(f"secreted + universal + Coccidioides-specific: {len(sl)}")
with open("cocci_antigens/shortlist_secreted.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(sl[0]) + ["fungi5k", "has_TM"], delimiter="\t")
    w.writeheader()
    for x in sl:
        w.writerow({**x, "fungi5k": best[x["protein"]][0], "has_TM": best[x["protein"]][0] in tm})
print("wrote cocci_antigens/shortlist_secreted.tsv")
