"""Build candidates.tsv and stats.txt from the saved API outputs (ids and metadata only, no sequences).

Usage: python3 -I build_candidates.py WORKDIR MEMBERS_TSV OUTDIR
WORKDIR holds pan.tsv (fetch_pan.py), pan_map.tsv (fetch_idmap.py), pan_tax.json (fetch_tax.py).
"""

import collections
import csv
import json
import re
import sys
from pathlib import Path

work, members, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
with members.open(newline="") as fh:
    mem = {r["hit_id"] for r in csv.DictReader(fh, delimiter="\t")}
mp = dict(
    line.rstrip("\n").split("\t")
    for line in (work / "pan_map.tsv").open()
    if not line.startswith("accession")
)
tax = json.load((work / "pan_tax.json").open())
with (work / "pan.tsv").open(newline="") as fh:
    panther = list(csv.DictReader(fh, delimiter="\t"))
named = re.compile(
    r"ecm33|sps2|sps22|sporulation-specific protein 2|pst1|gpi-anchored cell wall organization|cell surface gpi-anchored",
    re.I,
)
cols = [
    "accession",
    "organism",
    "taxid",
    "clade_fine",
    "source",
    "evidence",
    "tier",
    "panther_subfamily",
    "panther_evalue",
    "uniref50",
    "overlap_with_training",
]
rows = []
for x in panther:
    fine = tax[x["taxid"]][0] or "(no fine clade: %s)" % tax[x["taxid"]][1]
    u = mp.get(x["accession"], "")
    ov = "no_uniref50" if not u else ("in_training_cluster" if u in mem else "new_cluster")
    if x["reviewed"] == "reviewed":
        tier = "T1_swissprot"
    elif named.search(x["name"]):
        tier = "T2_named"
    else:
        tier = "T3_panther_only"
    ev = f"PANTHER {x['model']};uniprot={x['reviewed']};name={x['name']}"
    rows.append(
        [
            x["accession"],
            x["organism"],
            x["taxid"],
            fine,
            "InterPro API PANTHER PTHR31018 (IPR051648)",
            ev,
            tier,
            x["model"].split(":")[-1],
            x["score"],
            u,
            ov,
        ]
    )
with (out / "candidates.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    w.writerows(rows)

L = []
R = [dict(zip(cols, r, strict=False)) for r in rows]
L.append("proteins %d distinct accessions %d" % (len(R), len({r["accession"] for r in R})))
L.append("overlap %s" % dict(collections.Counter(r["overlap_with_training"] for r in R)))
L.append("tier %s" % dict(collections.Counter(r["tier"] for r in R)))
L.append("subfamily %s" % dict(collections.Counter(r["panther_subfamily"] for r in R)))
L.append(
    "distinct uniref50 in training %d new %d"
    % (
        len({r["uniref50"] for r in R if r["overlap_with_training"] == "in_training_cluster"}),
        len({r["uniref50"] for r in R if r["overlap_with_training"] == "new_cluster"}),
    )
)
L.append("clade\tproteins\tin_training\tnew_proteins\tnew_clusters\tnew_clusters_E<=1e-10")
clades = collections.Counter(r["clade_fine"] for r in R)
for c, n in clades.most_common():
    sel = [r for r in R if r["clade_fine"] == c]
    inn = sum(r["overlap_with_training"] == "in_training_cluster" for r in sel)
    nw = [r for r in sel if r["overlap_with_training"] == "new_cluster"]
    strong = {r["uniref50"] for r in nw if float(r["panther_evalue"]) <= 1e-10}
    L.append(
        "%s\t%d\t%d\t%d\t%d\t%d"
        % (c, n, inn, len(nw), len({r["uniref50"] for r in nw}), len(strong))
    )
for k in ("in_training_cluster", "new_cluster", "no_uniref50"):
    s = sorted(float(r["panther_evalue"]) for r in R if r["overlap_with_training"] == k)
    if not s:
        continue
    L.append(
        "%s n=%d median_E=%g E<=1e-10:%d E>1e-5:%d"
        % (k, len(s), s[len(s) // 2], sum(v <= 1e-10 for v in s), sum(v > 1e-5 for v in s))
    )
(out / "stats.txt").write_text("\n".join(L) + "\n")
print("\n".join(L))
