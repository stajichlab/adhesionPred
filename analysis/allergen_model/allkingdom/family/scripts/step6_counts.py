"""Write counts.md from measured outputs."""

import collections
import sys

sys.path.insert(
    0,
    "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/allkingdom/family/scripts",
)
from common import FAM, PLANT_TYPES, positives, read_tsv

pos = positives()
pfam = {
    r["accession"]: [x for x in r["pfam"].split(",") if x]
    for r in read_tsv(FAM + "/raw/positive_pfam.tsv")
}
neg = read_tsv(FAM + "/family_negatives.tsv")
stats = read_tsv(FAM + "/raw/pair_stats.tsv")
fname = {r["family"]: r["family_name"] for r in neg}
L = []
L.append("# Family-matched non-allergen negatives (all-kingdom allergens)\n")
L.append(
    "- Date: 2026-10-10; UniProt release 2026_03 (relnotes.txt, 02-Sep-2026); Swiss-Prot only (reviewed:true)"
)
L.append(
    "- Seed 20261010; cap 10 per family per group; MMseqs2 easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5\n"
)
L.append("## Positives\n")
L.append("| Group | Positives | With Pfam | Without Pfam | Pfam families |")
L.append("|---|---|---|---|---|")
pf_fam = collections.defaultdict(set)
for g in ("plant", "animal"):
    ps = [p for p in pos if p["group"] == g]
    nopf = [p for p in ps if not pfam.get(p["accession"])]
    for p in ps:
        for f in pfam.get(p["accession"], []):
            pf_fam[g].add(f)
    L.append(f"| {g} | {len(ps)} | {len(ps)-len(nopf)} | {len(nopf)} | {len(pf_fam[g])} |")
tot = len(pos)
totnopf = sum(1 for p in pos if not pfam.get(p["accession"]))
L.append(
    f"| total | {tot} | {tot-totnopf} | {totnopf} | {len(pf_fam['plant'] | pf_fam['animal'])} |\n"
)
both = pf_fam["plant"] & pf_fam["animal"]
L.append(f"Families used by both groups: {len(both)} (sampled separately per group).")
multi = [p for p in pos if len(set(p["type"].split(";"))) > 1]
L.append(
    f"Positives with more than one AllergenOnline type: {len(multi)}. Group rule: plant if any type is in {sorted(PLANT_TYPES)}; otherwise animal."
)
bact = [
    p
    for p in pos
    if p["group"] == "animal"
    and set(p["type"].split(";")) <= {"Bacteria skin", "Bacteria airway", "Protozoan", "Vaccine"}
]
L.append(
    f"Animal-group positives whose types are only Bacteria/Protozoan/Vaccine: {len(bact)} (sampled from Metazoa per the rule).\n"
)
L.append("## Negatives\n")
L.append("| Group | Negatives | Families with >=1 negative |")
L.append("|---|---|---|")
ng = collections.Counter(r["group"] for r in neg)
nf = collections.defaultdict(set)
for r in neg:
    nf[r["group"]].add(r["family"])
for g in ("plant", "animal"):
    L.append(f"| {g} | {ng[g]} | {len(nf[g])} |")
L.append(f"| total | {len(neg)} | {len({r['family'] for r in neg})} |\n")
L.append(f"Families with negatives in both groups: {len(nf['plant'] & nf['animal'])}.\n")
L.append("## Candidate pools\n")
L.append(
    "| Group | Family-group pairs | Pairs with empty pool | Pairs with 10 sampled | Pairs with <10 sampled |"
)
L.append("|---|---|---|---|---|")
for g in ("plant", "animal"):
    ss = [s for s in stats if s["group"] == g]
    L.append(
        f"| {g} | {len(ss)} | {sum(1 for s in ss if int(s['pool_after_exclusion'])==0)} | {sum(1 for s in ss if int(s['sampled'])==10)} | {sum(1 for s in ss if int(s['sampled'])<10)} |"
    )
L.append("")
L.append(
    "Candidate hits removed by exclusion (summed over pairs, a protein can count in several pairs): "
    f"in positives {sum(int(s['in_positives']) for s in stats)}, in kingdom negatives {sum(int(s['in_kingdom_negatives']) for s in stats)}.\n"
)
L.append("## Clusters (MMseqs2, positives and new negatives together)\n")
cl = collections.defaultdict(list)
for line in open(FAM + "/cluster/clu_cluster.tsv"):
    r, m = line.rstrip("\n").split("\t")[:2]
    cl[r].append(m)
po = no = mx = 0
for ms in cl.values():
    p = sum(1 for m in ms if m.startswith("P_"))
    n = len(ms) - p
    if p and n:
        mx += 1
    elif p:
        po += 1
    else:
        no += 1
L.append("| Cluster type | Count |")
L.append("|---|---|")
L.append(
    f"| positives only | {po} |\n| negatives only | {no} |\n| mixed | {mx} |\n| total | {len(cl)} |\n"
)
L.append("## Per-family counts\n")
L.append("| Group | Family | Family name | Pool after exclusion | Sampled |")
L.append("|---|---|---|---|---|")
for s in sorted(stats, key=lambda s: (s["family"], s["group"])):
    L.append(
        f"| {s['group']} | {s['family']} | {fname.get(s['family'], '')} | {s['pool_after_exclusion']} | {s['sampled']} |"
    )
open(FAM + "/counts.md", "w").write("\n".join(L) + "\n")
print("ok")
