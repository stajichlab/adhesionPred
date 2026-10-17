"""Measured counts and overlap checks for counts.md."""

import collections
import sys

sys.path.insert(
    0,
    "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/allkingdom/family/scripts",
)
from common import BASE, FAM, positives, read_tsv

pos = positives()
pfam = {r["accession"]: r["pfam"] for r in read_tsv(FAM + "/raw/positive_pfam.tsv")}
out = {}
pos_g = collections.Counter(p["group"] for p in pos)
nopf = collections.Counter(p["group"] for p in pos if not pfam.get(p["accession"]))
pf_fam = collections.defaultdict(set)
for p in pos:
    for f in [x for x in pfam.get(p["accession"], "").split(",") if x]:
        pf_fam[p["group"]].add(f)
neg = read_tsv(FAM + "/family_negatives.tsv")
neg_g = collections.Counter(r["group"] for r in neg)
neg_fam = collections.defaultdict(set)
for r in neg:
    neg_fam[r["group"]].add(r["family"])
print("pos", pos_g, "nopfam", nopf, "pos families", {g: len(v) for g, v in pf_fam.items()})
print("neg", neg_g, "neg families", {g: len(v) for g, v in neg_fam.items()})
print("pos families both groups", len(pf_fam["plant"] & pf_fam["animal"]))
print("neg families both groups", len(neg_fam["plant"] & neg_fam["animal"]))
# overlaps
kin = {r["accession"] for r in read_tsv(BASE + "/allkingdom/kingdom_negatives.tsv")}
pos_set = {p["accession"] for p in pos}
fung_types = {
    r["uniprot"]
    for r in read_tsv(BASE + "/allergenonline_all/allergenonline_fungi.tsv")
    if r["type"] in ("Aero Fungi", "Food Fungi")
}
fung_all = {r["uniprot"] for r in read_tsv(BASE + "/allergenonline_all/allergenonline_fungi.tsv")}
fpos = {r["accession"] for r in read_tsv(BASE + "/positives/positives.tsv")}
fcand = {r["accession"] for r in read_tsv(BASE + "/candidates/candidates.tsv")}
accs = [r["accession"] for r in neg]
print("neg in kingdom_negatives", sum(a in kin for a in accs))
print("neg in positives(all-kingdom)", sum(a in pos_set for a in accs))
print(
    "neg in fungal allergenonline types",
    sum(a in fung_types for a in accs),
    "fung_all",
    sum(a in fung_all for a in accs),
)
print(
    "neg in fungal positives",
    sum(a in fpos for a in accs),
    "fungal candidates",
    sum(a in fcand for a in accs),
)
print(
    "kingdom_negatives count",
    len(kin),
    "fungal positives",
    len(fpos),
    "fungal candidates",
    len(fcand),
)
