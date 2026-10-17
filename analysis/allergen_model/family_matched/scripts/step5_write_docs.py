"""Write counts.md from measured state files."""

import collections
import csv
import datetime
import json

OUT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/family_matched"
st = json.load(open(OUT + "/scripts/step2_state.json"))
s3 = json.load(open(OUT + "/scripts/step3_state.json"))
pp = list(csv.DictReader(open(OUT + "/positive_pfam.tsv"), delimiter="\t"))
found = [r for r in pp if r["in_uniprot"] == "true"]
with_pf = [r for r in found if r["pfam"]]
fams = sorted(st["pool_sizes"])
neg_per = collections.Counter(v["family"] for v in st["selected"].values())
pos_per = collections.Counter()
for r in found:
    for f in r["pfam"].split(","):
        if f:
            pos_per[f] += 1
clu = list(csv.DictReader(open(OUT + "/cluster/cluster_members.tsv"), delimiter="\t"))
n_pos_c = sum(1 for r in clu if int(r["n_pos"]) > 0)
n_neg_c = sum(1 for r in clu if int(r["n_neg"]) > 0)
n_sh = sum(1 for r in clu if r["class"] == "shared")
today = datetime.date.today().isoformat()
L = []
L.append("# Family-matched negatives: counts\n")
L.append(
    f"Date: {today}. UniProt release {st['release']}. Per-family cap {st['cap']}. Sampling seed {st['seed']}.\n"
)
L.append("## Totals\n")
L.append("| Item | Count |\n|---|---|")
L.append("| Positives in positives.tsv | 200 |")
L.append(f"| Positives with UniProt accession (queried) | {len(found)} |")
L.append(f"| UniProt positives with at least one Pfam domain | {len(with_pf)} |")
L.append(f"| Positives skipped (NCBI protein IDs, no UniProt mapping) | {200 - len(found)} |")
L.append(f"| Distinct Pfam families used by positives | {len(fams)} |")
L.append(f"| Families with at least one new negative | {sum(1 for f in fams if neg_per[f] > 0)} |")
L.append(f"| New negatives (family_negatives.faa) | {len(s3['seqlen'])} |")
L.append(f"| Distinct new negative accessions | {len(set(s3['seqlen']))} |")
L.append("| New negatives in candidates.tsv | 0 (excluded from pools) |")
L.append(
    f"| Pool rows (all families, with duplicates) that are candidates | {sum(st['pool_overlap_cand'].values())} |"
)
L.append(f"| Pool rows that are positives | {sum(st['pool_overlap_pos'].values())} |")
L.append(f"| Pool rows (all families, with duplicates) | {sum(st['pool_sizes'].values())} |")
L.append("")
L.append("## MMseqs2 clusters (easy-cluster, --min-seq-id 0.3 -c 0.5 --cov-mode 0)\n")
L.append("| Item | Count |\n|---|---|")
L.append(
    f"| Input sequences (200 positives + {len(s3['seqlen'])} new negatives) | {200 + len(s3['seqlen'])} |"
)
L.append(f"| Clusters containing positives | {n_pos_c} |")
L.append(f"| Clusters containing new negatives | {n_neg_c} |")
L.append(f"| Clusters containing both | {n_sh} |")
L.append(f"| Positive-only clusters | {n_pos_c - n_sh} |")
L.append(f"| Negative-only clusters | {n_neg_c - n_sh} |")
L.append("")
L.append("## Per family\n")
L.append(
    "Pool = reviewed fungal Swiss-Prot entries with the Pfam, after removing the KW-0020 and cc_allergen:allergic hits in the query; The pool column subtracts candidate and positive hits. Negatives already picked for an earlier family (families in sorted order) are removed before sampling, so New negatives can be below the cap.\n"
)
L.append(
    "| Pfam | Positives with this Pfam | Query hits | Pool (not positive/candidate) | Hits in candidates.tsv | Hits that are positives | New negatives |"
)
L.append("|---|---|---|---|---|---|---|")
for f in fams:
    L.append(
        f"| {f} | {pos_per[f]} | {st['pool_sizes'][f]} | {st['pool_sizes'][f] - st['pool_overlap_cand'][f] - st['pool_overlap_pos'][f]} | {st['pool_overlap_cand'][f]} | {st['pool_overlap_pos'][f]} | {neg_per[f]} |"
    )
L.append("")
L.append(
    "Pfam names: "
    + "; ".join(f"{f} {s3['fam_names'].get(f, '')}" for f in fams if neg_per[f])
    + ".\n"
)
L.append("## Skipped positives\n")
L.append(
    f"{200 - len(found)} positives have NCBI protein IDs (source ncbi_protein or compare, IDs like AAM77471.1 or XP_003030591.1). They were not assigned Pfam domains, so they add no family. The ID mapping step (scripts/step1b_idmap.py) used UniProt idmapping from EMBL-GenBank-DDBJ_CDS and RefSeq_Protein to UniProtKB. It mapped 0 of 26. The EMBL-to-UniParc step then returned HTTP 400, so no other route was tried. The 26 are kept in the clustering input.\n"
)
open(OUT + "/counts.md", "w").write("\n".join(L) + "\n")
print("ok")
