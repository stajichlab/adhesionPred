"""Step 4 summary: cluster counts, per-family counts, skipped positives."""

import collections

OUT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/family_matched"
clu = collections.defaultdict(list)
for ln in open(OUT + "/cluster/clu_cluster.tsv"):
    rep, mem = ln.rstrip("\n").split("\t")
    clu[rep].append(mem)
pos_c = neg_c = shared = 0
neg_in_shared = pos_in_shared = 0
with open(OUT + "/cluster/cluster_members.tsv", "w") as f:
    f.write("cluster_rep\tn_pos\tn_neg\tclass\n")
    for rep, mem in clu.items():
        np_ = sum(1 for m in mem if m.startswith("P_"))
        nn = sum(1 for m in mem if m.startswith("F_"))
        if np_ and nn:
            cls = "shared"
            shared += 1
            neg_in_shared += nn
            pos_in_shared += np_
        elif np_:
            cls = "positive_only"
        else:
            cls = "negative_only"
        pos_c += bool(np_)
        neg_c += bool(nn)
        f.write(f"{rep}\t{np_}\t{nn}\t{cls}\n")
pos_ids = [
    line[1:].strip() for line in open(OUT + "/../positives/positives.faa") if line.startswith(">")
]
neg_ids = [line[1:].strip() for line in open(OUT + "/family_negatives.faa") if line.startswith(">")]
# cluster stats
neg_clusters = neg_c
print("positive sequences", len(pos_ids), "new negatives", len(neg_ids))
print("clusters containing positives", pos_c)
print("clusters containing new negatives", neg_c)
print("clusters containing both", shared)
print(
    "new negatives in shared clusters", neg_in_shared, "positives in shared clusters", pos_in_shared
)
print("negative-only clusters", sum(1 for v in clu.values() if all(m.startswith("F_") for m in v)))
print("positive-only clusters", sum(1 for v in clu.values() if all(m.startswith("P_") for m in v)))
