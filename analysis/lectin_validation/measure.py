#!/usr/bin/env python3.12
"""Measure recall and false-hit rates of the lectin-like Pfam models on the validation truth set.

Inputs: truth_final.tsv, scan.hits.tsv (hmmsearch --cut_ga, 20 models). Optional: miss_subga.tsv and
miss_pfamscan.tsv (made by scan_misses.py) for the miss analysis.
Outputs (--out-dir): counts.tsv, per_row_hits.tsv, recall_by_group.tsv, false_hits.tsv, floors.tsv,
family_recall.tsv, leakage_split.tsv, misses.tsv, species_clusters.tsv.
It makes no call, no status and no threshold for the tool.
"""

import argparse
import csv
import os
from collections import defaultdict

import lectin_metrics as lm

M12 = [
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
]
M3 = ["PF03388", "PF00262", "PF07915", "PF02157", "PF11721"]
REV1 = ["PF07938", "PF21691", "PF00337", "PF07367", "PF21595"]  # five revision-1 models (spec F10)
FLOORS = [None, 50.0, 100.0, 200.0]
GA_PF07938 = 25.3

# A priori family of each positive: UniLectin3D class (and a few named accessions) -> expected model.
# Source: spec F11 class-to-model mapping (made from class names). "none" = the spec found no Pfam model.
CLASS_TO_MODEL = {
    "AAL-like (PropLec6A)": "PF07938",
    "Lyophyllum ginkbilobin-like": "PF21691",
    "Fungal fruit body lectin": "PF07367",
    "galectin-like": "PF00337",
    "Coprinus trefoil lectin-like": "PF21595",
    "cyanovirin-like": "PF08881",
    "LysM-like": "PF01476",
    "Ricin-like": "PF14200",
    "Sclerotinia trefoil lectin-like": "PF14200",
    "CBM67-like": "none",
    "CML-like": "none",
    "Monocot-lectin like": "none",
    "Fungal prism lectins": "none",
    # Review 1 (M4): ricin-B-type beta-trefoil lectins. Mapped to the RicinB_lectin_2 family. Five of the six
    # rows score 21.2 to 24.9 against the GA of 25.6. They are misses of an existing model, not of a missing model.
    "Laetiporus trefoil lectin-like": "PF14200",
    "Ig-like": "none",
    "PVL-like (PropLec7B)": "none",
    "Tectonin (PropLec6B)": "none",
}
ACC_TO_MODEL = {  # literature rows without a UniLectin3D class (expected model by protein family named in the paper)
    "F2ST24": "PF01476",
    "F2SIK4": "PF01476",
    "F9XMT4": "PF01476",
    "G2X4U8": "PF01476",
    "G4N906": "PF01476",
    "H1UVH3": "PF01476",
    "H1VSQ9": "PF01476",
    "J4URT3": "PF01476",
    "J5J7Y8": "PF01476",
    "Q00363": "PF01607",
    "B2ZRS9": "PF14200",
    "Q06100": "PF00337",
    "C1G842": "none",
    "P9WEM2": "PF01476",
    "G1X5Y2": "PF07938",
    "G3JGU5": "PF00652",
}
PRIMARY = [
    "structure_backed",
    "host_binding_h1a",
    "chitin_masking_h1b",
    "literature_binding",
    "known_gap",
]


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def expected_model(r):
    if r["accession"] in ACC_TO_MODEL:
        return ACC_TO_MODEL[r["accession"]]
    return CLASS_TO_MODEL.get(r["unilectin_class"], "unmapped")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--miss-dir", default=None)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    truth = read_tsv(a.truth)
    hit_rows = [(r["accession"], r["model"], float(r["score"])) for r in read_tsv(a.hits)]
    by_class = defaultdict(list)
    for r in truth:
        by_class[r["label_class"]].append(r)

    def members(rows):
        return [(r["accession"], r["cluster30"]) for r in rows]

    # ------------------------------------------------------------------ counts
    cnt = []
    for cls in PRIMARY + [
        "pa14_owned_elsewhere",
        "er_housekeeping",
        "enzyme_lectin_like_domain",
        "n1_shared",
        "assumed",
    ]:
        rows = by_class[cls]
        cnt.append(
            {
                "group": cls,
                "label": rows[0]["label"] if rows else "",
                "rows": len(rows),
                "sequences": len({r["accession"] for r in rows}),
                "clusters": len({r["cluster30"] for r in rows}),
                "species": len({r["taxon_id"] for r in rows}),
                "ascomycota_clusters": len(
                    {r["cluster30"] for r in rows if r["phylum"] == "Ascomycota"}
                ),
                "basidiomycota_clusters": len(
                    {r["cluster30"] for r in rows if r["phylum"] == "Basidiomycota"}
                ),
                "leakage_in_reference": sum(r["leakage"] == "in_reference" for r in rows),
                "leakage_none": sum(r["leakage"] == "none" for r in rows),
            }
        )
    prim = [r for c in PRIMARY for r in by_class[c]]
    cnt.append(
        {
            "group": "ALL_PRIMARY_POSITIVES",
            "label": "positive",
            "rows": len(prim),
            "sequences": len(prim),
            "clusters": len({r["cluster30"] for r in prim}),
            "species": len({r["taxon_id"] for r in prim}),
            "ascomycota_clusters": len(
                {r["cluster30"] for r in prim if r["phylum"] == "Ascomycota"}
            ),
            "basidiomycota_clusters": len(
                {r["cluster30"] for r in prim if r["phylum"] == "Basidiomycota"}
            ),
            "leakage_in_reference": sum(r["leakage"] == "in_reference" for r in prim),
            "leakage_none": sum(r["leakage"] == "none" for r in prim),
        }
    )
    reasoned = [
        r
        for c in ("er_housekeeping", "enzyme_lectin_like_domain", "n1_shared")
        for r in by_class[c]
    ]
    cnt.append(
        {
            "group": "ALL_REASONED_NEGATIVES",
            "label": "negative",
            "rows": len(reasoned),
            "sequences": len(reasoned),
            "clusters": len({r["cluster30"] for r in reasoned}),
            "species": len({r["taxon_id"] for r in reasoned}),
            "ascomycota_clusters": len(
                {r["cluster30"] for r in reasoned if r["phylum"] == "Ascomycota"}
            ),
            "basidiomycota_clusters": len(
                {r["cluster30"] for r in reasoned if r["phylum"] == "Basidiomycota"}
            ),
            "leakage_in_reference": sum(r["leakage"] == "in_reference" for r in reasoned),
            "leakage_none": sum(r["leakage"] == "none" for r in reasoned),
        }
    )
    ccols = list(cnt[0].keys())
    write_tsv(a.out_dir + "/counts.tsv", cnt, ccols)

    # mixed clusters (positive and negative members in one cluster)
    clu_labels = defaultdict(set)
    for r in truth:
        clu_labels[r["cluster30"]].add(r["label"])
    mixed = [c for c, s in clu_labels.items() if len(s) > 1]
    with open(a.out_dir + "/mixed_clusters.tsv", "w") as fh:
        fh.write("cluster\tmembers\n")
        for c in sorted(mixed):
            fh.write(
                c
                + "\t"
                + ",".join(
                    r["accession"] + ":" + r["label_class"] for r in truth if r["cluster30"] == c
                )
                + "\n"
            )

    # ------------------------------------------------------------------ per-row hits
    best = defaultdict(dict)
    for acc, m, s in hit_rows:
        best[acc][m] = max(s, best[acc].get(m, -1))
    pr = []
    union15 = lm.hit_accessions(hit_rows, M12)
    for r in truth:
        h = best.get(r["accession"], {})
        pr.append(
            {
                "accession": r["accession"],
                "label": r["label"],
                "label_class": r["label_class"],
                "cluster30": r["cluster30"],
                "species": r["species"],
                "leakage": r["leakage"],
                "seed_hits": r["seed_hits"],
                "expected_model": expected_model(r) if r["label"] == "positive" else "",
                "m12_models_hit": ";".join(f"{m}:{h[m]:.1f}" for m in M12 if m in h),
                "m3_models_hit": ";".join(f"{m}:{h[m]:.1f}" for m in M3 if m in h),
                "union15_hit": "yes" if r["accession"] in union15 else "no",
                "pf07938_score": f"{h['PF07938']:.1f}" if "PF07938" in h else "",
            }
        )
    write_tsv(a.out_dir + "/per_row_hits.tsv", pr, list(pr[0].keys()))

    # ------------------------------------------------------------------ recall by group and model
    def row_for(group_name, rows, set_name, hitset):
        s = lm.summarize(members(rows), hitset)
        return {
            "group": group_name,
            "model_set": set_name,
            "n_seq": s["n_seq"],
            "k_seq": s["k_seq"],
            "seq_ci": lm.fmt(s["k_seq"], s["n_seq"], s["seq_lo"], s["seq_hi"]),
            "n_clu": s["n_clu"],
            "k_clu": s["k_clu"],
            "clu_ci": lm.fmt(s["k_clu"], s["n_clu"], s["clu_lo"], s["clu_hi"]),
            "recall_clu": f"{s['k_clu'] / s['n_clu']:.3f}" if s["n_clu"] else "",
        }

    sets = {
        "union15": union15,
        "rev1_union5": lm.hit_accessions(hit_rows, REV1),
        "m3_union5": lm.hit_accessions(hit_rows, M3),
    }
    for m in M12 + M3:
        sets[m] = lm.hit_accessions(hit_rows, [m])
    pos_groups = {c: by_class[c] for c in PRIMARY + ["pa14_owned_elsewhere"]}
    pos_groups["ALL_PRIMARY_POSITIVES"] = prim
    pos_groups["ALL_PRIMARY_leakage_none"] = [r for r in prim if r["leakage"] == "none"]
    pos_groups["ALL_PRIMARY_in_reference"] = [r for r in prim if r["leakage"] == "in_reference"]
    pos_groups["ALL_PRIMARY_independent_route"] = [
        r for r in prim if r["selection_route"] != "model_gated_pfam_xref"
    ]
    pos_groups["ALL_PRIMARY_model_gated_route"] = [
        r for r in prim if r["selection_route"] == "model_gated_pfam_xref"
    ]
    pos_groups["ALL_PRIMARY_T1"] = [r for r in prim if r["tier"] == "T1"]
    pos_groups["ALL_PRIMARY_ascomycota"] = [r for r in prim if r["phylum"] == "Ascomycota"]
    pos_groups["ALL_PRIMARY_basidiomycota"] = [r for r in prim if r["phylum"] == "Basidiomycota"]
    rec = []
    for g, rows in pos_groups.items():
        for sname, hs in sets.items():
            rec.append(row_for(g, rows, sname, hs))
    write_tsv(a.out_dir + "/recall_by_group.tsv", rec, list(rec[0].keys()))

    # ------------------------------------------------------------------ false hits
    neg_groups = {
        c: by_class[c] for c in ("er_housekeeping", "enzyme_lectin_like_domain", "n1_shared")
    }
    neg_groups["ALL_REASONED_NEGATIVES"] = reasoned
    neg_groups["er_housekeeping_experimental_only"] = [
        r for r in by_class["er_housekeeping"] if r["negative_evidence"] == "experimental_eco269"
    ]
    fh_rows = []
    for g, rows in neg_groups.items():
        for sname, hs in sets.items():
            x = row_for(g, rows, sname, hs)
            x["fhr_clu"] = x.pop("recall_clu")
            fh_rows.append(x)
    write_tsv(a.out_dir + "/false_hits.tsv", fh_rows, list(fh_rows[0].keys()))
    # construction check (review 1, M3): the assumed set was filtered by UniProt Pfam cross-references
    cc = []
    for sname in ("union15", "PF07938"):
        x = row_for("assumed_construction_check", by_class["assumed"], sname, sets[sname])
        x["fhr_clu"] = x.pop("recall_clu")
        cc.append(x)
    write_tsv(a.out_dir + "/construction_check.tsv", cc, list(cc[0].keys()))

    # ------------------------------------------------------------------ PF07938 floors
    fl = []
    for floor in FLOORS:
        fm = {"PF07938": floor} if floor else {}
        pf = lm.hit_accessions(hit_rows, ["PF07938"], fm)
        others = [m for m in M12 if m != "PF07938"]
        un = lm.hit_accessions(hit_rows, ["PF07938"], fm) | lm.hit_accessions(hit_rows, others)
        label = "GA (25.3)" if floor is None else f"score>={int(floor)}"
        fam = [r for r in prim if expected_model(r) == "PF07938"]
        for g, rows in (
            ("PF07938_family_positives", fam),
            ("ALL_PRIMARY_POSITIVES", prim),
            ("ALL_REASONED_NEGATIVES", reasoned),
        ):
            for sname, hs in (("PF07938_only", pf), ("union15_with_floor", un)):
                x = row_for(g, rows, sname, hs)
                x["floor"] = label
                fl.append(x)
    write_tsv(a.out_dir + "/floors.tsv", fl, ["floor"] + [k for k in fl[0] if k != "floor"])

    # ------------------------------------------------------------------ family recall + seed split (cluster level)
    fam_rows = []
    SPLIT_SIM = "cluster has a member with 30% identity or more to a seed of this model"
    SPLIT_IND = "cluster has no such member (seed-independent)"

    def fam_row(model, split, rows, hs, variant="all"):
        x = row_for("family", rows, model, hs)
        x.pop("group")
        x.pop("model_set")
        return {"model": model, "variant": variant, "split": split, **x}

    for m in M12:
        fam = [r for r in prim if expected_model(r) == m]
        if not fam:
            continue
        hs = sets[m]
        variants = [("all", fam)]
        if m == "PF14200":
            variants.append(
                (
                    "excluding_laetiporus_trefoil_rows",
                    [r for r in fam if r["unilectin_class"] != "Laetiporus trefoil lectin-like"],
                )
            )
        for variant, rows_v in variants:
            fam_rows.append(fam_row(m, "all", rows_v, hs, variant))
            flagged = {r["accession"] for r in rows_v if m + ":" in r["seed_hits"]}
            sim, ind = lm.cluster_seed_split(members(rows_v), flagged)
            sim_acc = {a_ for a_, _ in sim}
            ind_acc = {a_ for a_, _ in ind}
            fam_rows.append(
                fam_row(m, SPLIT_SIM, [r for r in rows_v if r["accession"] in sim_acc], hs, variant)
            )
            fam_rows.append(
                fam_row(m, SPLIT_IND, [r for r in rows_v if r["accession"] in ind_acc], hs, variant)
            )
    write_tsv(a.out_dir + "/family_recall.tsv", fam_rows, list(fam_rows[0].keys()))
    # families with no model
    nomodel = [r for r in prim if expected_model(r) in ("none", "unmapped")]
    s = row_for("no_expected_model", nomodel, "union15", union15)
    with open(a.out_dir + "/no_model_families.tsv", "w") as fh:
        fh.write("n_seq\tk_seq_hit_by_any_of_15\tn_clu\tk_clu\n")
        fh.write(f"{s['n_seq']}\t{s['k_seq']}\t{s['n_clu']}\t{s['k_clu']}\n")

    # ------------------------------------------------------------------ leakage split for the union
    # Row level: a row is in_reference or none. Clusters can span both, so cluster counts of the two
    # groups can add up to more than the clusters of the group. Cluster level: a cluster goes to one group.
    ls = []
    for g in (
        "structure_backed",
        "host_binding_h1a",
        "chitin_masking_h1b",
        "literature_binding",
        "ALL_PRIMARY_POSITIVES",
    ):
        rows = pos_groups[g]
        for lk in ("none", "in_reference"):
            sub = [r for r in rows if r["leakage"] == lk]
            x = row_for(g, sub, "union15", union15)
            x["leakage"] = lk
            x["basis"] = "row_level (clusters can be counted in both groups)"
            ls.append(x)
        flagged = {r["accession"] for r in rows if r["leakage"] == "in_reference"}
        sim, ind = lm.cluster_seed_split(members(rows), flagged)
        for lk, grp in (("cluster_has_in_reference_member", sim), ("cluster_all_none", ind)):
            acc_set = {a_ for a_, _ in grp}
            x = row_for(g, [r for r in rows if r["accession"] in acc_set], "union15", union15)
            x["leakage"] = lk
            x["basis"] = "cluster_level (each cluster in one group)"
            ls.append(x)
    write_tsv(
        a.out_dir + "/leakage_split.tsv",
        ls,
        ["basis", "leakage"] + [k for k in ls[0] if k not in ("basis", "leakage")],
    )

    # identity to the nearest seed sequence (any of the 20 models), primary positives
    def best_ident(r):
        v = [float(x.split(":")[1]) for x in r["seed_hits"].split(";") if x]
        return max(v) if v else None

    bins = defaultdict(list)
    for r in prim:
        bins[lm.identity_bin(best_ident(r))].append(r)
    order = ["98% or more", "90 to 97.9%", "50 to 89.9%", "30 to 49.9%", "no seed hit (<30%)"]
    ib = []
    for b in order:
        x = row_for(b, bins.get(b, []), "union15", union15)
        x["identity_bin"] = x.pop("group")
        x.pop("model_set")
        ib.append(x)
    write_tsv(
        a.out_dir + "/identity_bins.tsv",
        ib,
        ["identity_bin"] + [k for k in ib[0] if k != "identity_bin"],
    )
    # strata that share clusters
    clu_strata = defaultdict(set)
    for r in prim:
        clu_strata[r["cluster30"]].add(r["label_class"])
    with open(a.out_dir + "/strata_shared_clusters.tsv", "w") as fh:
        fh.write("cluster\tstrata\n")
        for c, st in sorted(clu_strata.items()):
            if len(st) > 1:
                fh.write(c + "\t" + ",".join(sorted(st)) + "\n")
        fh.write(
            f"#sum_of_stratum_cluster_counts\t{sum(len(v) for v in clu_strata.values())}\t"
            f"distinct_clusters\t{len(clu_strata)}\n"
        )

    # ------------------------------------------------------------------ misses
    sub = {}
    pf = {}
    if a.miss_dir:
        p1 = a.miss_dir + "/miss_subga.tsv"
        p2 = a.miss_dir + "/miss_pfamscan.tsv"
        if os.path.exists(p1):
            for r in read_tsv(p1):
                sub.setdefault(r["accession"], []).append(
                    f"{r['model']}:{r['score']}(GA {r['ga']})"
                )
        if os.path.exists(p2):
            for r in read_tsv(p2):
                pf.setdefault(r["accession"], []).append(f"{r['family']}:{r['score']}")
    ms = []
    for r in truth:
        if (
            r["label"] != "positive"
            or r["label_class"] == "pa14_owned_elsewhere"
            or r["accession"] in union15
        ):
            continue
        em = expected_model(r)
        ms.append(
            {
                "accession": r["accession"],
                "species": r["species"],
                "label_class": r["label_class"],
                "unilectin_class": r["unilectin_class"],
                "protein_name": r["protein_name"][:50],
                "length": r["length"],
                "expected_model": em,
                "reason": (
                    "no model among the 15"
                    if em in ("none", "unmapped")
                    else (
                        "below GA of the expected model"
                        if any(x.startswith(em + ":") for x in sub.get(r["accession"], []))
                        else "no hit to the expected model (E<=10)"
                    )
                ),
                "sub_ga_hits_to_the_20_models": ";".join(sub.get(r["accession"], []))
                or "none (E<=10)",
                "full_pfam_ga_hits": ";".join(pf.get(r["accession"], [])) or "none",
            }
        )
    write_tsv(a.out_dir + "/misses.tsv", ms, list(ms[0].keys()) if ms else ["accession"])

    # ------------------------------------------------------------------ species clusters (20-cluster floor)
    sp = defaultdict(lambda: defaultdict(set))
    for r in truth:
        if r["label"] == "positive":
            key = "pos_all" if r["label_class"] == "pa14_owned_elsewhere" else "pos_primary"
        else:
            key = "neg_assumed" if r["label_class"] == "assumed" else "neg_reasoned"
            if r["label_class"] == "n1_shared":
                sp[r["species"]]["neg_n1"].add(r["cluster30"])
        sp[r["species"]][key].add(r["cluster30"])
    rows = [
        {
            "species": k,
            "taxon_id": next(x["taxon_id"] for x in truth if x["species"] == k),
            "primary_positive_clusters": len(v["pos_primary"]),
            "pa14_positive_clusters": len(v["pos_all"]),
            "reasoned_negative_clusters": len(v["neg_reasoned"]),
            "n1_negative_clusters": len(v["neg_n1"]),
            "assumed_negative_clusters": len(v["neg_assumed"]),
        }
        for k, v in sp.items()
    ]
    rows.sort(
        key=lambda x: (
            -x["primary_positive_clusters"],
            -x["n1_negative_clusters"],
            x["species"],
        )
    )
    write_tsv(a.out_dir + "/species_clusters.tsv", rows, list(rows[0].keys()))
    print("ok; union15 hits", len(union15), "of", len(truth))


main()
