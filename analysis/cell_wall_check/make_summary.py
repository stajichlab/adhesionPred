#!/usr/bin/env python3.12
"""Join scan.hits.tsv, full_pfam_ga_hits.tsv, uniprot_features.tsv, the committed shipped-model scores
(analysis/model_review/s288c_scores.tsv) and the committed S288C sorting-hat rows
(docs/reports/data/sorting_hat/classifier_check/curated_label_check.rows.tsv) into summary.tsv.
Nothing is re-scored here. Run from this directory."""
import csv, collections
R = "../../"
T = list(csv.DictReader(open("test_set.tsv"), delimiter="\t"))
fam = {r["pfam_acc"]: r for r in csv.DictReader(open(R + "data/sorting_hat/family_table.tsv"), delimiter="\t")}
hits = collections.defaultdict(lambda: collections.defaultdict(list))
for r in csv.DictReader(open("scan.hits.tsv"), delimiter="\t"):
    g = r["protein"].split("|")[1]
    hits[g][(r["model"], r["model_name"], r["seq_score"])].append(r["dom_score"])
feat = {r["gene"]: r for r in csv.DictReader(open("uniprot_features.tsv"), delimiter="\t")}
sc = {r["id"]: r["p_shipped_t6_8M"] for r in csv.DictReader(open(R + "analysis/model_review/s288c_scores.tsv"), delimiter="\t")}
rows = {}
for r in csv.DictReader(open(R + "docs/reports/data/sorting_hat/classifier_check/curated_label_check.rows.tsv"), delimiter="\t"):
    if r["accession"] == "E9P9G2": continue
    rows[r["accession"]] = r
other = collections.defaultdict(list)
for r in csv.DictReader(open("full_pfam_ga_hits.tsv"), delimiter="\t"):
    if r["pfam"] not in fam and r["name"] not in ("Glyco_hydro_72","Glyco_hydro_16","Glyco_hydro_17","Glyco_hydro_76","Cellulase","Glyco_hydro_18","SUN","Kre9_KNH","X8"):
        other[r["gene"]].append(f"{r['name']}({r['pfam']}) {r['best_domain_score']}")
with open("summary.tsv", "w") as o:
    o.write("gene\tuniprot\tsystematic_id\tgroup\tfamily_table_hits\tcontext_enzyme_hits\tother_pfam_hits_full_scan\tuniprot_gpi_site\tshipped_esm_p\tR0_signal_peptide_committed_run\ttandem_repeat_committed_run\n")
    for t in T:
        ft, cx = [], []
        for (m, n, s), d in hits[t["gene"]].items():
            txt = f"{n}({m}) seq {s}, {len(d)} dom" if len(d) > 1 else f"{n}({m}) {s}"
            if m in fam: ft.append(txt + f" [{'active' if fam[m]['active']=='yes' else 'inactive'}" + (f", cond {fam[m]['second_condition']}]" if fam[m]['second_condition'] else "]"))
            else: cx.append(txt)
        c = rows.get(t["uniprot"], {})
        o.write("\t".join([t["gene"], t["uniprot"], t["systematic_id"], t["group"], "; ".join(ft) or "none", "; ".join(cx) or "none",
                           "; ".join(other[t["gene"]]) or "none", feat[t["gene"]]["gpi_site"], sc.get(t["systematic_id"], "NA"),
                           c.get("signal_peptide_protein[R0]", "not in committed rows"), c.get("tandem_repeat_protein", "not in committed rows")]) + "\n")
