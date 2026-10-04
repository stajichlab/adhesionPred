#!/usr/bin/env python3
"""Summarise the follow-up checks (11_followup_job.sh) for the 45 follow-up genes.

Per gene: ortholog hits by group (phmmer), gene-model agreement with the pangenome RS annotation
(diamond), Pfam domains (hmmscan), Phobius calls. Also counts over all spherule-up genes for Pfam
and Phobius. Usage: 12_summarize_followup.py --followup <dir> --table <table.tsv.gz> --out <tsv>
Standard library only.
"""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

DIMORPHIC = {"Histoplasma", "Blastomyces", "Paracoccidioides", "Emergomyces", "Emmonsia",
             "Emmonsiellopsis", "Ajellomyces", "Lomentospora"}  # fmt: skip
STRICT = 1e-5
GROUPS = ("coccidioides", "dimorphic", "other_onygenales", "outgroup")


def group_of(row):
    g = row["genus"]
    if g == "Coccidioides":
        return "coccidioides"
    if g in DIMORPHIC:
        return "dimorphic"
    return "other_onygenales" if row["order"] == "Onygenales" else "outgroup"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--followup", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fu = Path(a.followup)
    genes = list(csv.DictReader(open(fu / "followup_genes.tsv"), delimiter="\t"))
    table = {r["gene_id"]: r for r in csv.DictReader(gzip.open(a.table, "rt"), delimiter="\t")}

    stem_group, stem_species = {}, {}
    for r in csv.DictReader(open(fu / "search_proteomes.tsv"), delimiter="\t"):
        stem = Path(r["file"]).name[: -len(".proteins.fa")]
        stem_group[stem] = group_of(r)
        stem_species[stem] = r["species"]

    # phmmer: best E-value per (query gene, group), species with a hit per group
    best = defaultdict(lambda: 1e9)
    species_hit = defaultdict(set)
    species_strict = defaultdict(set)
    in_fungi5k_rs = defaultdict(bool)
    for line in open(fu / "results" / "phmmer.tblout"):
        if line.startswith("#"):
            continue
        f = line.split()
        stem = f[0].split("|")[0]
        gene = f[2].split("|")[0]
        ev = float(f[4])
        if stem == "Coccidioides_immitis_RS":  # Fungi5k annotation of the same genome: model check
            in_fungi5k_rs[gene] = True
            continue
        grp = stem_group.get(stem)
        if grp is None:
            continue
        best[(gene, grp)] = min(best[(gene, grp)], ev)
        species_hit[(gene, grp)].add(stem)
        if ev <= STRICT:
            species_strict[(gene, grp)].add(stem)

    # diamond against the pangenome RS proteome: first (best) hit per query
    model = {}
    for line in open(fu / "results" / "diamond_rs_pangenome.tsv"):
        f = line.rstrip("\n").split("\t")
        q = f[0].split("|")[0]
        if q in model:
            continue
        pid, length, qlen, slen = float(f[2]), int(f[3]), int(f[4]), int(f[5])
        model[q] = (pid, length / qlen, length / slen)

    # Pfam domains
    pfam = defaultdict(list)
    for line in open(fu / "results" / "pfam.domtblout"):
        if line.startswith("#"):
            continue
        f = line.split()
        gene = f[3].split("|")[0]
        pfam[gene].append((f[0], f[1].split(".")[0], float(f[12])))

    # Phobius
    phob = {}
    for line in open(fu / "results" / "phobius_short.txt"):
        f = line.split()
        if f[0] == "SEQENCE":
            continue
        phob[f[0].split("|")[0]] = (int(f[1]), f[2] == "Y")

    up_genes = {r["gene_id"] for r in table.values() if r["flag_spherule_up"] == "yes"}
    sp_signalp = {g for g in up_genes if table[g]["signalp_call"] == "SP"}
    sp_phobius = {g for g in up_genes if g in phob and phob[g][1]}
    cfem = sorted(g for g in up_genes if any(d[1] == "PF05730" for d in pfam[g]))
    with_domain = sum(1 for g in up_genes if pfam[g])
    with open(Path(a.out).with_name("phobius_spherule_up.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "phobius_sp", "phobius_tm", "signalp_call", "signalp_sp_prob"])
        for g in sorted(up_genes):
            w.writerow([g, "yes" if g in sp_phobius else "no", phob.get(g, ("", False))[0],
                        table[g]["signalp_call"], table[g]["signalp_sp_prob"]])  # fmt: skip
    print(
        f"spherule-up genes: {len(up_genes)}; Phobius SP: {len(sp_phobius)}; SignalP SP: {len(sp_signalp)}; "
        f"both: {len(sp_signalp & sp_phobius)}; Phobius only: {len(sp_phobius - sp_signalp)}"
    )
    print(f"genes with a Pfam domain: {with_domain}; CFEM (PF05730): {cfem}")

    cols = ["gene_id", "set", "length", "log2fc_48h", "cys_frac", "ranking_specific",
            "model_pident", "model_qcov", "model_scov", "model_call", "in_fungi5k_RS_model",
            "coccidioides_best_E", "coccidioides_nsp", "dimorphic_best_E", "dimorphic_nsp",
            "other_ony_best_E", "other_ony_nsp", "outgroup_best_E", "outgroup_nsp",
            "pfam_domains", "phobius_sp", "phobius_tm", "signalp_sp_prob", "tmhmm_helices"]  # fmt: skip
    out_rows = []
    for g in genes:
        gid = g["gene_id"]
        t = table[gid]
        m = model.get(gid)
        if m is None:
            call = "no_hit"
        elif m[0] >= 95 and m[1] >= 0.9 and m[2] >= 0.9:
            call = "agree"
        else:
            call = "partial_or_different"
        row = {
            "gene_id": gid, "set": g["set"], "length": g["length"], "log2fc_48h": g["log2fc_48h"],
            "cys_frac": g["cys_frac"], "ranking_specific": t["flag_cocci_specific"],
            "model_pident": f"{m[0]:.1f}" if m else "", "model_qcov": f"{m[1]:.2f}" if m else "",
            "model_scov": f"{m[2]:.2f}" if m else "", "model_call": call,
            "in_fungi5k_RS_model": "yes" if in_fungi5k_rs[gid] else "no",
            "pfam_domains": ";".join(sorted({f"{d[0]}({d[1]})" for d in pfam[gid]})),
            "phobius_sp": "yes" if phob.get(gid, (0, False))[1] else "no",
            "phobius_tm": phob.get(gid, ("", False))[0], "signalp_sp_prob": t["signalp_sp_prob"],
            "tmhmm_helices": t["tm_helices"],
        }  # fmt: skip
        for grp, key in zip(
            GROUPS, ("coccidioides", "dimorphic", "other_ony", "outgroup"), strict=True
        ):
            b = best.get((gid, grp))
            row[f"{key}_best_E"] = f"{b:.0e}" if b is not None else ""
            row[f"{key}_nsp"] = len(species_strict.get((gid, grp), ()))
        out_rows.append(row)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)

    n = len(out_rows)

    def count(pred):
        return sum(1 for r in out_rows if pred(r))

    n_agree = count(lambda r: r["model_call"] == "agree")
    n_partial = count(lambda r: r["model_call"] == "partial_or_different")
    n_nohit = count(lambda r: r["model_call"] == "no_hit")
    n_absent = count(lambda r: r["in_fungi5k_RS_model"] == "no")
    print(
        f"\n{n} follow-up genes. Gene model vs pangenome RS (diamond): agree {n_agree}, "
        f"partial/different {n_partial}, no hit {n_nohit}; absent from Fungi5k RS model {n_absent}"
    )
    for key, name in (
        ("dimorphic", "dimorphic relatives"),
        ("other_ony", "other Onygenales"),
        ("outgroup", "outgroups"),
    ):
        n_strict = count(lambda r, k=key: int(r[f"{k}_nsp"]) > 0)
        n_none = count(lambda r, k=key: r[f"{k}_best_E"] == "")
        print(f"  strict hit (E<={STRICT:g}) in {name}: {n_strict}; none at E<=1e-3: {n_none}")
    n_cocc = count(lambda r: int(r["coccidioides_nsp"]) > 0)
    n_pho = count(lambda r: r["phobius_sp"] == "yes")
    n_dom = count(lambda r: bool(r["pfam_domains"]))
    n_flag_wrong = count(
        lambda r: (
            r["ranking_specific"] == "yes"
            and int(r["dimorphic_nsp"]) + int(r["other_ony_nsp"]) + int(r["outgroup_nsp"]) > 0
        )
    )
    print(f"  strict hit in another Coccidioides proteome: {n_cocc}")
    print(f"  Phobius SP: {n_pho}; Pfam domain: {n_dom}")
    print(
        f"  flagged specific, with a strict hit in dimorphic, other Onygenales or outgroup: {n_flag_wrong}"
    )


if __name__ == "__main__":
    main()
