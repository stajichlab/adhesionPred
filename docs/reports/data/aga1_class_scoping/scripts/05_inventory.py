#!/usr/bin/env python3
"""Inventory counts (proteins and MMseqs2 clusters) per role, source table and clade. Reads proteins.tsv, features.tsv,
clusters.tsv. Writes source_counts.tsv, sequence_availability.tsv, class_clade_counts.tsv, lectin_overlap.tsv (counts only), kiss_adhesion_rows.tsv.

Negative subsets (pre-declared before any model was run):
  hn_all        every hard_negative protein
  hn_st30       hard_negative with max_ST30 >= 0.5 (some 30-residue window at least half Ser/Thr)
  hn_stfrac25   hard_negative with whole-protein S+T fraction >= 0.25
  hn_repeat     hard_negative with >= 2 UniProt repeat features
  bg            background_assumed_negative (S288C / C. albicans proteome members with no curated class; negatives ASSUMED)
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import OUT, ROOT, read_tsv, write_tsv  # noqa: E402


def load():
    P = {r["accession"]: r for r in read_tsv(OUT / "proteins.tsv")}
    F = {r["accession"]: r for r in read_tsv(OUT / "features.tsv")}
    C = {r["accession"]: r for r in read_tsv(OUT / "clusters.tsv")}
    for a, p in P.items():
        p.update(F[a])
        p["cluster"] = C[a]["mmseqs_cluster"]
    return P


def subsets(p):
    """Names of the evaluation groups a protein belongs to."""
    out = []
    if p["role"] == "class_positive":
        out.append("pos_all")
        if p["pos_tier"] == "A":
            out.append("pos_tierA")
        if p["uniprot_gpi_feature"] == "yes":
            out.append("pos_gpi_feature")
    elif p["role"] == "hard_negative":
        out.append("hn_all")
        if float(p["max_ST30"]) >= 0.5:
            out.append("hn_st30")
        if float(p["frac_ST"]) >= 0.25:
            out.append("hn_stfrac25")
        if p["uniprot_repeat_n"] not in ("", "0", "1"):
            out.append("hn_repeat")
    elif p["role"] == "background_assumed_negative":
        out.append("bg")
    return out


def main():
    P = load()
    # 1. source x role
    c = Counter()
    for p in P.values():
        for s in p["sources"].split(";"):
            c[(s, p["role"])] += 1
    write_tsv(
        OUT / "source_counts.tsv",
        [{"source": s, "role": r, "n_proteins": n} for (s, r), n in sorted(c.items())],
        ["source", "role", "n_proteins"],
    )
    sa = Counter((p["role"], p["sequence_source"]) for p in P.values())
    write_tsv(OUT / "sequence_availability.tsv", [{"role": r, "sequence_source": s_, "n_proteins": n} for (r, s_), n in sorted(sa.items())], ["role", "sequence_source", "n_proteins"])  # fmt: skip
    # 2. subset x clade
    prot, clus = defaultdict(int), defaultdict(set)
    for p in P.values():
        for s in subsets(p):
            for clade in (p["clade_fine"], "ALL"):
                prot[(s, clade)] += 1
                clus[(s, clade)].add(p["cluster"])
    rows = [
        {"subset": s, "clade_fine": cl, "n_proteins": n, "n_clusters": len(clus[(s, cl)])}
        for (s, cl), n in sorted(prot.items())
    ]
    write_tsv(
        OUT / "class_clade_counts.tsv", rows, ["subset", "clade_fine", "n_proteins", "n_clusters"]
    )
    # 3. lectin controls overlap (counts only)
    lect = read_tsv(ROOT / "analysis/lectin_validation/controls/controls.tsv")
    rows = []
    cnt = Counter()
    for r in lect:
        cnt[
            (
                r["label"],
                r["stratum"].split(":")[0]
                if r["stratum"].startswith("n1_shared")
                else r["stratum"],
                P[r["accession"]]["role"] if r["accession"] in P else "not_in_scoping_table",
            )
        ] += 1
    for (lab, strat, role), n in sorted(cnt.items()):
        rows.append(
            {"lectin_label": lab, "lectin_stratum": strat, "role_in_scoping_table": role, "n": n}
        )
    write_tsv(
        OUT / "lectin_overlap.tsv",
        rows,
        ["lectin_label", "lectin_stratum", "role_in_scoping_table", "n"],
    )
    # 4. Kiss 2019 adhesion-slice rows (wall-fusion-kiss2019/controls.tsv; the JGI-only rows have no committed sequence)
    krows = []
    for r in read_tsv(ROOT / "data/controls/wall-fusion-kiss2019/controls.tsv"):
        if r["stratum"] != "adhesion":
            continue
        p = P[r["accession"]]
        krows.append({"accession": r["accession"], "gene": r["gene"], "clade_fine": r["clade_fine"], "role_in_scoping_table": p["role"], "length": p["length"], "frac_ST": p["frac_ST"], "max_ST30": p["max_ST30"], "call_repeat": p["call_repeat"], "kiss_cluster": r["kiss_cluster_number"]})  # fmt: skip
    write_tsv(OUT / "kiss_adhesion_rows.tsv", krows, list(krows[0]))
    print(len(P), "proteins")


if __name__ == "__main__":
    main()
