#!/usr/bin/env python3
"""Rank Coccidioides proteins as serodiagnostic antigen candidates.

Rewritten after an independent review found four defects in the first version:
  1. SOWgp and Ag2/PRA were absent from the ranking entirely -- the C. posadasii reference
     proteome used was strain "C735 delta SOWgp", a SOWgp DELETION strain.
  2. The complement-fixation antigen actually used in clinical serology (endochitinase-1,
     CiX1) scored 4/8, tied with 678 other proteins.
  3. The score penalized human homology but NOT cross-reactivity with the fungi that actually
     confound coccidioidomycosis serology: Histoplasma, Blastomyces, Paracoccidioides,
     Aspergillus. The old top hits (Gel1, pepsins, subtilisins, chitinases) are pan-fungal
     conserved, i.e. maximally cross-reactive.
  4. 63% of rows tied at a single integer score.

This version scores the C. immitis RS reference proteome, uses a continuous score, adds the
fungal cross-reactivity penalty, and adds pangenome prevalence across 489 Coccidioides
proteomes. It is held to an ACCEPTANCE TEST: the known antigens must rank near the top, or
the ranking is not trustworthy on unknown proteins.
"""

import argparse
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis._common.paths import COCCI_PANGENOME  # noqa: E402

PAN = COCCI_PANGENOME
OG = PAN / "results/Cocci_496_OG2_5_5/Orthogroups"
REF_COL = "CimmitisRS_FungiDB"
CONFOUNDERS = ["Histoplasma", "Blastomyces", "Paracoccidioides", "Aspergillus_fumigatus"]


def read_m8(path, cols=("query", "target", "pident", "qcov")):
    out = defaultdict(list)
    if not Path(path).exists():
        return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        out[f[0]].append(dict(zip(cols, f)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--work", required=True, help="dir holding idmap.m8, crossreact.m8, human.m8, anchors.m8"
    )
    ap.add_argument("--antigens", required=True, help="curated IEDB antigen table (antigens.tsv)")
    ap.add_argument("--iedb-hits", default=None, help="optional m8 of reference vs IEDB antigens")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    W = Path(args.work)

    # --- orthogroup membership and per-strain gene counts
    prot2og, og_members = {}, {}
    with open(OG / "Orthogroups.tsv") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ref_i = hdr.index(REF_COL)
        for row in rd:
            og_members[row[0]] = row
            for p in (row[ref_i] or "").split(", "):
                if p.strip():
                    prot2og[p.strip()] = row[0]
    print(f"{len(prot2og)} reference proteins assigned to orthogroups")

    counts = {}
    with open(OG / "Orthogroups.GeneCount.tsv") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        strain_idx = [i for i, h in enumerate(hdr) if i and h != "Total"]
        for row in rd:
            v = [int(row[i]) for i in strain_idx]
            counts[row[0]] = v
    n_strains = len(strain_idx)
    print(f"{len(counts)} orthogroups x {n_strains} proteomes")

    # --- evidence tables
    idmap = read_m8(W / "idmap.m8")
    human = read_m8(W / "human.m8")
    anchors = read_m8(W / "anchors.m8")
    cross = defaultdict(dict)
    if (W / "crossreact.m8").exists():
        for line in open(W / "crossreact.m8"):
            f = line.rstrip("\n").split("\t")
            genus, q, pid = f[0], f[1], float(f[3])
            cross[q][genus] = max(cross[q].get(genus, 0), pid)

    anchor_of = {}
    for q, hits in anchors.items():
        best = max(hits, key=lambda h: float(h["pident"]))
        anchor_of.setdefault(best["target"], []).append(q)

    iedb = read_m8(args.iedb_hits) if args.iedb_hits else {}

    rows = []
    for prot, og in prot2og.items():
        c = counts.get(og, [])
        present = sum(1 for x in c if x > 0)
        prevalence = present / n_strains if n_strains else 0
        cn = [x for x in c if x > 0]
        cn_mean = st.mean(cn) if cn else 0
        cn_cv = (st.pstdev(cn) / cn_mean) if cn and cn_mean else 0

        xr = cross.get(prot, {})
        max_xr = max(xr.values()) / 100 if xr else 0.0  # 0..1, higher = worse
        n_xr = len(xr)
        hs = max((float(h["pident"]) for h in human.get(prot, [])), default=0) / 100
        ag = max((float(h["pident"]) for h in iedb.get(prot, [])), default=0) / 100

        # --- TWO reported axes, kept separate, because they answer different questions and
        # the first version wrongly merged them.
        #
        # antigenicity: is this protein likely to be seen by the immune system at all?
        # specificity : would an assay built on it distinguish Coccidioides from the fungi that
        #               actually confound coccidioidomycosis serology?
        #
        # The complement-fixation antigen in clinical use (endochitinase CiX1) is the reason
        # this split matters: it is genuinely immunogenic AND genuinely cross-reactive
        # (64% identical to confounder chitinases in our own data). Scoring it on one merged
        # axis either rewards a cross-reactive antigen or punishes a real one. It is therefore
        # a positive control for antigenicity and a NEGATIVE control for specificity.
        antigenicity = 2.5 * prevalence + 2.0 * ag
        specificity = -3.0 * max_xr - 0.5 * (n_xr / len(CONFOUNDERS)) - 1.5 * hs
        # Copy-number variability is NOT part of the antigen score. It is interesting biology
        # and it answers the separate pangenome-variability question, but for a serodiagnostic
        # a variable multi-copy family is if anything a liability. Reported, not scored.
        score = antigenicity + specificity
        rows.append(
            {
                "protein": prot,
                "orthogroup": og,
                "anchor": ";".join(anchor_of.get(prot, [])),
                "score": round(score, 4),
                "antigenicity": round(antigenicity, 4),
                "specificity": round(specificity, 4),
                "prevalence": round(prevalence, 4),
                "n_strains_present": present,
                "copy_number_mean": round(cn_mean, 2),
                "copy_number_cv": round(cn_cv, 3),
                "max_fungal_crossreact_pid": round(max_xr * 100, 1),
                "n_confounder_genera": n_xr,
                "human_homolog_pid": round(hs * 100, 1),
                "iedb_antigen_pid": round(ag * 100, 1),
                "fungi5k_id": (idmap.get(prot) or [{}])[0].get("target", ""),
            }
        )

    # One entry per orthogroup. Without this a single expanded family floods the top: the
    # 14 reference members of OG0000001 all carry an identical pangenome profile and occupied
    # 14 of the first 50 rows.
    best = {}
    for r in rows:
        k = r["orthogroup"]
        if (
            k not in best
            or r["score"] > best[k]["score"]
            or (r["anchor"] and not best[k]["anchor"])
        ):
            best[k] = r
    for r in rows:
        r["is_orthogroup_representative"] = "yes" if best.get(r["orthogroup"]) is r else "no"
    rows.sort(key=lambda r: -r["score"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
        r["percentile"] = round(100 * i / len(rows), 2)
    reps = [r for r in rows if r["is_orthogroup_representative"] == "yes"]
    for i, r in enumerate(reps, 1):
        r["rank_dedup"] = i
        r["percentile_dedup"] = round(100 * i / len(reps), 2)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(
        f"wrote {len(rows)} ranked proteins to {args.out} "
        f"({sum(1 for r in rows if r['is_orthogroup_representative'] == 'yes')} orthogroup representatives)"
    )

    for i, r in enumerate(reps, 1):
        r["rank_dedup"] = i
        r["percentile_dedup"] = round(100 * i / len(reps), 2)
    for key in ("antigenicity", "specificity"):
        ordered = sorted(rows, key=lambda r: -r[key])
        for i, r in enumerate(ordered, 1):
            r[f"{key}_percentile"] = round(100 * i / len(rows), 2)

    print("\n=== ACCEPTANCE TEST ===")
    print("Coccidioides-SPECIFIC antigens (SOWgp, PRA family) must reach the top decile on the")
    print("combined score. The CF antigen is a positive control for antigenicity and a NEGATIVE")
    print("control for specificity -- it is the known cross-reactive one.\n")
    print(
        f"{'anchor':<24}{'protein':<20}{'comb%':>7}{'antig%':>8}{'spec%':>7}{'xreact':>8}  verdict"
    )
    anchored = [r for r in rows if r["anchor"]]
    specific, ok = [], 0
    for r in sorted(anchored, key=lambda r: r["rank"]):
        is_cf = "CFantigen" in r["anchor"]
        if is_cf:
            # expected: high antigenicity, LOW specificity
            good = r["antigenicity_percentile"] <= 25 and r["specificity_percentile"] > 50
            verdict = "PASS (cross-reactive, as expected)" if good else "unexpected"
        else:
            specific.append(r)
            good = r["percentile"] <= 10
            ok += good
            verdict = "PASS" if good else "FAIL"
        print(
            f"{r['anchor'][:23]:<24}{r['protein'][:19]:<20}{r['percentile']:>6.1f}%"
            f"{r['antigenicity_percentile']:>7.1f}%{r['specificity_percentile']:>6.1f}%"
            f"{r['max_fungal_crossreact_pid']:>7.0f}%  {verdict}"
        )
    print(f"\n{ok}/{len(specific)} Coccidioides-specific anchors in the top decile")
    if specific and ok < len(specific):
        print("NOT CALIBRATED on specific antigens -- do not trust the ranking on unknowns.")
    else:
        print(
            f"Calibrated on the available controls. Still only {len(anchored)} controls; "
            "treat accordingly."
        )


if __name__ == "__main__":
    main()
