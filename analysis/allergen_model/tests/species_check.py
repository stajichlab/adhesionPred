#!/usr/bin/env python3.12
"""Species-level checks on the candidate and positive sets (stdlib only).

Inputs (from the full runs):
  --candidates  analysis/allergen_model/candidates/candidates.tsv   (accession, organism_id, organism_name, genus, ...)
  --positives   analysis/allergen_model/positives/positives.tsv     (accession, in_iuis, in_kw0020, evidence_level, organism, ...)
  --pairs       full pairs table (zstd): query, target, fident, alnlen, qcov, tcov, evalue, bits   (query N_, target P_)
  --outdir      output folder

Checks (each prints PASS or FAIL; the script exits 1 if a structural check fails):
  S1  No accession is both a positive and a candidate (the KW-0020 and IUIS exclusion works).
  S2  Every candidate has a genus.
  S3  Pairing counts at T_C 0.6 and 0.7 with the O7 coverage rule (qcov and tcov >= 0.5), per genus.
      A candidate that pairs with a known allergen is a possible unannotated allergen (label noise, spec 3.2).
  S4  Per target genus: positives (IUIS and KW-0020), candidates, and candidates that pair (T_C 0.6).

Writes species_check.tsv and species_check.md to --outdir.
"""

import argparse
import collections
import csv
import os
import subprocess
import sys


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_pairs(path):
    out = subprocess.run(["zstd", "-dc", path], capture_output=True, text=True, check=True).stdout
    lines = out.strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=False)) for ln in lines[1:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--positives", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cands = read_tsv(args.candidates)
    poss = read_tsv(args.positives)
    pairs = read_pairs(args.pairs)
    fails = []

    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + (f" [{detail}]" if detail else ""))
        if not ok:
            fails.append(name)

    cand_acc = {r["accession"] for r in cands}
    pos_acc = {r["accession"] for r in poss}
    overlap = cand_acc & pos_acc
    check("S1 no accession is both candidate and positive", not overlap, str(sorted(overlap)[:5]))
    check("S2 every candidate has a genus", all(r["genus"] for r in cands))

    pos_genus = {
        r["accession"]: (r["organism"].split(" ")[0] if r["organism"] else "unknown") for r in poss
    }
    cand_genus = {r["accession"]: r["genus"] for r in cands}

    # S3: best pairing per candidate at each T_C, under the O7 coverage rule.
    best = {}
    for p in pairs:
        q, t = p["query"][2:], p["target"][2:]  # strip N_ / P_
        fid, qc, tc = float(p["fident"]), float(p["qcov"]), float(p["tcov"])
        if qc < 0.5 or tc < 0.5:
            continue
        if q not in best or fid > best[q][0]:
            best[q] = (fid, t)
    check("S3a every pair query is a candidate", all(q in cand_acc for q in best), "")
    check("S3b every pair target is a positive", all(t in pos_acc for _, t in best.values()), "")

    rows = []
    summary = collections.defaultdict(lambda: collections.Counter())
    for acc in sorted(cand_acc):
        g = cand_genus[acc]
        summary[g]["candidates"] += 1
        if acc in best:
            fid, tgt = best[acc]
            summary[g]["any_pair_cov_0.5"] += 1
            if fid >= 0.6:
                summary[g]["pair_T_C_0.6"] += 1
            if fid >= 0.7:
                summary[g]["pair_T_C_0.7"] += 1
            rows.append((acc, g, f"{fid:.3f}", tgt, pos_genus.get(tgt, "unknown")))
    for acc in sorted(pos_acc):
        summary[pos_genus[acc]]["positives"] += 1

    # S4 table, one row per genus that has positives or candidates of interest
    focus = [
        "Aspergillus",
        "Alternaria",
        "Candida",
        "Coccidioides",
        "Histoplasma",
        "Blastomyces",
        "Paracoccidioides",
        "Trichophyton",
        "Arthroderma",
        "Ajellomyces",
        "Emmonsia",
        "Malassezia",
        "Penicillium",
        "Curvularia",
    ]
    genera = sorted(set(focus) | set(summary))
    with open(os.path.join(args.outdir, "species_check.tsv"), "w") as fh:
        fh.write("genus\tpositives\tcandidates\tany_pair_cov_0.5\tpair_T_C_0.6\tpair_T_C_0.7\n")
        for g in genera:
            c = summary[g]
            fh.write(
                f"{g}\t{c['positives']}\t{c['candidates']}\t{c['any_pair_cov_0.5']}\t{c['pair_T_C_0.6']}\t{c['pair_T_C_0.7']}\n"
            )
    with open(os.path.join(args.outdir, "species_check.md"), "w") as fh:
        fh.write(
            "# Species check: positives, candidates, and candidates that pair with a known allergen\n\n"
        )
        fh.write(
            "Pairing uses the O7 coverage rule (query and target coverage >= 0.5). A candidate that pairs at T_C 0.6 or 0.7 "
            "may be an unannotated allergen (label noise, spec 3.2).\n\n"
        )
        fh.write(
            "| Genus | Positives | Candidates | Any pair | Pair T_C 0.6 | Pair T_C 0.7 |\n|---|---|---|---|---|---|\n"
        )
        for g in genera:
            c = summary[g]
            fh.write(
                f"| {g} | {c['positives']} | {c['candidates']} | {c['any_pair_cov_0.5']} | {c['pair_T_C_0.6']} | {c['pair_T_C_0.7']} |\n"
            )
        fh.write("\n## Candidates that pair with a positive (best identity per candidate)\n\n")
        fh.write(
            "| Candidate | Genus | Best fident | Positive | Positive genus |\n|---|---|---|---|---|\n"
        )
        for acc, g, fid, tgt, pg in sorted(rows, key=lambda r: -float(r[2])):
            fh.write(f"| {acc} | {g} | {fid} | {tgt} | {pg} |\n")
    print(f"candidates={len(cand_acc)} positives={len(pos_acc)} pairing_candidates={len(best)}")
    if fails:
        print(f"structural checks failed: {fails}")
        sys.exit(1)


if __name__ == "__main__":
    main()
