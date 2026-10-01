#!/usr/bin/env python3.12
"""Summarise the class 2b structural results, with artifact checks.

Three tables:
  1. the all-vs-all TM-score matrix over the confident cores;
  2. the best PDB hit per protein (3Di E-value and TM-align TM-score);
  3. a composition check -- Cys, Pro and low-complexity content of each confident core.

The composition check exists because small disulfide-rich domains and low-complexity regions
both generate structural matches that are not fold homology. A hit is only called a fold
relationship here if TM-score >= 0.5 AND the aligned length covers most of the query core.

    /usr/bin/python3.12 05_summarize.py --workdir <workdir>
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent

ALLVALL_COLS = (
    "query target fident alnlen evalue bits alntmscore qtmscore ttmscore " "qlen tlen lddt prob"
).split()
PDB_COLS = (
    "query target fident alnlen evalue bits qtmscore ttmscore alntmscore "
    "qlen tlen lddt prob taxname"
).split()


def read_tsv(path, cols):
    if not Path(path).exists():
        return []
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < len(cols):
            continue
        out.append(dict(zip(cols, f)))
    return out


def core_sequence(per_res, label, segs):
    """Reconstruct the confident-core sequence from the per-residue table."""
    keep = set()
    for s in segs.split(";"):
        if s in ("none", ""):
            continue
        a, b = s.split("-")
        keep.update(range(int(a), int(b) + 1))
    return "".join(r["aa"] for r in per_res if r["label"] == label and int(r["resnum"]) in keep)


def lowcomplexity(seq, k=12):
    """Fraction of k-mer windows whose 3 most common residues cover >= 70% of the window."""
    if len(seq) < k:
        return 0.0
    hits = 0
    for i in range(len(seq) - k + 1):
        w = seq[i : i + k]
        top3 = sum(c for _, c in Counter(w).most_common(3))
        if top3 / k >= 0.7:
            hits += 1
    return hits / (len(seq) - k + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--candidates", default=str(HERE / "candidates.tsv"))
    args = ap.parse_args()
    wd = Path(args.workdir)

    cand = {r["label"]: r for r in csv.DictReader(open(args.candidates), delimiter="\t")}
    conf = {
        r["label"]: r for r in csv.DictReader(open(wd / "confidence_summary.tsv"), delimiter="\t")
    }
    per_res = list(csv.DictReader(open(wd / "plddt_per_residue.tsv"), delimiter="\t"))

    # ---- 1. all-vs-all -------------------------------------------------------
    av = read_tsv(wd / "foldseek" / "allvall.tsv", ALLVALL_COLS)
    m = {}
    for r in av:
        q = r["query"].replace(".pdb", "")
        t = r["target"].replace(".pdb", "")
        m[(q, t)] = r
    labs = [lbl for lbl in cand if lbl in conf and conf[lbl]["core_residues"] != "0"]

    print("## 1. All-vs-all TM-score over pLDDT>=70 cores (normalised by query)\n")
    print("%-9s" % "query", "".join("%-8s" % lab[:7] for lab in labs))
    for q in labs:
        cells = []
        for t in labs:
            r = m.get((q, t))
            cells.append("%-8s" % ("%.3f" % float(r["qtmscore"]) if r else "  .   "))
        print("%-9s" % q[:9], "".join(cells))

    print("\noff-diagonal pairs, TM-score >= 0.4:")
    seen = set()
    any_pair = False
    for (q, t), r in sorted(m.items(), key=lambda kv: -float(kv[1]["qtmscore"])):
        if q == t or (t, q) in seen:
            continue
        seen.add((q, t))
        tm = float(r["qtmscore"])
        if tm >= 0.4:
            any_pair = True
            print(
                f"  {q:10s} {t:10s} TM={tm:.3f}  E={float(r['evalue']):.2g}  "
                f"seqid={float(r['fident']):.2f}  alnlen={r['alnlen']}  "
                f"lddt={float(r['lddt']):.2f}"
            )
    if not any_pair:
        print("  none")

    # ---- 2. best PDB hit -----------------------------------------------------
    for fn, cols, tag in (
        ("vs_pdb_3di.tsv", PDB_COLS, "3Di+AA"),
        ("vs_pdb_tmalign.tsv", PDB_COLS, "TM-align"),
    ):
        hits = read_tsv(wd / "foldseek" / fn, cols)
        if not hits:
            print(f"\n## 2. PDB search ({tag}): no result file, skipped")
            continue
        print(f"\n## 2. Best PDB hits per protein ({tag})\n")
        print(
            "%-9s %-24s %-9s %-7s %-7s %-6s %s"
            % ("query", "top_pdb_hit", "evalue", "TMq", "seqid", "alnlen", "taxname")
        )
        best = {}
        for r in hits:
            q = r["query"].replace(".pdb", "")
            key = float(r["evalue"]) if tag == "3Di+AA" else -float(r["qtmscore"])
            if q not in best or key < best[q][0]:
                best[q] = (key, r)
        for q in labs:
            if q not in best:
                print("%-9s %-24s %s" % (q, "-", "no hit at the reported threshold"))
                continue
            r = best[q][1]
            print(
                "%-9s %-24s %-9.2g %-7.3f %-7.2f %-6s %s"
                % (
                    q,
                    r["target"][:24],
                    float(r["evalue"]),
                    float(r["qtmscore"]),
                    float(r["fident"]),
                    r["alnlen"],
                    r.get("taxname", "")[:30],
                )
            )

    # ---- 2b. targeted pairwise vs experimental references --------------------
    ref = read_tsv(wd / "foldseek" / "vs_ref_tmalign.tsv", ALLVALL_COLS)
    if ref:
        print(
            "\n## 2b. Targeted pairwise TM-align vs experimental PDB references"
            " (prefilter disabled)\n"
        )
        print(
            "%-9s %-22s %-7s %-6s %-7s %s"
            % ("query", "reference", "TMq", "LDDT", "seqid", "alnlen")
        )
        best = {}
        for r in ref:
            q = r["query"].replace(".pdb", "")
            pid = r["target"].split(".pdb")[0].split("_")[0]
            k = (q, pid)
            if k not in best or float(r["qtmscore"]) > float(best[k]["qtmscore"]):
                best[k] = r
        for q in labs:
            for (qq, _pid), r in sorted(best.items(), key=lambda kv: -float(kv[1]["qtmscore"])):
                if qq != q:
                    continue
                print(
                    "%-9s %-22s %-7.3f %-6.2f %-7.2f %s"
                    % (
                        q,
                        r["target"],
                        float(r["qtmscore"]),
                        float(r["lddt"]),
                        float(r["fident"]),
                        r["alnlen"],
                    )
                )

    # ---- 3. composition ------------------------------------------------------
    print("\n## 3. Composition of the confident core (artifact check)\n")
    print(
        "%-9s %-6s %-7s %-7s %-7s %-9s %s"
        % ("label", "core", "%Cys", "%Pro", "%Ser+Thr", "%lowcplx", "role")
    )
    for lab in labs:
        seq = core_sequence(per_res, lab, conf[lab]["core_segments"])
        n = max(len(seq), 1)
        c = Counter(seq)
        print(
            "%-9s %-6d %-7.1f %-7.1f %-7.1f %-9.2f %s"
            % (
                lab,
                len(seq),
                100 * c["C"] / n,
                100 * c["P"] / n,
                100 * (c["S"] + c["T"]) / n,
                lowcomplexity(seq),
                cand[lab]["role"],
            )
        )


if __name__ == "__main__":
    main()
