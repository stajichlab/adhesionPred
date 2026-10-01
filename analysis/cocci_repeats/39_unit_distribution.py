#!/usr/bin/python3.12
"""Where is the SOWgp repeat unit found? Tables from the 37 search.

Reads   unit_genomes.tsv (35), unit_hmm_hits.tsv.gz, unit_anchor_hits.tsv.gz (37),
        sowgp_unit_modal.fa, sowgp_seed.fa (34, 05).
Writes  unit_hits_by_protein.tsv     one row per protein with >= 1 HMM domain at E <= 1
        unit_copynumber.tsv          copy-number distribution, by group
and prints the tables quoted in REPORT_2026-09-30_sowgp_unit_distribution.md.

DEFINITIONS (all thresholds fixed here, none tuned to the result)
  unit hit       hmmsearch domain with i-evalue <= 1 at Z = 6e7 (36: 0 hits in 30,715 shuffled
                 proteins at E <= 1000; chance expectation 0.5).
  copy number    number of such domains in the protein, by HMM. For SOWgp proteins it is
                 compared with the anchored count (30/37) as a check on the HMM count.
  SOWgp protein  the N-terminal flank ref[:first PTDCYGDC] of SOWgp58 (86 aa) aligns to the
                 protein at >= 60% identity over >= 30 aa. The flank carries no repeat, so this
                 separates "SOWgp" from "another protein carrying the unit".
  unit identity  ungapped identity, over the aligned columns, of a BLOSUM62 local alignment of
                 the HMM-domain residues to the SOWgp modal unit (47 aa); also the fraction of
                 the 47 modal positions identical.
  tier           best-domain bit score; see the sensitivity table (36): 75 ~ 80% identity to the
                 modal unit, 50 ~ 60%, 40 ~ 50%, 26 ~ 40%.
"""

import csv
import gzip
import importlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_common"))
sys.path.insert(0, str(HERE))
import paths  # noqa: E402  -- from ../_common

fam = importlib.import_module("30_anchor_family_search")
det = importlib.import_module("14_repeat_detect_general")

E_HIT = 1.0
FLANK_ID, FLANK_AA = fam.MIN_IDENTITY, fam.MIN_ALN_LEN


def tier(score):
    return (
        "A>=75" if score >= 75 else "B50-75" if score >= 50 else "C40-50" if score >= 40 else "D<40"
    )


def main():
    genomes = {}
    for r in csv.DictReader(open(HERE / "unit_genomes.tsv"), delimiter="\t"):
        genomes[(r["source"], r["label"])] = r
    by_label = defaultdict(list)
    for (_src, lab), r in genomes.items():
        by_label[lab].append(r)

    modal = list(fam.read_fasta(HERE / "sowgp_unit_modal.fa").values())[0]
    sow = fam.read_fasta(HERE / "sowgp_seed.fa")["SOWgp58_Cocci_immitis_published"]
    flank = sow[: sow.index("PTDCYGDC")]
    print(f"SOWgp58 N-terminal flank: {len(flank)} aa")
    al = fam.make_aligner()

    # ---- domain hits
    dom = defaultdict(list)  # (label, protein) -> [(score, E, ali_from, ali_to, tlen)]
    n_all = 0
    with gzip.open(HERE / "unit_hmm_hits.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_all += 1
            if float(r["dom_ievalue"]) <= E_HIT:
                dom[(r["proteome"], r["protein"])].append(
                    (
                        float(r["dom_score"]),
                        float(r["dom_ievalue"]),
                        int(r["ali_from"]),
                        int(r["ali_to"]),
                        int(r["length"]),
                    )
                )
    print(
        f"{n_all} domain rows at E <= 1000; {sum(len(v) for v in dom.values())} at E <= {E_HIT} in "
        f"{len(dom)} proteins"
    )

    # anchored members, for the cross-check of copy number
    anch = {}
    with gzip.open(HERE / "unit_anchor_hits.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            anch[(r["proteome"], r["protein"])] = r

    # ---- fetch sequences for hit proteins, one proteome at a time
    want = defaultdict(set)
    for lab, p in dom:
        want[lab].add(p)
    try:
        import duckdb

        con = duckdb.connect(str(paths.FUNGI5K_DUCKDB), read_only=True)
    except Exception as e:  # noqa: BLE001
        con = None
        print("no duckdb, SignalP not joined:", e)

    rows = []
    for lab, prots in sorted(want.items()):
        recs = by_label.get(lab, [])
        if not recs:
            continue
        g = recs[0]
        seqs = fam.read_fasta(g["path"])
        for p in sorted(prots):
            s = seqs.get(p)
            if s is None:
                continue
            ds = sorted(dom[(lab, p)], key=lambda d: d[2])
            fl = fam.align_to_ref(al, s, flank)
            is_sow = fl["identity_pct"] >= FLANK_ID and fl["aln_len"] >= FLANK_AA
            ids47, idal = [], []
            for _sc, _e, a, b, _ in ds:
                seg = s[max(a - 1, 0) : b]
                aln = al.align(seg, modal)[0]
                qb, rb = aln.aligned
                ident = sum(
                    1
                    for (q0, q1), (r0, r1) in zip(qb, rb, strict=False)
                    for i in range(q1 - q0)
                    if seg[q0 + i] == modal[r0 + i]
                )
                alen = sum(int(q1 - q0) for q0, q1 in qb)
                ids47.append(ident / len(modal))
                idal.append(ident / max(alen, 1))
            comp = det.composition(s)
            best = max(d[0] for d in ds)
            rows.append(
                {
                    "proteome": lab,
                    "source": g["source"],
                    "species": g["species"],
                    "genus": g["genus"],
                    "order": g["order"],
                    "protein": p,
                    "length": len(s),
                    "n_units_hmm": len(ds),
                    "best_score": round(best, 1),
                    "tier": tier(best),
                    "median_unit_id47": round(float(np.median(ids47)), 3),
                    "median_unit_id_aligned": round(float(np.median(idal)), 3),
                    "flank_id": fl["identity_pct"],
                    "flank_aln": fl["aln_len"],
                    "sowgp": "yes" if is_sow else "no",
                    "n_anchor": anch.get((lab, p), {}).get("n_anchor", 0),
                    "pct_pro": comp["pct_pro"],
                    "pct_cys": comp["pct_cys"],
                    "pct_ser_thr": comp["pct_ser_thr"],
                    "pct_charged": comp["pct_charged"],
                    "scores": ",".join(f"{d[0]:.0f}" for d in ds),
                    "starts": ",".join(str(d[2]) for d in ds),
                }
            )
    if con is not None:
        ids = tuple({r["protein"] for r in rows if r["source"] == "fungi5k"})
        if ids:
            q = ",".join("'%s'" % i for i in ids)
            sp = dict(
                con.execute(
                    f"SELECT protein_id, probability FROM signalp WHERE protein_id IN ({q})"
                ).fetchall()
            )
            for r in rows:
                r["signalp_prob"] = sp.get(r["protein"], "") if r["source"] == "fungi5k" else "n/a"
    with open(HERE / "unit_hits_by_protein.tsv", "w") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=list(rows[0]) + (["signalp_prob"] if "signalp_prob" not in rows[0] else []),
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} proteins with >= 1 unit domain at E <= {E_HIT}\n")

    # ---- summaries
    def grp(r):
        if r["source"] in ("cocci_longread", "cocci_pangenome"):
            return "Coccidioides (7 long-read + 493 pangenome)"
        if r["genus"] == "Coccidioides":
            return "Coccidioides (Fungi_5k RS/Silveira; duplicates of the 7 long-read)"
        if r["order"] == "Onygenales":
            return "Onygenales, not Coccidioides"
        return "other Fungi_5k: " + (r["order"] or "unassigned")

    print("== SOWgp vs other proteins, by group ==")
    tab = defaultdict(Counter)
    for r in rows:
        tab[grp(r)]["sowgp" if r["sowgp"] == "yes" else "other"] += 1
    for k in sorted(tab):
        print(f"{k:<70}{tab[k]['sowgp']:>6}{tab[k]['other']:>6}")

    print("\n== copy number of SOWgp proteins: HMM count vs anchored count ==")
    cc = Counter(
        (r["n_units_hmm"], int(r["n_anchor"] or 0))
        for r in rows
        if r["sowgp"] == "yes" and r["source"] == "cocci_pangenome"
    )
    for (h, a), n in sorted(cc.items()):
        print(f"  HMM units {h}  anchored {a}  proteins {n}")
    sow_rows = [r for r in rows if r["sowgp"] == "yes"]
    print(
        f"  SOWgp proteins total {len(sow_rows)}; HMM copy number distribution:",
        dict(sorted(Counter(r["n_units_hmm"] for r in sow_rows).items())),
    )

    oth = [r for r in rows if r["sowgp"] == "no"]
    print(f"\n== {len(oth)} non-SOWgp proteins with a unit domain at E <= {E_HIT} ==")
    print("by tier and group:")
    tt = defaultdict(Counter)
    for r in oth:
        tt[grp(r)][r["tier"]] += 1
    for k in sorted(tt):
        print(
            f"  {k:<70}"
            + "  ".join(f"{t}:{tt[k][t]}" for t in ("A>=75", "B50-75", "C40-50", "D<40"))
        )
    print("\nnon-SOWgp hits, best score first (top 60):")
    hdr = f"{'genus':<18}{'species':<34}{'protein':<24}{'len':>5}{'n':>3}{'best':>6}{'id47':>6}{'P%':>5}{'C%':>5}{'SP':>6}  order"
    print(hdr)
    for r in sorted(oth, key=lambda r: -r["best_score"])[:60]:
        print(
            f"{r['genus'][:17]:<18}{r['species'][:33]:<34}{r['protein'][:23]:<24}{r['length']:>5}"
            f"{r['n_units_hmm']:>3}{r['best_score']:>6}{r['median_unit_id47']:>6}{r['pct_pro']:>5}"
            f"{r['pct_cys']:>5}{str(r.get('signalp_prob',''))[:5]:>6}  {r['order']}"
        )
    with open(HERE / "unit_copynumber.tsv", "w") as fh:
        fh.write("group\tsowgp\tn_units_hmm\tproteins\n")
        c = Counter((grp(r), r["sowgp"], r["n_units_hmm"]) for r in rows)
        for (g, s, n), v in sorted(c.items()):
            fh.write(f"{g}\t{s}\t{n}\t{v}\n")


if __name__ == "__main__":
    main()
