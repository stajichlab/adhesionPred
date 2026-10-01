#!/usr/bin/python3.12
"""Positive control for the SOWgp unit search: BAD1 of Blastomyces dermatitidis.

BAD1 (UniProt A4D962, 1,146 aa, 24 aa Cys/Trp/His/Asp-rich tandem repeat) is a known
repeat surface adhesin. If the search is calibrated it must show, for BAD1 against SOWgp:
  (i)  ARCHITECTURE: a tandem repeat by the periodicity detector (14), with a composition
       comparable to SOWgp's;
  (ii) NO UNIT: no hit of the SOWgp unit HMM, no KKYGDC/PTDCYGDC anchor, and a local
       alignment of its repeat unit to the SOWgp modal unit that does not beat shuffled
       controls.
It then lists BAD1-like proteins in the Onygenales proteomes that the architecture scan (37,
unit_architecture_onygenales.tsv.gz) called, aligned to BAD1 by the same >=60% identity over
>=30 aa rule as 30_anchor_family_search.py.

Needs hmmer/3.4 on PATH. Run:  sbatch 38_unit_bad1_control.sh    (prints to the log)
"""

import csv
import gzip
import importlib
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
det = importlib.import_module("14_repeat_detect_general")
fam = importlib.import_module("30_anchor_family_search")


def local_score(al, a, b):
    return al.score(a, b)


def main():
    bad1 = list(fam.read_fasta(HERE / "bad1_A4D962.fa").values())[0]
    modal = list(fam.read_fasta(HERE / "sowgp_unit_modal.fa").values())[0]
    sow = fam.read_fasta(HERE / "sowgp_seed.fa")["SOWgp58_Cocci_immitis_published"]

    print("== (i) architecture: detector 14 + composition ==")
    print(
        f"{'protein':<22}{'len':>6}{'period':>8}{'copies':>8}{'cov':>7}{'z_reg':>7}"
        f"{'%P':>6}{'%C':>6}{'%S+T':>6}{'%chg':>6}"
    )
    for name, s in (("BAD1 A4D962", bad1), ("SOWgp58", sow)):
        d = det.detect(s)
        c = det.composition(s)
        print(
            f"{name:<22}{len(s):>6}{d['period']:>8}{d['n_copies']:>8}{d['coverage']:>7}"
            f"{d['z_region']:>7}{c['pct_pro']:>6}{c['pct_cys']:>6}{c['pct_ser_thr']:>6}"
            f"{c['pct_charged']:>6}"
        )
        if name.startswith("BAD1"):
            bad1_unit = d["unit"]
            print(f"  BAD1 detected unit ({len(bad1_unit)} aa): {bad1_unit}")

    print("\n== (ii) unit-level homology to SOWgp ==")
    # anchors
    for m in ("KKYGDC", "PTDCYGDC"):
        print(f"anchor {m}: {bad1.count(m)} occurrences in BAD1")
    # HMM
    scratch = Path(os.environ.get("SCRATCH", "/tmp")) / "bad1ctl"
    scratch.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "hmmsearch",
            "-Z",
            "60000000",
            "--domZ",
            "60000000",
            "-E",
            "1000",
            "--domE",
            "1000",
            "--noali",
            "--domtblout",
            str(scratch / "b.dom"),
            "-o",
            os.devnull,
            str(HERE / "sowgp_unit.hmm"),
            str(HERE / "bad1_A4D962.fa"),
        ],
        check=True,
    )
    hits = [ln.split() for ln in open(scratch / "b.dom") if not ln.startswith("#")]
    print(
        f"unit HMM on BAD1: {len(hits)} domain hits at E <= 1000 (Z = 6e7)"
        + (f"; best score {max(float(h[13]) for h in hits):.1f}" if hits else "")
    )
    # unit vs modal unit, local alignment, against shuffled modal units
    al = fam.make_aligner()
    rng = random.Random(20260930)
    obs = al.score(bad1_unit, modal)
    # the detector's unit is a phase-arbitrary 24-mer; try every rotation of it
    rot = max(al.score(bad1_unit[i:] + bad1_unit[:i], modal) for i in range(len(bad1_unit)))
    null = []
    for _ in range(2000):
        ls = list(modal)
        rng.shuffle(ls)
        null.append(al.score(bad1_unit, "".join(ls)))
    null = np.array(null)
    print(
        f"BAD1 unit vs SOWgp modal unit, BLOSUM62 local score {obs:.0f} "
        f"(best rotation {rot:.0f}); 2000 shuffles of the SOWgp unit: mean {null.mean():.1f}, "
        f"max {null.max():.0f}; z = {(obs - null.mean()) / null.std():.2f}"
    )
    a = al.align(bad1_unit, modal)[0]
    print(f"alignment of the two units:\n{a}")
    # whole-protein alignment identity, for comparison with the 60%/30 aa membership rule
    r = fam.align_to_ref(al, bad1, sow)
    print(
        f"BAD1 vs SOWgp58 whole-protein local alignment: {r['identity_pct']}% identity over "
        f"{r['aln_len']} aa (membership rule: >= 60% over >= 30 aa)"
    )

    print("\n== BAD1-like proteins in Onygenales (architecture scan, 14) ==")
    genomes = {
        r["label"]: r for r in csv.DictReader(open(HERE / "unit_genomes.tsv"), delimiter="\t")
    }
    arch = HERE / "unit_architecture_onygenales.tsv.gz"
    rows = [
        r
        for r in csv.DictReader(gzip.open(arch, "rt"), delimiter="\t")
        if r["period"] and int(r["period"]) > 0
    ]
    print(f"{len(rows)} Onygenales proteins with a detected tandem repeat (all periods)")
    want = {}
    for r in rows:
        if 20 <= int(r["period"]) <= 28:
            want.setdefault(r["strain"], set()).add(r["protein"])
    n_aln = 0
    out = []
    for lab, prots in sorted(want.items()):
        g = genomes.get(lab)
        if not g:
            continue
        seqs = fam.read_fasta(g["path"])
        for p in prots:
            s = seqs.get(p)
            if not s:
                continue
            res = fam.align_to_ref(al, s, bad1)
            n_aln += 1
            if res["identity_pct"] >= fam.MIN_IDENTITY and res["aln_len"] >= fam.MIN_ALN_LEN:
                rr = next(x for x in rows if x["strain"] == lab and x["protein"] == p)
                out.append(
                    (
                        lab,
                        g["species"],
                        p,
                        len(s),
                        rr["period"],
                        rr["n_copies"],
                        res["identity_pct"],
                        res["aln_len"],
                        res["ref_cov"],
                    )
                )
    print(f"{n_aln} proteins with period 20-28 aligned to BAD1")
    print(
        f"{'proteome':<40}{'protein':<24}{'len':>6}{'per':>5}{'copies':>7}{'id%':>6}{'aln':>6}{'BAD1cov':>8}"
    )
    for o in out:
        print(f"{o[0]:<40}{o[2]:<24}{o[3]:>6}{o[4]:>5}{o[5]:>7}{o[6]:>6}{o[7]:>6}{o[8]:>8}")
    with open(HERE / "unit_bad1_like.tsv", "w") as fh:
        fh.write(
            "proteome\tspecies\tprotein\tlength\tperiod\tcopies\tidentity_pct\taln_len\tbad1_cov\n"
        )
        for o in out:
            fh.write("\t".join(str(x) for x in o) + "\n")


if __name__ == "__main__":
    main()
