#!/usr/bin/env python
"""Class 2a retrain test (docs/HANDOFF-HPCC.md §4): does adding the 41 curated
Coccidioides class2a candidates as training positives improve cross-clade
repeat-structure recovery, specifically for SOWgp58/SOWgp66?

Baseline (repeat_structure_transfer.py, unchanged): trains on Saccharomycotina
adhesins only, recovers 3/5 of {SOWgp58, SOWgp66, SOWgp82, BAD1, CspA} with
repeat-structure features. SOWgp58 (p=0.334) and SOWgp66 (p=0.450) are the
misses closest to the 0.5 threshold.

This script adds analysis/cocci_repeats/class2a_candidates.tsv sequences as
extra training positives, in two variants:
  - "+all_candidates": all 41 candidates minus any that ARE a TEST protein
  - "+procys_only": only the Pro/Cys-rich subset minus TEST proteins

Held-out-integrity check (the caution from the handoff): candidate sequences
are aligned against all 5 TEST proteins (Biopython pairwise local alignment,
BLOSUM62); any candidate >60% identity to a TEST protein is EXCLUDED from
training as the same protein/ortholog, not a genuinely independent example.
This is done by direct alignment, not by trusting the locus name or the
candidate table's own comp_class label (one excluded candidate, Silveira
QVM10648.1, was labeled "other" but is 63.6% identical to SOWgp58).

Caveat carried forward from the handoff: these 41 candidates are
computational predictions (repeat detector + composition profile) with no
functional validation. A null result here is as reportable as a positive one.
"""

import csv
import sys
import urllib.parse
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "analysis" / "curation"))
from curation_lib import get  # noqa: E402

sys.path.insert(0, str(REPO / "analysis" / "model_review"))
from repeat_structure_transfer import (  # noqa: E402
    SACCH,
    TEST,
    comp_desc,
    repeat_desc,
)

# Default: the 41 candidates from 02_repeat_profile.py + 03 (the 2026-09-28 run).
# Override with --candidates to test a different candidate set, e.g. the 58 produced by
# 14_repeat_detect_general.py + 03 after the 2026-09-30 SignalP re-run. The file must have
# the same columns (strain, protein, comp_class).
CANDIDATES_TSV = REPO / "analysis" / "cocci_repeats" / "class2a_candidates.tsv"
STRAIN_FASTA = {
    "CimmitisRS_FungiDB": "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input_run2/CimmitisRS_FungiDB.fasta",
    "CposadasiiSilveira2022_FungiDB": "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input_run2/CposadasiiSilveira2022_FungiDB.fasta",
    "Coccidioides_immitis_CiB10637": "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc/CiB10637/Coccidioides_immitis_CiB10637.proteins.fa",
    "Coccidioides_immitis_CiB10992": "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc/CiB10992/Coccidioides_immitis_CiB10992.proteins.fa",
    "Coccidioides_immitis_VFC140": "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc/VFC140/Coccidioides_immitis_VFC140.proteins.fa",
    "Coccidioides_posadasii_Cpos1038": "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc/Cpos1038/Coccidioides_posadasii_Cpos1038.proteins.fa",
    "Coccidioides_posadasii_Cpos3700": "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc/Cpos3700/Coccidioides_posadasii_Cpos3700.proteins.fa",
}
IDENTITY_EXCLUDE_THRESHOLD = 60.0


def read_fasta(path):
    d, hdr, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if hdr is not None:
                d[hdr] = "".join(buf)
            hdr = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if hdr is not None:
        d[hdr] = "".join(buf)
    return d


def load_candidates():
    rows = list(csv.DictReader(open(CANDIDATES_TSV), delimiter="\t"))
    by_strain = {}
    for r in rows:
        by_strain.setdefault(r["strain"], []).append(r["protein"])
    seqs, comp_class = {}, {}
    for strain in by_strain:
        strain_seqs = read_fasta(STRAIN_FASTA[strain])
        for r in rows:
            if r["strain"] != strain:
                continue
            key = f"{strain}|{r['protein']}"
            seqs[key] = strain_seqs[r["protein"]]
            comp_class[key] = r["comp_class"]
    return seqs, comp_class


def fetch_test_seqs():
    accs = sorted(TEST)
    q = urllib.parse.urlencode({"accessions": ",".join(accs), "format": "fasta"})
    seqs = {}
    for blk in get(f"https://rest.uniprot.org/uniprotkb/accessions?{q}", "text/plain").split(">"):
        if blk.strip():
            h, *rest = blk.split("\n")
            seqs[h.split("|")[1]] = "".join(rest)
    return seqs


def pct_identity(aligner, seq1, seq2):
    aln = aligner.align(seq1, seq2)[0]
    a1, a2 = aln.aligned
    matches = alnlen = 0
    for (s1, e1), (s2, e2) in zip(a1, a2):
        x, y = seq1[s1:e1], seq2[s2:e2]
        alnlen += len(x)
        matches += sum(1 for p, q in zip(x, y) if p == q)
    return matches / alnlen * 100 if alnlen else 0.0


def exclude_test_orthologs(cand_seqs, test_seqs):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"
    keep, excluded = {}, []
    for key, seq in cand_seqs.items():
        best_pid = max(pct_identity(aligner, seq, tseq) for tseq in test_seqs.values())
        if best_pid > IDENTITY_EXCLUDE_THRESHOLD:
            excluded.append((key, best_pid))
        else:
            keep[key] = seq
    return keep, excluded


def evaluate(fn, keys, train_accs, seqs, y_extra_keys, y):
    X_train = np.array([[fn(seqs[a]).get(k, 0.0) for k in keys] for a in train_accs])
    m = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")
    ).fit(X_train, y)
    found = 0
    rows = []
    for acc, name in TEST.items():
        if acc not in seqs:
            rows.append((name, None))
            continue
        p = m.predict_proba(np.array([[fn(seqs[acc]).get(k, 0.0) for k in keys]]))[0, 1]
        found += p > 0.5
        rows.append((name, p))
    return found, rows


def main():
    import argparse

    global CANDIDATES_TSV
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--candidates",
        default=None,
        help="candidate TSV to use instead of the default 41-candidate file",
    )
    args = ap.parse_args()
    if args.candidates:
        CANDIDATES_TSV = Path(args.candidates)
    print(f"candidate table: {CANDIDATES_TSV}")

    with open(REPO / "data" / "curated" / "surface" / "surface.tsv") as f:
        meta = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
    lab = {
        a: (1 if r["adhesion_status"] == "adhesin" else 0)
        for a, r in meta.items()
        if r["adhesion_status"] in ("adhesin", "non_adhesin")
    }
    sacch_train = [a for a in lab if any(k in meta[a]["genome"] for k in SACCH)]

    test_seqs = fetch_test_seqs()
    q = urllib.parse.urlencode({"accessions": ",".join(sacch_train), "format": "fasta"})
    sacch_seqs = {}
    for blk in get(f"https://rest.uniprot.org/uniprotkb/accessions?{q}", "text/plain").split(">"):
        if blk.strip():
            h, *rest = blk.split("\n")
            sacch_seqs[h.split("|")[1]] = "".join(rest)
    sacch_train = [a for a in sacch_train if a in sacch_seqs]

    cand_seqs, comp_class = load_candidates()
    print(f"class2a candidates: {len(cand_seqs)} total")
    kept, excluded = exclude_test_orthologs(cand_seqs, test_seqs)
    print(f"  excluded as TEST-protein orthologs (>{IDENTITY_EXCLUDE_THRESHOLD}% identity):")
    for key, pid in excluded:
        print(f"    {key}  ({pid:.1f}%)")
    print(f"  {len(kept)} remain as candidate extra positives\n")

    procys_only = {k: v for k, v in kept.items() if "ProCys_rich" in comp_class[k]}
    print(f"  of which {len(procys_only)} are Pro/Cys-rich (the architecture SOWgp/BAD1 need)\n")

    all_seqs = {**sacch_seqs, **kept, **test_seqs}

    variants = {
        "baseline (Saccharomycotina only, unchanged)": [],
        # Counts are computed, not hardcoded. They were literal "36"/"16" until
        # 2026-09-30, which silently mislabelled every run on a different candidate set.
        f"+all_candidates ({len(kept)} non-SOWgp class2a candidates)": list(kept),
        f"+procys_only ({len(procys_only)} non-SOWgp Pro/Cys-rich candidates)": list(procys_only),
    }

    for variant_name, extra_keys in variants.items():
        train_accs = sacch_train + extra_keys
        y = np.array([lab[a] for a in sacch_train] + [1] * len(extra_keys))
        print(f"=== {variant_name} ===")
        print(f"train: {int(y.sum())} positives, {int((y == 0).sum())} negatives")
        for feat_name, fn in [
            ("composition (20 aa freq + length)", comp_desc),
            ("repeat STRUCTURE (composition-agnostic)", repeat_desc),
        ]:
            keys = sorted(fn("ACDEFGHIKLMNPQRSTVWY" * 10))
            found, rows = evaluate(fn, keys, train_accs, all_seqs, extra_keys, y)
            print(f"  {feat_name}:")
            for name, p in rows:
                if p is None:
                    print(f"     {name:<10} (sequence not retrieved)")
                else:
                    print(f"     {name:<10} p={p:.3f}  {'FOUND' if p > 0.5 else 'missed'}")
            print(f"     -> {found}/{len(TEST)} recovered")
        print()


if __name__ == "__main__":
    main()
