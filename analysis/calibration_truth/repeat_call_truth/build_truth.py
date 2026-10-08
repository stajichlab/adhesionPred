#!/usr/bin/env python3
"""Build per-species truth tables for the call `tandem_repeat_protein` ("has a tandem repeat region").

Input: repeat_truth_candidates.tsv (label 1, 0 or empty per curated protein; see make_candidates.py)
and the proteome FASTA of each run (to cluster the truth proteins at 30% identity, coverage 0.5, as
analysis/calibration_truth/c1_truth_count.py does). Output per species: truth.<species>.tsv with the
columns id, label, cluster (what `cellsurface_sorting_hat_calibrate truth` reads) and
truth.<species>.annotated.tsv with the basis of each label.

Labels: 1 = paper statement or at least two UniProt Repeat features. 0 = no UniProt Repeat feature
(absence of annotation, so an ASSUMED negative: UniProt lacks features for some repeat proteins).
Proteins with one feature, a family-inference-only claim, or no record are left out.

Usage: build_truth.py --candidates FILE --fasta Scer_S288C=PATH --fasta Calb_SC5314=PATH --out DIR
Needs `mmseqs` on PATH (module load MMseqs2/17-b804f).
"""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).rstrip("*")
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).rstrip("*")
    return seqs


def cluster(seqs, threads=4):
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "in.faa"
        with open(fa, "w") as f:
            for k, s in seqs.items():
                f.write(f">{k}\n{s}\n")
        subprocess.run(
            [
                "mmseqs",
                "easy-cluster",
                str(fa),
                f"{tmp}/clu",
                f"{tmp}/tmp",
                "--min-seq-id",
                "0.3",
                "-c",
                "0.5",
                "-v",
                "1",
                "--threads",
                str(threads),
            ],
            check=True,
            capture_output=True,
        )
        rep = {}
        for line in open(f"{tmp}/clu_cluster.tsv"):
            a, b = line.rstrip("\n").split("\t")
            rep[b] = a
    return rep


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--fasta", action="append", required=True, help="PROTEOME=PATH, repeatable")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.candidates), delimiter="\t"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for spec in a.fasta:
        prot, path = spec.split("=", 1)
        seqs = read_fasta(path)
        keep = [r for r in rows if r["proteome"] == prot and r["label"] in ("0", "1")]
        missing = [r["protein"] for r in keep if r["protein"] not in seqs]
        if missing:
            raise SystemExit(f"{prot}: {len(missing)} proteins not in the FASTA, e.g. {missing[0]}")
        rep = cluster({r["protein"]: seqs[r["protein"]] for r in keep})
        with open(out / f"truth.{prot}.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["id", "label", "cluster"])
            for r in keep:
                w.writerow([r["protein"], r["label"], rep[r["protein"]]])
        with open(out / f"truth.{prot}.annotated.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(
                [
                    "id",
                    "gene",
                    "accession",
                    "label",
                    "cluster",
                    "basis",
                    "curated_class",
                    "evidence",
                    "tuned_or_homolog",
                ]
            )
            for r in keep:
                w.writerow(
                    [
                        r["protein"],
                        r["gene"],
                        r["accession"],
                        r["label"],
                        rep[r["protein"]],
                        r["basis"],
                        r["curated_class"],
                        r["evidence"],
                        r["tuned_or_homolog"],
                    ]
                )
        pos = [r for r in keep if r["label"] == "1"]
        neg = [r for r in keep if r["label"] == "0"]
        print(
            f"{prot}: {len(pos)} positives in {len({rep[r['protein']] for r in pos})} clusters; "
            f"{len(neg)} negatives in {len({rep[r['protein']] for r in neg})} clusters",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
