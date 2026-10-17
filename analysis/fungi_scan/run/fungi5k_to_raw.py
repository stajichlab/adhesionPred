#!/usr/bin/python3.12
"""Write raw SignalP 6 and TMHMM tables in the scan format from the precomputed Fungi_5k results (B2a of the discovery-pass spec).

B2a showed that these results agree with the scan run on 56 proteomes (31,977 SP proteins identical; 489,473 TMHMM rows identical).
SignalP: Fungi_5k keeps a GFF3 with the signal_peptide features only. Every protein of the FASTA that is not in the GFF3 gets
Prediction OTHER with OTHER probability 1 and SP probability 0 (the true SP probability is not stored; it is an evidence field,
not a call). Usage: fungi5k_to_raw.py --name NAME --fasta FASTA --workdir DIR
"""

import argparse
import gzip
import sys
from pathlib import Path

RES = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/results/function")


def fasta_ids(path):
    with open(path) as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def signalp_rows(ids, gff_lines):
    """Return (header lines, rows) in the prediction_results.txt layout."""
    sp = {}
    for line in gff_lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        sp[f[0].split()[0]] = (float(f[5]), int(f[4]))
    rows = []
    for i in ids:
        if i in sp:
            p, end = sp[i]
            rows.append(f"{i}\tSP\t{1 - p:.6f}\t{p:.6f}\tCS pos: {end}-{end + 1}.")
        else:
            rows.append(f"{i}\tOTHER\t1.000000\t0.000000\t")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--workdir", required=True)
    a = ap.parse_args()
    ids = fasta_ids(a.fasta)
    raw = Path(a.workdir) / "raw"
    (raw / "signalp").mkdir(parents=True, exist_ok=True)
    (raw / "tmhmm").mkdir(parents=True, exist_ok=True)
    with gzip.open(RES / "signalp" / f"{a.name}.signalp.gff3.gz", "rt") as fh:
        rows = signalp_rows(ids, fh)
    (raw / "signalp" / "prediction_results.txt").write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: from Fungi_5k results/function/signalp\n"
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n" + "\n".join(rows) + "\n"
    )
    (raw / "signalp" / "version.txt").write_text(
        "module=Fungi_5k results/function/signalp (B2a: identical to signalp/6-gpu on 56 proteomes)\n"
        "signalp6_version=SignalP 6.0 Signal peptide prediction tool 6.0h\n"
    )
    out = ["protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology"]
    with gzip.open(RES / "tmhmm" / f"{a.name}.tmhmm_short.tsv.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append("\t".join([f[0]] + [x.split("=", 1)[1] for x in f[1:6]]))
    (raw / "tmhmm" / "tmhmm.tsv").write_text("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
