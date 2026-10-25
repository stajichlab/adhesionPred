"""Truth tables (id, label, cluster) for the status files of the AGA1 calls, from the labelled 850 proteins (stdlib only; accessions only, no sequence).

Two truths, because a call is measured against the property it claims:
  class   : label 1 = tier A positive; label 0 = hard negative or background assumed negative. The other roles (tier B, excluded, other mechanism,
            Kiss rows, uncurated repeat-bearing) are not truth. Used for `gpi_long_disorder_protein`: its status says how well the property combination
            agrees with the AGA1 class labels, NOT that every GPI protein with a long disordered run is in the class.
  uniprot_gpi : label 1 / 0 = the UniProt GPI-anchor feature is present / absent (a protein without the column value is dropped); all roles.
            Used for `gpi_anchored_protein`. Absence of a UniProt feature is an assumed negative.
The cluster is the MMseqs cluster of the scoping study (a protein without one is its own cluster, SINGLETON_<accession>). Only proteins that are in the
proteome of the species (given by its proteins.tsv of stage_proteomes.py) are written.

  python3 -I analysis/sorting_hat_run/aga1_status_truth.py --scoping PROTEINS_TSV --proteome-proteins PROTEINS_TSV2 --species SPECIES --out-class F --out-gpi F

Refuses an existing output. FATAL: and exit 1 if a table would have no positive or no negative.
"""

import argparse
import csv
import sys
from pathlib import Path

NEG_ROLES = ("hard_negative", "background_assumed_negative")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cluster_of(r):
    return r["cluster_id"].strip() or f"SINGLETON_{r['accession']}"


def class_truth(rows):
    out = []
    for r in rows:
        if r["role"] == "class_positive" and r["pos_tier"] == "A":
            out.append((r["accession"], 1, cluster_of(r)))
        elif r["role"] in NEG_ROLES:
            out.append((r["accession"], 0, cluster_of(r)))
    return out


def gpi_truth(rows):
    out = []
    for r in rows:
        v = r["uniprot_gpi_feature"].strip().lower()
        if v in ("yes", "no"):
            out.append((r["accession"], 1 if v == "yes" else 0, cluster_of(r)))
    return out


def write(path, rows):
    if Path(path).exists():
        raise SystemExit(f"FATAL: {path} exists; nothing is overwritten")
    labels = {lab for _, lab, _ in rows}
    if labels != {0, 1}:
        raise SystemExit(
            f"FATAL: {path} would have labels {sorted(labels)}; it needs positives and negatives"
        )
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "label", "cluster"])
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scoping", required=True)
    ap.add_argument("--proteome-proteins", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out-class", required=True)
    ap.add_argument("--out-gpi", required=True)
    a = ap.parse_args(argv)
    members = {r["accession"] for r in read_tsv(a.proteome_proteins) if r["species"] == a.species}
    if not members:
        print(
            f"FATAL: no proteins of species {a.species!r} in {a.proteome_proteins}", file=sys.stderr
        )
        return 1
    rows = [r for r in read_tsv(a.scoping) if r["accession"] in members]
    for path in (a.out_class, a.out_gpi):
        if Path(path).exists():
            print(f"FATAL: {path} exists; nothing is overwritten", file=sys.stderr)
            return 1
    try:
        write(a.out_class, class_truth(rows))
        write(a.out_gpi, gpi_truth(rows))
    except SystemExit as e:
        print(e, file=sys.stderr)
        return 1
    print(
        f"{a.species}: {len(rows)} labelled proteins in the proteome; class truth and UniProt GPI truth written"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
