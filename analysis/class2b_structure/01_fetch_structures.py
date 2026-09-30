#!/usr/bin/env python3.12
"""Fetch AlphaFold DB models for the class 2b candidate set.

Class 2b in docs/TOOL-ARCHITECTURE.md is "small receptor-binding invasins" (CalA, Ag2/PRA,
PRA3). Every member already has an AlphaFold DB model, so no structure prediction is needed.
This script downloads the v6 mmCIF for each accession in candidates.tsv.

Run on the login node (network only, no compute):
    /usr/bin/python3.12 01_fetch_structures.py --outdir <workdir>/structures
"""

import argparse
import csv
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
AFDB = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def fetch(url, dest, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                dest.write_bytes(r.read())
            return True
        except urllib.error.HTTPError as e:
            print(f"  HTTP {e.code} for {url}", file=sys.stderr)
            return False
        except Exception as e:  # transient network
            print(f"  attempt {i + 1} failed: {e}", file=sys.stderr)
            time.sleep(5)
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--candidates", default=str(HERE / "candidates.tsv"))
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "cif").mkdir(parents=True, exist_ok=True)
    (outdir / "fasta").mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    n_ok = 0
    for row in rows:
        acc, label = row["accession"], row["label"]
        cif = outdir / "cif" / f"{label}_{acc}.cif"
        faa = outdir / "fasta" / f"{label}_{acc}.fasta"
        if not cif.exists():
            print(f"{label} ({acc}): fetching model")
            if not fetch(AFDB.format(acc=acc), cif):
                print(f"  NO AlphaFold DB model for {acc}", file=sys.stderr)
                cif.unlink(missing_ok=True)
                continue
        if not faa.exists():
            fetch(UNIPROT_FASTA.format(acc=acc), faa)
        n_ok += 1
    print(f"\n{n_ok}/{len(rows)} models present in {outdir / 'cif'}")


if __name__ == "__main__":
    main()
