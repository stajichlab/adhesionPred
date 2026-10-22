#!/usr/bin/env python3
"""Fetch the Pfam seed alignments of the Task 5.3 models from the InterPro API (LOCAL, needs network).

DRAFT, not run. Not for HPCC compute nodes (they may have no network): run it where the network is
open, then copy `seeds.faa.gz` and its sha256 to the HPCC input directory of cluster_split.sh.

Models = the 21 wall and fusion models of analysis/wall_fusion_report/wall_fusion_models.tsv plus the
active adhesion family rows of data/sorting_hat/family_table.tsv (column `active` = yes). The union
is the seed reference for the `seed-like` / `seed-independent` label (analysis/lectin_validation
used the 20 lectin models for its own set; this one follows the models of the Kiss slices).

Header of the FASTA: `<PFxxxxx>|<seed sequence name>` (the form that seed_split.py reads).
InterPro serves the CURRENT Pfam release. The summary records the version of each alignment; compare
it with the Pfam 38.2 file used for the scan (docs/superpowers/plans: Pfam 38.2, HMMER 3.4) and report
a difference.

Usage: python3 -I analysis/wall_fusion_truth/fetch_seeds_models.py --out-dir DIR [--list-only]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WALL_MODELS = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
FAMILY_TABLE = ROOT / "data/sorting_hat/family_table.tsv"
URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=alignment:seed&download"


def model_list(wall: Path = WALL_MODELS, family: Path = FAMILY_TABLE) -> list[str]:
    with wall.open(encoding="utf-8-sig", newline="") as fh:
        ids = [r["id"] for r in csv.DictReader(fh, delimiter="\t")]
    with family.open(encoding="utf-8-sig", newline="") as fh:
        ids += [
            r["pfam_acc"]
            for r in csv.DictReader(fh, delimiter="\t")
            if r["active"] == "yes" and r["pfam_acc"]
        ]
    return sorted(dict.fromkeys(ids))


def parse_stockholm(text: str) -> tuple[str, str, dict[str, str]]:
    """(version AC, name ID, {sequence name: ungapped upper-case sequence}) of a Stockholm seed file."""
    ver = name = ""
    seqs: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith("#=GF AC"):
            ver = line.split()[2]
        elif line.startswith("#=GF ID"):
            name = line.split()[2]
        elif line and not line.startswith("#") and not line.startswith("//"):
            f = line.split()
            if len(f) == 2:
                seqs[f[0]] = seqs.get(f[0], "") + f[1]
    clean = {k: "".join(c for c in s if c.isalpha()).upper() for k, s in seqs.items()}
    return ver, name, clean


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--list-only", action="store_true", help="print the model list and stop")
    a = ap.parse_args(argv)
    models = model_list()
    if a.list_only:
        print("\n".join(models))
        return 0
    a.out_dir.mkdir(parents=True, exist_ok=True)
    summary = ["accession\tversion\tname\tn_seed_sequences"]
    with gzip.open(a.out_dir / "seeds.faa.gz", "wt") as fa:
        for acc in models:
            for _attempt in range(3):
                try:
                    req = urllib.request.Request(
                        URL.format(acc=acc), headers={"User-Agent": "adhesionPred-kiss"}
                    )
                    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310 (https only)
                        txt = gzip.open(io.BytesIO(r.read()), "rt").read()
                    break
                except Exception as exc:  # network errors are retried, then fatal
                    print("retry", acc, exc, file=sys.stderr)
                    time.sleep(3)
            else:
                raise SystemExit(f"failed {acc}")
            ver, name, seqs = parse_stockholm(txt)
            for k in sorted(seqs):
                fa.write(f">{acc}|{k}\n{seqs[k]}\n")
            summary.append(f"{acc}\t{ver}\t{name}\t{len(seqs)}")
            print(acc, ver, name, len(seqs), file=sys.stderr)
    (a.out_dir / "seeds_summary.tsv").write_text("\n".join(summary) + "\n")
    p = a.out_dir / "seeds.faa.gz"
    (a.out_dir / "seeds.faa.gz.sha256").write_text(
        hashlib.sha256(p.read_bytes()).hexdigest() + "\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
