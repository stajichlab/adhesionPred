#!/usr/bin/python3.12
"""Map truth proteins (UniProt sequences) to proteome protein IDs by BLASTP (task H2b).

A truth protein maps to a proteome protein at >= 95% identity and >= 90% query coverage (best bitscore).
A proteome protein is used by one truth protein only; conflicts are listed as dropped.
Proteome IDs differ from UniProt IDs (Fungi_5k uses ...-T1), so no ID matching is done.

Usage: map_to_proteomes.py --truth truth_all.tsv --proteome NAME=FASTA [--proteome ...] --out DIR
Needs ncbi-blast+ on PATH (module load ncbi-blast/2.14.0+). Reads plain or .gz FASTA.
"""

import argparse
import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

FIELDS = ["qseqid", "sseqid", "pident", "qcovs", "bitscore"]


def filter_hits(hits, min_pident=95.0, min_qcov=90.0):
    keep = [h for h in hits if float(h["pident"]) >= min_pident and float(h["qcovs"]) >= min_qcov]
    dropped = [h for h in hits if h not in keep]
    best = {}
    for h in keep:
        q = h["qseqid"]
        if q not in best or float(h["bitscore"]) > float(best[q]["bitscore"]):
            best[q] = h
    mapped, taken = [], {}
    for h in sorted(best.values(), key=lambda r: -float(r["bitscore"])):
        s = h["sseqid"]
        if s in taken:
            dropped.append({**h, "reason": f"subject taken by {taken[s]}"})
        else:
            taken[s] = h["qseqid"]
            mapped.append(h)
    return mapped, dropped


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def run_blast(query_faa, proteome_faa, tmp):
    db = Path(tmp) / "db"
    plain = Path(tmp) / "proteome.faa"
    with _open(proteome_faa) as src, open(plain, "w") as dst:
        dst.write(src.read())
    subprocess.run(
        ["makeblastdb", "-in", str(plain), "-dbtype", "prot", "-out", str(db)],
        check=True,
        capture_output=True,
    )
    out = subprocess.run(
        [
            "blastp",
            "-query",
            str(query_faa),
            "-db",
            str(db),
            "-evalue",
            "1e-10",
            "-max_target_seqs",
            "5",
            "-outfmt",
            "6 " + " ".join(FIELDS),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return [dict(zip(FIELDS, line.split("\t"), strict=True)) for line in out.splitlines() if line]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--truth", required=True)
    ap.add_argument("--proteome", action="append", required=True, help="NAME=FASTA")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = [
        r for r in csv.DictReader(open(a.truth), delimiter="\t") if r["tier"] in ("T1", "T2", "T3")
    ]
    with tempfile.TemporaryDirectory(dir=None) as tmp:
        q = Path(tmp) / "truth.faa"
        with open(q, "w") as f:
            for r in rows:
                f.write(f">{r['accession']}\n{r['sequence']}\n")
        with (
            open(out / "proteome_map.tsv", "w", newline="") as fm,
            open(out / "proteome_map_dropped.tsv", "w", newline="") as fd,
        ):
            wm = csv.writer(fm, delimiter="\t")
            wm.writerow(
                ["proteome", "truth_accession", "proteome_id", "pident", "qcov", "bitscore"]
            )
            wd = csv.writer(fd, delimiter="\t")
            wd.writerow(["proteome", "truth_accession", "proteome_id", "pident", "qcov", "reason"])
            for spec in a.proteome:
                name, fasta = spec.split("=", 1)
                hits = run_blast(q, fasta, tmp)
                mapped, dropped = filter_hits(hits)
                for h in mapped:
                    wm.writerow(
                        [name, h["qseqid"], h["sseqid"], h["pident"], h["qcovs"], h["bitscore"]]
                    )
                # near misses only: best hit below threshold, for the unmatched list
                for h in dropped:
                    wd.writerow(
                        [
                            name,
                            h["qseqid"],
                            h["sseqid"],
                            h["pident"],
                            h["qcovs"],
                            h.get("reason", "below threshold"),
                        ]
                    )
                print(f"{name}: {len(mapped)} truth proteins mapped", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
