#!/usr/bin/env python3.12
"""One-off check: PRA3 orthologs and weak paralogs in each proteome, with tier and features.

Runs the mmseqs2 non-AVX2 binary (easy-search, PRA3 UniProt Q2TVJ9 as query) against each
proteome FASTA. Then joins the hits to the per-proteome tables written by cys_candidates.py.
Not part of cys_candidates.py. Python 3.12, standard library only.

An ortholog row is a hit with fident >= --min-ortholog-fident (default 0.90) and alignment
length >= --min-ortholog-alnlen (default 100). Other hits with e-value <= --max-evalue
(default 1e-3) are listed as weak paralog-like hits. These cut-offs are options, not facts.

STOP (stderr, exit 2) when the binary, a FASTA or a table is missing.
"""

import argparse
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cys_candidates as cc  # noqa: E402

MMSEQS_DEFAULT = "/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs"
FMT = "query,target,fident,alnlen,evalue,bits,qcov,tcov"


def run_mmseqs(mmseqs, query, target, tmp_dir, threads):
    out = Path(tmp_dir) / "hits.tsv"
    cmd = [
        mmseqs, "easy-search", str(query), str(target), str(out), str(Path(tmp_dir) / "work"),
        "--format-output", FMT, "-e", "1e-3", "-s", "7.5", "--threads", str(threads),
    ]  # fmt: skip
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise cc.Stop(f"mmseqs failed ({res.returncode}): {res.stderr.strip()[-300:]}")
    return [line.rstrip("\n").split("\t") for line in out.read_text().splitlines() if line]


def load_table(path):
    if not Path(path).is_file():
        raise cc.Stop(f"table not found: {path}")
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return {
            r[1]: dict(zip(head, r, strict=True))
            for r in (ln.rstrip("\n").split("\t") for ln in fh)
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="name, fasta, signalp[, domtbl] TSV")
    ap.add_argument("--table-dir", required=True, help="output dir of cys_candidates.py")
    ap.add_argument("--query", required=True, help="PRA3 UniProt FASTA (Q2TVJ9)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mmseqs", default=MMSEQS_DEFAULT)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-ortholog-fident", type=float, default=0.90)
    ap.add_argument("--min-ortholog-alnlen", type=int, default=100)
    ap.add_argument("--max-evalue", type=float, default=1e-3)
    a = ap.parse_args()
    try:
        if not Path(a.mmseqs).is_file():
            raise cc.Stop(f"mmseqs binary not found: {a.mmseqs}")
        if not Path(a.query).is_file():
            raise cc.Stop(f"query FASTA not found: {a.query}")
        jobs = cc.read_manifest(a.manifest)
        cols = ["proteome", "protein_id", "group", "fident", "alnlen", "evalue", "length",
                "sp_call", "mature_length", "cys_count", "max_cys_window", "tier"]  # fmt: skip
        lines = ["\t".join(cols)]
        for name, fasta, _sp, _dom in jobs:
            table = load_table(Path(a.table_dir) / f"{name}.tsv.gz")
            with tempfile.TemporaryDirectory() as tmp:
                hits = run_mmseqs(a.mmseqs, a.query, fasta, tmp, a.threads)
            best = {}
            for q, t, fid, aln, ev, *_ in hits:
                t = t.split()[0]
                if t not in best or float(ev) < float(best[t][3]):
                    best[t] = (q, float(fid), int(aln), ev)
            for t, (_q, fid, aln, ev) in sorted(best.items(), key=lambda kv: -kv[1][1]):
                is_orth = fid >= a.min_ortholog_fident and aln >= a.min_ortholog_alnlen
                if not is_orth and float(ev) > a.max_evalue:
                    continue
                if t not in table:
                    raise cc.Stop(f"{name}: mmseqs target {t!r} not in the table")
                r = table[t]
                lines.append(
                    "\t".join(
                        [
                            name,
                            t,
                            "ortholog" if is_orth else "weak_hit",
                            f"{fid:.3f}",
                            str(aln),
                            ev,
                            r["length"],
                            r["sp_call"],
                            r["mature_length"],
                            r["cys_count"],
                            r["max_cys_window"],
                            r["tier"],
                        ]
                    )  # fmt: skip
                )
        tmp_out = Path(a.out + ".tmp")
        tmp_out.write_text("\n".join(lines) + "\n")
        tmp_out.replace(a.out)
    except cc.Stop as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
