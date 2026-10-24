"""Stage the input FASTA of the AGA1-class prediction step (local, deterministic, stdlib only).

Reads the 850 accessions of docs/reports/data/aga1_class_scoping/proteins.tsv, takes each sequence from the committed
data/controls/*/sequences.faa or the OUT-OF-REPO UniProt cache of the scoping study (same helper as the scoping scripts) and
writes, into a NEW directory outside the repo:
  aga1_proteins.faa (+ .sha256), lengths.tsv (accession, length, n_nonstandard, source), count.txt.
The committed manifest (manifest_accessions.txt, accessions only, sorted) is written with --write-manifest or checked with
the default mode. No sequence is ever written inside the repo. Downloaded/cached data is untrusted: run with `python3 -I`.

  python3 -I analysis/aga1_class/stage_inputs.py --out /path/NEW_DIR_OUTSIDE_REPO [--write-manifest]
"""

import argparse
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTEINS = ROOT / "docs/reports/data/aga1_class_scoping/proteins.tsv"
MANIFEST = HERE / "manifest_accessions.txt"
SCOPING_SCRIPTS = ROOT / "docs/reports/data/aga1_class_scoping/scripts"
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
N_EXPECTED = 850


def accessions_from_proteins(path=PROTEINS):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    col = rows[0].index("accession")
    acc = [r[col] for r in rows[1:]]
    if len(set(acc)) != len(acc):
        raise SystemExit("FATAL: duplicate accessions in proteins.tsv")
    if any(not a for a in acc):
        raise SystemExit("FATAL: empty accession in proteins.tsv")
    return sorted(acc)


def clean(seq):
    return seq.strip().upper().rstrip("*")


def build(accessions, seqs, src):
    """Return (fasta text, lengths rows). Fails on a missing or empty sequence."""
    missing = [a for a in accessions if not seqs.get(a)]
    if missing:
        raise SystemExit(f"FATAL: {len(missing)} accessions have no sequence: {missing[:10]}")
    lines, rows = [], []
    for a in accessions:
        s = clean(seqs[a])
        if not s:
            raise SystemExit(f"FATAL: empty sequence for {a}")
        lines.append(f">{a}")
        lines.extend(s[i : i + 60] for i in range(0, len(s), 60))
        rows.append((a, len(s), sum(c not in STANDARD for c in s), src.get(a, "")))
    return "\n".join(lines) + "\n", rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", help="NEW directory, outside the repo, that must not exist")
    ap.add_argument(
        "--write-manifest", action="store_true", help="(re)write the committed accession manifest"
    )
    args = ap.parse_args(argv)
    acc = accessions_from_proteins()
    if len(acc) != N_EXPECTED:
        raise SystemExit(f"FATAL: {len(acc)} accessions in proteins.tsv, expected {N_EXPECTED}")
    text = "\n".join(acc) + "\n"
    if args.write_manifest:
        MANIFEST.write_text(text)
    elif not MANIFEST.is_file() or MANIFEST.read_text() != text:
        raise SystemExit(
            "FATAL: manifest_accessions.txt differs from proteins.tsv (use --write-manifest)"
        )
    if not args.out:
        print(f"manifest ok: {len(acc)} accessions")
        return 0
    out = Path(args.out).resolve()
    if out == ROOT or ROOT in out.parents:
        raise SystemExit("FATAL: --out must be outside the repo (no sequence is committed)")
    if out.exists():
        raise SystemExit(f"FATAL: {out} exists; nothing is overwritten")
    sys.path.insert(0, str(SCOPING_SCRIPTS))
    from common import all_sequences  # noqa: E402, I001

    seqs, src = all_sequences()
    fasta, rows = build(acc, seqs, src)
    out.mkdir(parents=True)
    fa = out / "aga1_proteins.faa"
    fa.write_text(fasta)
    sha = hashlib.sha256(fasta.encode()).hexdigest()
    (out / "aga1_proteins.faa.sha256").write_text(sha + "\n")
    (out / "count.txt").write_text(f"{len(rows)}\n")
    with open(out / "lengths.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "length", "n_nonstandard", "source"])
        w.writerows(rows)
    over = sum(r[1] > 1022 for r in rows)
    print(f"staged {len(rows)} sequences in {out}; sha256 {sha}; {over} longer than 1022 aa; "
          f"{sum(r[2] > 0 for r in rows)} with non-standard residues")  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
