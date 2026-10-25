"""Stage UniProt reference-proteome FASTA files for the AGA1-class proteome test (local, deterministic, stdlib only).

Reads one FASTA per species (UniProt `sp|ACC|ENTRY_SPECIES description` or `tr|` headers), writes into a NEW directory outside the
repo: chunks/<species>_NNN.faa (bare-accession headers, at most --chunk-size proteins each, input order), chunks.tsv (chunk,
species, n_proteins, sha256), proteins.tsv (accession, species, chunk, length, n_nonstandard, entry_name, gene), manifest.txt
(all accessions, sorted) and SHA256SUMS. No sequence is ever written inside the repo. Downloaded data is untrusted: run with
`python3 -I`. A trailing `*` is removed; nothing else in a sequence is changed (X, U, B, Z are kept and counted).

  python3 -I analysis/aga1_class/stage_proteomes.py --out /path/NEW_DIR_OUTSIDE_REPO [--chunk-size 1000] \
      Scer_S288C=/path/Scer.faa Calb_SC5314=/path/Calb.faa ...

Fails loudly (FATAL:, exit 1) on an existing --out, a repeated species name, an empty or header-less file, a duplicate
accession (within or across files), an empty sequence or a record without an accession.
"""

import argparse
import hashlib
import re
import sys
from pathlib import Path

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


class StageError(Exception):
    pass


def parse_fasta(path):
    """Yield (accession, entry_name, gene, sequence) in file order."""
    acc, entry, gene = None, "", ""
    seq = []
    out = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.startswith(">"):
            if acc is not None:
                out.append((acc, entry, gene, "".join(seq)))
            m = re.match(r">(?:sp|tr)\|([A-Za-z0-9_-]+)\|(\S+)", line)
            if not m:
                raise StageError(f"{path}:{n}: header without a UniProt accession: {line[:60]!r}")
            acc, entry = m.group(1), m.group(2)
            g = re.search(r"\sGN=(\S+)", line)
            gene = g.group(1) if g else ""
            seq = []
        elif line.strip():
            if acc is None:
                raise StageError(f"{path}:{n}: sequence before any header")
            seq.append(line.strip().upper().rstrip("*"))
    if acc is not None:
        out.append((acc, entry, gene, "".join(seq)))
    if not out:
        raise StageError(f"{path}: no records")
    empty = [a for a, _, _, s in out if not s]
    if empty:
        raise StageError(f"{path}: {len(empty)} records without a sequence, first {empty[0]}")
    return out


def stage(species_files, out_dir, chunk_size=1000):
    if chunk_size < 1:
        raise StageError("--chunk-size must be at least 1")
    names = [n for n, _ in species_files]
    if len(set(names)) != len(names):
        raise StageError("a species name is given twice")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise StageError(f"{out_dir} exists; nothing is overwritten")
    seen, chunk_rows, protein_rows = {}, [], []
    files = {}
    for species, path in species_files:
        recs = parse_fasta(path)
        for i in range(0, len(recs), chunk_size):
            chunk = f"{species}_{i // chunk_size + 1:03d}"
            part = recs[i : i + chunk_size]
            text = "".join(f">{a}\n{s}\n" for a, _, _, s in part)
            files[chunk] = text
            chunk_rows.append(
                (chunk, species, len(part), hashlib.sha256(text.encode()).hexdigest())
            )
            for a, entry, gene, s in part:
                if a in seen:
                    raise StageError(f"accession {a} occurs in {seen[a]} and {species}")
                seen[a] = species
                nonstd = sum(c not in STANDARD for c in s)
                protein_rows.append((a, species, chunk, len(s), nonstd, entry, gene))
    (out_dir / "chunks").mkdir(parents=True)
    for chunk, text in files.items():
        (out_dir / "chunks" / f"{chunk}.faa").write_text(text)
    (out_dir / "chunks.tsv").write_text(
        "chunk\tspecies\tn_proteins\tsha256\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in chunk_rows)
    )
    (out_dir / "proteins.tsv").write_text(
        "accession\tspecies\tchunk\tlength\tn_nonstandard\tentry_name\tgene\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in protein_rows)
    )
    (out_dir / "manifest.txt").write_text("".join(f"{a}\n" for a in sorted(seen)))
    sums = [f"{r[3]}  chunks/{r[0]}.faa" for r in chunk_rows]
    for name in ("chunks.tsv", "proteins.tsv", "manifest.txt"):
        sums.append(f"{hashlib.sha256((out_dir / name).read_bytes()).hexdigest()}  {name}")
    (out_dir / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    return chunk_rows, protein_rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk-size", type=int, default=1000)
    ap.add_argument("species", nargs="+", help="NAME=FASTA")
    a = ap.parse_args(argv)
    try:
        pairs = []
        for item in a.species:
            if "=" not in item:
                raise StageError(f"expected NAME=FASTA, got {item!r}")
            n, p = item.split("=", 1)
            pairs.append((n, p))
        chunks, prots = stage(pairs, a.out, a.chunk_size)
    except (StageError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    print(f"wrote {a.out}: {len(prots)} proteins in {len(chunks)} chunks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
