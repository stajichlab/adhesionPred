"""Read and check a protein FASTA (plain, .gz or .zst).

Rules (spec section 3.8): a trailing ``*`` is stripped; an internal ``*`` or a character that is not
a residue makes the protein ``na_invalid``; ``X B Z U J O`` are allowed and counted; duplicate IDs,
an empty file and a file with no valid protein stop the run.
"""

import gzip
import hashlib
import subprocess
import zlib
from dataclasses import dataclass
from pathlib import Path

OK = "ok"
NA_INVALID = "na_invalid"
STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS = frozenset("XBZUJO")


class FastaError(ValueError):
    """The FASTA cannot be used; the run must stop."""


@dataclass(frozen=True)
class Protein:
    id: str
    sequence: str
    sha256: str
    state: str
    note: str
    ambiguous_fraction: float
    trailing_stop: bool = False


def _lines(path):
    path = Path(path)
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8-sig") as fh:
                yield from fh
        elif path.suffix == ".zst":
            out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True)
            yield from out.stdout.decode("utf-8-sig").splitlines()
        else:
            with open(path, encoding="utf-8-sig") as fh:
                yield from fh
    except (UnicodeDecodeError, EOFError, zlib.error) as err:
        raise FastaError(f"{path}: cannot be read as text ({err.__class__.__name__})") from err


def _records(path):
    header, chunks = None, []
    for raw in _lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].strip(), []
        elif header is not None:
            chunks.append(line)
        else:
            raise FastaError("sequence lines before the first header")
    if header is not None:
        yield header, "".join(chunks)


def _check(sequence):
    if not sequence:
        return NA_INVALID, "empty sequence", 0.0
    if "*" in sequence:
        return NA_INVALID, "internal stop codon", 0.0
    bad = sorted(set(sequence) - STANDARD - AMBIGUOUS)
    if bad:
        return NA_INVALID, "not residues: " + "".join(bad), 0.0
    ambiguous = sum(1 for c in sequence if c in AMBIGUOUS) / len(sequence)
    return OK, "", ambiguous


def read_fasta(path):
    proteins, seen = [], set()
    for header, raw_seq in _records(path):
        if not header:
            raise FastaError("a record has an empty header")
        pid = header.split()[0]
        if pid in seen:
            raise FastaError(f"duplicate ID: {pid}")
        seen.add(pid)
        sequence = raw_seq.upper()
        trailing = sequence.endswith("*")
        if trailing:
            sequence = sequence[:-1]
        state, note, ambiguous = _check(sequence)
        sha = hashlib.sha256(sequence.encode()).hexdigest()
        proteins.append(Protein(pid, sequence, sha, state, note, ambiguous, trailing))
    if not proteins:
        raise FastaError("the FASTA has no records")
    if all(p.state != OK for p in proteins):
        raise FastaError("the FASTA has no valid proteins")
    return proteins
