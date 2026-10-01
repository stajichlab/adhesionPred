"""Map truth-table gene IDs to protein FASTA records (spec 3.3 table). Standard library only.

Mapping types (column `id_mapping` of species.tsv):
- sgd:     GAF column 2 `S000...`  <-> FASTA header field `SGDID:S000...`
- cgd:     GAF column 11 first synonym (`C1_00010W_A`) <-> FASTA first token. The CGD FASTA
           headers carry no CAL... ID (checked 2026-10-01).
- pombase: GAF column 2 `SPAC1002.01` <-> FASTA `SPAC1002.01.1:pep` (transcript suffix removed)
- uniprot: GAF column 2 accession <-> FASTA `sp|ACC|NAME` or `tr|ACC|NAME`
Keys are compared after `normalize_id`: whitespace stripped, upper-cased, and a UniProt
isoform suffix (`-2`) removed for the uniprot mapping.
"""

import re
from collections.abc import Iterator
from pathlib import Path

from gaf import open_text
from seqhash import clean

_SGD = re.compile(r"SGDID:(S\d+)")
_POMBASE = re.compile(r"^(\S+)\.\d+:pep\b")
_UNIPROT = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
MAPPINGS = ("sgd", "cgd", "pombase", "uniprot")


class MappingError(ValueError):
    """Two FASTA records give the same key, or the mapping type is unknown."""


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    header, chunks = None, []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def normalize_id(value: str, mapping: str) -> str:
    key = value.strip().upper()
    if mapping == "uniprot":
        key = re.sub(r"-\d+$", "", key)
    return key


def fasta_key(header: str, mapping: str) -> str | None:
    if mapping == "sgd":
        match = _SGD.search(header)
        raw = match.group(1) if match else None
    elif mapping == "cgd":
        raw = header.split()[0] if header.split() else None
    elif mapping == "pombase":
        match = _POMBASE.match(header)
        raw = match.group(1) if match else None
    elif mapping == "uniprot":
        match = _UNIPROT.match(header)
        raw = match.group(1) if match else None
    else:
        raise MappingError(f"unknown id_mapping {mapping!r}")
    return normalize_id(raw, mapping) if raw else None


def gene_key(row: dict[str, str], mapping: str) -> str:
    value = row["synonym1"] if mapping == "cgd" else row["gene_id"]
    return normalize_id(value, mapping)


def index_fasta(path: str | Path, mapping: str) -> tuple[dict[str, tuple[str, str]], int]:
    """Return (key -> (fasta_id, sequence), number of identical duplicate records skipped).

    A key that occurs twice with the same cleaned sequence is kept once (the CGD Assembly 22
    file has 38 such records on 2026-10-01). A key that occurs twice with different sequences
    raises MappingError.
    """
    index: dict[str, tuple[str, str]] = {}
    duplicates = 0
    for header, seq in read_fasta(path):
        key = fasta_key(header, mapping)
        if key is None:
            continue
        if key in index:
            if clean(index[key][1]) != clean(seq):
                raise MappingError(f"{Path(path).name}: key {key} has two different sequences")
            duplicates += 1
            continue
        index[key] = (header.split()[0], seq)
    return index, duplicates


def attach(rows: list[dict[str, str]], index: dict[str, tuple[str, str]], mapping: str):
    """Return (matched, unmatched). matched rows gain fasta_id and sequence."""
    matched, unmatched = [], []
    for row in rows:
        hit = index.get(gene_key(row, mapping))
        if hit is None:
            unmatched.append(row)
        else:
            matched.append({**row, "fasta_id": hit[0], "sequence": hit[1]})
    return matched, unmatched
