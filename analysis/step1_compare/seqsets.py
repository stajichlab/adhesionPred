"""Phase B sequence universe: members of every input set and the unique sequences.

Standard library only. A member is one protein of one input set (one row per set_id,
source_id, gene_id). Members with the same cleaned sequence share one seq_sha256 and are
embedded and scored once. Unique sequences are sorted by seq_sha256; `row` is the position in
that order and is the row of the N-terminal embedding matrix. `cterm_row` numbers the
sequences longer than MAX_RESIDUES in the same order and is the row of the C-terminal matrix.
"""

from collections.abc import Iterable
from pathlib import Path

import seqhash
import sequences
import seqwindow
import truth_table

SET_COLUMNS = ("set_id", "kind", "location", "note")
KINDS = ("truth", "keyword", "download", "site")
MEMBER_COLUMNS = ("set_id", "source_id", "gene_id", "seq_sha256", "length")
UNIQUE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", "sequence")


class SequenceSetError(ValueError):
    """An input set is malformed, missing, or holds a sequence ESM-2 cannot take as it is."""


def read_sets(path: str | Path) -> list[dict[str, str]]:
    rows = truth_table.read_tsv(path)
    seen = set()
    for row in rows:
        missing = [c for c in SET_COLUMNS if c not in row]
        if missing:
            raise SequenceSetError(f"{path}: row {row} lacks {missing}")
        if row["kind"] not in KINDS:
            raise SequenceSetError(f"{path}: set {row['set_id']} has unknown kind {row['kind']}")
        if row["set_id"] in seen:
            raise SequenceSetError(f"{path}: set_id {row['set_id']} occurs twice")
        seen.add(row["set_id"])
    return rows


def resolve(row: dict[str, str], work: Path, downloads: Path, site_value) -> Path:
    """Path of a set's input. `site` locations are `<site.yaml key>:<relative path>`."""
    kind, location = row["kind"], row["location"]
    if kind in ("truth", "keyword"):
        return Path(work) / location
    if kind == "download":
        return Path(downloads) / location
    key, _, rel = location.partition(":")
    if not rel:
        raise SequenceSetError(f"set {row['set_id']}: site location needs 'key:path'")
    return Path(site_value(key)) / rel


def _check(seq: str, where: str) -> None:
    bad = seqwindow.bad_characters(seq)
    if bad:
        raise SequenceSetError(f"{where}: characters {sorted(bad)} are not ESM-2 residues")


def truth_members(rows: Iterable[dict[str, str]], set_id: str = "truth"):
    """Members from truth_sequences.tsv.gz rows. The stored seq_sha256 is checked."""
    members, seqs = [], {}
    for r in rows:
        seq = r["sequence"]
        digest = seqhash.seq_sha256(seq)
        if digest != r["seq_sha256"] or seq != seqhash.clean(seq):
            raise SequenceSetError(
                f"truth_sequences: {r['source_id']}:{r['gene_id']} has a stored seq_sha256 or "
                "sequence that does not match the cleaned sequence; re-run 02_attach_sequences.py"
            )
        _check(seq, f"{set_id} {r['source_id']}:{r['gene_id']}")
        members.append(_member(set_id, r["source_id"], r["gene_id"], digest, seq))
        seqs[digest] = seq
    return members, seqs


def fasta_members(set_id: str, path: str | Path, kind: str):
    """Members from a FASTA file. Return (members, sequences by hash, empty records).

    gene_id is the UniProt accession for kind `keyword`, else the first header token.
    A gene_id that occurs twice with the same cleaned sequence is kept once; with a different
    sequence it raises SequenceSetError."""
    members, seqs, by_gene = [], {}, {}
    empty = 0
    for header, raw in sequences.read_fasta(path):
        seq = seqhash.clean(raw)
        if not seq:
            empty += 1
            continue
        if kind == "keyword":
            gene_id = sequences.fasta_key(header, "uniprot")
            if gene_id is None:
                raise SequenceSetError(f"{set_id}: header {header[:40]!r} is not UniProt FASTA")
        else:
            gene_id = header.split()[0]
        digest = seqhash.seq_sha256(seq)
        if gene_id in by_gene:
            if by_gene[gene_id] != digest:
                if kind == "keyword":
                    raise SequenceSetError(
                        f"{set_id}: two records give accession {gene_id} after the isoform "
                        "suffix ('-<n>') is removed, and their sequences differ"
                    )
                raise SequenceSetError(f"{set_id}: {gene_id} has two different sequences")
            continue
        _check(seq, f"{set_id} {gene_id}")
        by_gene[gene_id] = digest
        members.append(_member(set_id, set_id, gene_id, digest, seq))
        seqs[digest] = seq
    return members, seqs, empty


def _member(set_id, source_id, gene_id, digest, seq) -> dict[str, str]:
    return {
        "set_id": set_id,
        "source_id": source_id,
        "gene_id": gene_id,
        "seq_sha256": digest,
        "length": str(len(seq)),
    }


def unique_rows(seqs: dict[str, str]) -> list[dict[str, str]]:
    """One row per unique sequence, sorted by seq_sha256, with row and cterm_row."""
    rows, cterm = [], 0
    for i, digest in enumerate(sorted(seqs)):
        seq = seqs[digest]
        crow = ""
        if seqwindow.needs_cterm(seq):
            crow = str(cterm)
            cterm += 1
        rows.append(
            {
                "row": str(i),
                "seq_sha256": digest,
                "length": str(len(seq)),
                "cterm_row": crow,
                "sequence": seq,
            }
        )
    return rows
