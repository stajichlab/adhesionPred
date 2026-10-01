#!/usr/bin/env python3
"""Phase B step 1: collect every input set and dedupe the sequences by seq_sha256.

Reads the sets in sequence_sets.tsv: truth_sequences.tsv.gz (02), keyword_sequences.fasta.gz
(04), the proteome FASTA files in $STEP1_WORKDIR/downloads, and the C. immitis RS FASTA (path
from config/site.yaml). Writes to $STEP1_WORKDIR/phaseb/:

  sequence_members.tsv.gz   one row per (set_id, source_id, gene_id)
  unique_sequences.tsv.gz   one row per unique sequence (row, seq_sha256, length, cterm_row,
                            sequence), sorted by seq_sha256
  unique_sequences.fasta.gz the same sequences as FASTA; the header is the seq_sha256
  prepare_run.json          input hashes, counts, all_sources, git commit, arguments

STOP (exit 2, no output written): a missing input; a run log (sequence_run.json,
keyword_tier_run.json) that does not say all_sources: true, unless --allow-partial-truth-set;
a stored seq_sha256 that does not match its sequence; a character that ESM-2 does not accept;
one gene_id with two different sequences in one set. Columns are listed in COLUMNS.md.
"""

import argparse
import gzip
import sys
from pathlib import Path

import manifest
import paths
import runinfo
import seqsets
import truth_table

OUTPUT_NAMES = (
    "sequence_members.tsv.gz",
    "unique_sequences.tsv.gz",
    "unique_sequences.fasta.gz",
    "prepare_run.json",
)


def phaseb_dir(work: Path) -> Path:
    return Path(work) / "phaseb"


def write_fasta_gz(path: Path, rows) -> None:
    """FASTA with the seq_sha256 as header, gzip with a fixed mtime (byte-stable)."""
    data = "".join(f">{r['seq_sha256']}\n{r['sequence']}\n" for r in rows).encode("ascii")
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(data)


def run(set_rows, work: Path, downloads: Path, site_value, provenance: dict | None = None):
    members, seqs, inputs = [], {}, {}
    for row in set_rows:
        path = seqsets.resolve(row, work, downloads, site_value)
        if not path.exists():
            raise seqsets.SequenceSetError(f"set {row['set_id']}: {path} not found")
        if row["kind"] == "truth":
            got, got_seqs = seqsets.truth_members(truth_table.read_tsv(path), row["set_id"])
            empty = 0
        else:
            got, got_seqs, empty = seqsets.fasta_members(row["set_id"], path, row["kind"])
        if not got:
            raise seqsets.SequenceSetError(f"set {row['set_id']}: {path} gave no sequence")
        members += got
        seqs.update(got_seqs)
        inputs[row["set_id"]] = {
            "file": str(path),
            "sha256": manifest.sha256_file(path),
            "members": len(got),
            "empty_records": empty,
        }
    unique = seqsets.unique_rows(seqs)
    lengths = [int(r["length"]) for r in unique]
    log = {
        "inputs": inputs,
        "members": len(members),
        "unique_sequences": len(unique),
        "unique_residues": sum(lengths),
        "over_max_residues": sum(r["cterm_row"] != "" for r in unique),
        "all_sources": None,
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
        "arguments": [],
    }
    log.update(provenance or {})
    out = phaseb_dir(work)
    runinfo.atomic_write_all(
        out,
        {
            "sequence_members.tsv.gz": lambda p: truth_table.write_tsv(
                p, seqsets.MEMBER_COLUMNS, members
            ),
            "unique_sequences.tsv.gz": lambda p: truth_table.write_tsv(
                p, seqsets.UNIQUE_COLUMNS, unique
            ),
            "unique_sequences.fasta.gz": lambda p: write_fasta_gz(p, unique),
            "prepare_run.json": lambda p: runinfo.write_json(p, log),
        },
    )
    return members, unique, log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sets", default=str(paths.STEP1_DIR / "sequence_sets.tsv"))
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument(
        "--allow-partial-truth-set",
        action="store_true",
        help="use inputs whose run logs do not say all_sources: true",
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    downloads = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    try:
        set_rows = seqsets.read_sets(args.sets)
        kinds = {r["kind"] for r in set_rows}
        allow = args.allow_partial_truth_set
        logs = []
        if "truth" in kinds:
            logs.append("sequence_run.json")
        if "keyword" in kinds:
            logs.append("keyword_tier_run.json")
        for name in logs:
            runinfo.require_full(work / name, name, allow)
        provenance = {
            "all_sources": all(runinfo.says_all_sources(work / n) for n in logs),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        _, unique, log = run(set_rows, work, downloads, paths.site_value, provenance)
    except (seqsets.SequenceSetError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(
        f"members={log['members']} unique={log['unique_sequences']} "
        f"residues={log['unique_residues']} over_1022={log['over_max_residues']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
