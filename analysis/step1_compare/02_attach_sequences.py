#!/usr/bin/env python3
"""D1 step 2: attach protein sequences to truth-table genes and report unmatched IDs.

Reads $STEP1_WORKDIR/truth_set.tsv.gz and the FASTA files named in species.tsv. Writes
truth_sequences.tsv.gz, unmatched_ids.tsv and sequence_counts.tsv to the work directory.
Stops (exit 2, no truth_sequences.tsv.gz) if a source matches no gene, or if a P-ext or
ambiguous gene has no sequence; unmatched_ids.tsv and sequence_counts.tsv are still written.
"""

import argparse
import os
import sys
from pathlib import Path

import labels
import manifest
import paths
import seqhash
import sequences
import truth_table

SEQUENCE_COLUMNS = (
    "source_id",
    "gene_id",
    "label",
    "fasta_id",
    "length",
    "seq_sha256",
    "sequence",
)
OUTPUT_NAMES = ("unmatched_ids.tsv", "sequence_counts.tsv", "truth_sequences.tsv.gz")
UNMATCHED_COLUMNS = ("source_id", "gene_id", "symbol", "synonym1", "label")
SEQ_COUNT_COLUMNS = (
    "source_id",
    "fasta_file",
    "fasta_sha256",
    "fasta_duplicate_records",
    "genes",
    "matched",
    "unmatched",
    "unmatched_p_ext",
    "unmatched_ambiguous",
    "unmatched_n_int",
    "unmatched_n_sec",
)


def run(truth_rows, species_rows, manifest_rows, input_dir: Path, out_dir: Path):
    by_file = {r["file"]: r for r in manifest_rows}
    indexes: dict[str, tuple[dict, int, str]] = {}
    seq_rows, unmatched_rows, counts = [], [], []
    for sp in species_rows:
        rows = [r for r in truth_rows if r["source_id"] == sp["source_id"]]
        if not rows:
            continue
        fasta = sp["fasta_file"]
        if fasta not in indexes:
            path = input_dir / fasta
            if fasta not in by_file:
                raise manifest.DownloadError(f"{fasta}: no manifest row")
            sha = manifest.verify_against_manifest(path, by_file[fasta])
            index, duplicates = sequences.index_fasta(path, sp["id_mapping"])
            indexes[fasta] = (index, duplicates, sha)
        index, duplicates, sha = indexes[fasta]
        matched, unmatched = sequences.attach(rows, index, sp["id_mapping"])
        for m in matched:
            seq = seqhash.clean(m["sequence"])
            seq_rows.append(
                {
                    "source_id": m["source_id"],
                    "gene_id": m["gene_id"],
                    "label": m["label"],
                    "fasta_id": m["fasta_id"],
                    "length": str(len(seq)),
                    "seq_sha256": seqhash.seq_sha256(seq),
                    "sequence": seq,
                }
            )
        unmatched_rows.extend({c: u[c] for c in UNMATCHED_COLUMNS} for u in unmatched)
        counts.append(
            {
                "source_id": sp["source_id"],
                "fasta_file": fasta,
                "fasta_sha256": sha,
                "fasta_duplicate_records": str(duplicates),
                "genes": str(len(rows)),
                "matched": str(len(matched)),
                "unmatched": str(len(unmatched)),
                "unmatched_p_ext": str(sum(u["label"] == labels.P_EXT for u in unmatched)),
                "unmatched_ambiguous": str(sum(u["label"] == labels.AMBIGUOUS for u in unmatched)),
                "unmatched_n_int": str(sum(u["label"] == labels.N_INT for u in unmatched)),
                "unmatched_n_sec": str(sum(u["label"] == labels.N_SEC for u in unmatched)),
            }
        )
    problems = [f"{c['source_id']}: no gene matched" for c in counts if c["matched"] == "0"]
    problems += [
        f"{u['source_id']}: {u['label']} gene {u['gene_id']} has no sequence"
        for u in unmatched_rows
        if u["label"] in (labels.P_EXT, labels.AMBIGUOUS)
    ]
    _write_outputs(out_dir, unmatched_rows, counts, None if problems else seq_rows)
    if problems:
        raise sequences.MappingError("; ".join(problems[:10]))
    return seq_rows, unmatched_rows, counts


def _write_outputs(out_dir: Path, unmatched_rows, counts, seq_rows) -> None:
    """Write to temp names, then os.replace. seq_rows is None on a STOP: only the two
    diagnostic files are replaced, and an earlier truth_sequences.tsv.gz stays untouched."""
    out_dir.mkdir(parents=True, exist_ok=True)
    names = OUTPUT_NAMES if seq_rows is not None else OUTPUT_NAMES[:2]
    temps = {name: out_dir / f".tmp.{name}" for name in names}
    try:
        truth_table.write_tsv(temps["unmatched_ids.tsv"], UNMATCHED_COLUMNS, unmatched_rows)
        truth_table.write_tsv(temps["sequence_counts.tsv"], SEQ_COUNT_COLUMNS, counts)
        if seq_rows is not None:
            truth_table.write_tsv(temps["truth_sequences.tsv.gz"], SEQUENCE_COLUMNS, seq_rows)
        for name in names:
            os.replace(temps[name], out_dir / name)
    finally:
        for temp in temps.values():
            temp.unlink(missing_ok=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None, help="source_id values to run")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    input_dir = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    try:
        species_rows = truth_table.read_tsv(args.species)
        if args.sources is not None:
            valid = [r["source_id"] for r in species_rows]
            for source_id in args.sources:
                if source_id not in valid:
                    raise manifest.DownloadError(
                        f"unknown source {source_id}; valid ids: {', '.join(valid)}"
                    )
            species_rows = [r for r in species_rows if r["source_id"] in args.sources]
        if not species_rows:
            raise manifest.DownloadError("no sources selected")
        truth_rows = truth_table.read_tsv(work / "truth_set.tsv.gz")
        _, _, counts = run(
            truth_rows, species_rows, manifest.read_manifest(args.manifest), input_dir, work
        )
    except (manifest.DownloadError, sequences.MappingError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in counts:
        print("\t".join(f"{k}={c[k]}" for k in ("source_id", "genes", "matched", "unmatched")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
