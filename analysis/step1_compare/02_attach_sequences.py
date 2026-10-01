#!/usr/bin/env python3
"""D1 step 2: attach protein sequences to truth-table genes and report unmatched IDs.

Reads $STEP1_WORKDIR/truth_set.tsv.gz and the FASTA files named in species.tsv. Writes
truth_sequences.tsv.gz, unmatched_ids.tsv and sequence_counts.tsv to the work directory.
Stops (exit 2, no truth_sequences.tsv.gz) if a source matches no gene, or if a P-ext or
ambiguous gene has no sequence, or if a sequence has a non-ASCII character. On a STOP an old
truth_sequences.tsv.gz is deleted; unmatched_ids.tsv, sequence_counts.tsv and sequence_run.json
are still written, so the cause can be read.

sequence_run.json records the selected sources, whether they are all sources (`all_sources`),
the SHA-256 of the truth set that was read and the FASTA hashes. A run with --sources is partial
and replaces the full outputs: check `all_sources` before you use the files. The script refuses a
truth_set.tsv.gz whose extract_log.json does not say `all_sources: true`, unless
--allow-partial-truth-set is given. Columns are listed in COLUMNS.md.
"""

import argparse
import json
import os
import platform
import subprocess
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
OUTPUT_NAMES = (
    "unmatched_ids.tsv",
    "sequence_counts.tsv",
    "sequence_run.json",
    "truth_sequences.tsv.gz",
)
UNMATCHED_COLUMNS = ("source_id", "gene_id", "symbol", "synonym1", "label", "reason")
SEQ_COUNT_COLUMNS = (
    "source_id",
    "truth_set_sha256",
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


def git_commit() -> str:
    """HEAD commit of the repository, or 'unknown'."""
    try:
        done = subprocess.run(
            ["git", "-C", str(paths.repo_root()), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return done.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def run(
    truth_rows,
    species_rows,
    manifest_rows,
    input_dir: Path,
    out_dir: Path,
    truth_set_sha256: str = "",
    provenance: dict | None = None,
):
    by_file = {r["file"]: r for r in manifest_rows}
    indexes: dict[str, tuple[dict, int, str]] = {}
    seq_rows, unmatched_rows, counts = [], [], []
    non_ascii: list[str] = []
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
            if not seq.isascii():
                non_ascii.append(f"{m['source_id']}: gene {m['gene_id']} has a non-ASCII sequence")
                continue
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
                "truth_set_sha256": truth_set_sha256,
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
    problems += non_ascii
    problems += [
        f"{u['source_id']}: {u['label']} gene {u['gene_id']} has no sequence"
        for u in unmatched_rows
        if u["label"] in (labels.P_EXT, labels.AMBIGUOUS)
    ]
    log = {
        "sources": [c["source_id"] for c in counts],
        "all_sources": None,
        "truth_set_sha256": truth_set_sha256,
        "fasta_sha256": {c["source_id"]: c["fasta_sha256"] for c in counts},
        "git_commit": git_commit(),
        "python": platform.python_version(),
        "arguments": [],
    }
    log.update(provenance or {})
    _write_outputs(out_dir, unmatched_rows, counts, log, None if problems else seq_rows)
    if problems:
        raise sequences.MappingError("; ".join(problems[:10]))
    return seq_rows, unmatched_rows, counts


def _write_outputs(out_dir: Path, unmatched_rows, counts, log, seq_rows) -> None:
    """Write to temp names, then os.replace. seq_rows is None on a STOP: the diagnostic files
    are replaced and an earlier truth_sequences.tsv.gz is deleted (the contract is: after a
    STOP there is no truth_sequences.tsv.gz)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    names = OUTPUT_NAMES if seq_rows is not None else OUTPUT_NAMES[:3]
    temps = {name: out_dir / f".tmp.{name}" for name in names}
    try:
        truth_table.write_tsv(temps["unmatched_ids.tsv"], UNMATCHED_COLUMNS, unmatched_rows)
        truth_table.write_tsv(temps["sequence_counts.tsv"], SEQ_COUNT_COLUMNS, counts)
        temps["sequence_run.json"].write_text(json.dumps(log, indent=2, sort_keys=True) + "\n")
        if seq_rows is not None:
            truth_table.write_tsv(temps["truth_sequences.tsv.gz"], SEQUENCE_COLUMNS, seq_rows)
        for name in names:
            os.replace(temps[name], out_dir / name)
        if seq_rows is None:
            (out_dir / "truth_sequences.tsv.gz").unlink(missing_ok=True)
    finally:
        for temp in temps.values():
            temp.unlink(missing_ok=True)


def _require_full_truth_set(log_path: Path) -> None:
    """Stop unless extract_log.json says that truth_set.tsv.gz covers all sources."""
    try:
        all_sources = json.loads(log_path.read_text()).get("all_sources")
    except (OSError, ValueError) as exc:
        raise manifest.DownloadError(
            f"cannot read {log_path} ({exc}); use --allow-partial-truth-set to skip this check"
        ) from exc
    if all_sources is not True:
        raise manifest.DownloadError(
            f"{log_path.name} does not say all_sources: true (value: {all_sources!r}), so "
            "truth_set.tsv.gz may hold only some sources; rerun 01_extract_go_truth.py without "
            "--sources or use --allow-partial-truth-set"
        )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None, help="source_id values to run")
    parser.add_argument(
        "--allow-partial-truth-set",
        action="store_true",
        help="use a truth_set.tsv.gz whose extract_log.json does not say all_sources: true",
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    input_dir = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    try:
        species_rows = truth_table.read_tsv(args.species)
        all_source_ids = [r["source_id"] for r in species_rows]
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
        truth_path = work / "truth_set.tsv.gz"
        if not args.allow_partial_truth_set:
            _require_full_truth_set(work / "extract_log.json")
        truth_rows = truth_table.read_tsv(truth_path)
        provenance = {
            "sources": [r["source_id"] for r in species_rows],
            "all_sources": [r["source_id"] for r in species_rows] == all_source_ids,
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        _, _, counts = run(
            truth_rows,
            species_rows,
            manifest.read_manifest(args.manifest),
            input_dir,
            work,
            truth_set_sha256=manifest.sha256_file(truth_path),
            provenance=provenance,
        )
    except (manifest.DownloadError, sequences.MappingError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in counts:
        print("\t".join(f"{k}={c[k]}" for k in ("source_id", "genes", "matched", "unmatched")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
