#!/usr/bin/env python3
"""
Generate a random subset of sequence IDs from *_AnnotatedProteins.fasta files.

Steps:
1) Read all *.FASTA.tab files and build an IgnoredIDs dictionary from a column.
2) Parse FASTA IDs, remove any in IgnoredIDs, then randomly sample a fixed number.
"""

import argparse
import glob
import os
import random
import sys
from typing import Dict, Iterable, List


def read_ignored_ids(tab_files: Iterable[str], column_index: int) -> Dict[str, bool]:
    ignored_ids: Dict[str, bool] = {}
    for tab_file in tab_files:
        try:
            with open(tab_file, encoding="utf-8") as handle:
                for line_num, line in enumerate(handle, 1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) < column_index:
                        print(
                            f"Warning: {tab_file}:{line_num} has fewer than {column_index} columns, skipping",
                            file=sys.stderr,
                        )
                        continue
                    ignored_ids[fields[column_index - 1]] = True
        except FileNotFoundError:
            print(f"Warning: File not found: {tab_file}", file=sys.stderr)
        except OSError as exc:
            print(f"Warning: Could not read {tab_file}: {exc}", file=sys.stderr)
    return ignored_ids


def read_fasta_ids(fasta_file: str) -> List[str]:
    ids: List[str] = []
    seen = set()
    with open(fasta_file, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            seq_id = line[1:].strip().split()[0]
            if not seq_id or seq_id in seen:
                continue
            seen.add(seq_id)
            ids.append(seq_id)
    return ids


def find_files(input_dir: str, pattern: str) -> List[str]:
    search_pattern = os.path.join(input_dir, "**", pattern)
    return sorted(glob.glob(search_pattern, recursive=True))


def write_sampled_ids(output_dir: str, fasta_file: str, sampled_ids: List[str]) -> str:
    base = os.path.basename(fasta_file)
    if base.endswith(".fasta"):
        base = base[: -len(".fasta")]
    output_path = os.path.join(output_dir, f"{base}.sampled_ids.txt")
    with open(output_path, "w", encoding="utf-8") as handle:
        for seq_id in sampled_ids:
            handle.write(f"{seq_id}\n")
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Randomly sample IDs from *_AnnotatedProteins.fasta files after removing IDs found in *.FASTA.tab files."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Root directory to search for input files (default: current directory)",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to write sampled ID files (default: current directory)",
    )
    parser.add_argument(
        "--tab-pattern",
        default="*.FASTA.tab",
        help="Glob pattern for tab files (default: *.FASTA.tab)",
    )
    parser.add_argument(
        "--fasta-pattern",
        default="*_AnnotatedProteins.fasta",
        help="Glob pattern for FASTA files (default: *_AnnotatedProteins.fasta)",
    )
    parser.add_argument(
        "--tab-column",
        type=int,
        default=2,
        help="1-based column index for IDs in *.FASTA.tab files (default: 2)",
    )
    parser.add_argument(
        "--num-ids",
        type=int,
        default=100,
        help="Number of IDs to sample per FASTA file (default: 100)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible sampling (default: random)",
    )

    args = parser.parse_args()

    if args.tab_column < 1:
        print("Error: --tab-column must be 1 or greater", file=sys.stderr)
        return 2
    if args.num_ids < 1:
        print("Error: --num-ids must be 1 or greater", file=sys.stderr)
        return 2

    tab_files = find_files(args.input_dir, args.tab_pattern)
    fasta_files = find_files(args.input_dir, args.fasta_pattern)

    if not fasta_files:
        print("Error: No FASTA files found", file=sys.stderr)
        return 1

    ignored_ids = read_ignored_ids(tab_files, args.tab_column)
    rng = random.Random(args.seed)

    os.makedirs(args.output_dir, exist_ok=True)

    for fasta_file in fasta_files:
        fasta_ids = read_fasta_ids(fasta_file)
        filtered_ids = [seq_id for seq_id in fasta_ids if seq_id not in ignored_ids]
        rng.shuffle(filtered_ids)
        if len(filtered_ids) < args.num_ids:
            print(
                f"Warning: {fasta_file} has {len(filtered_ids)} IDs after filtering; sampling all available",
                file=sys.stderr,
            )
        sampled_ids = filtered_ids[: args.num_ids]
        output_path = write_sampled_ids(args.output_dir, fasta_file, sampled_ids)
        print(f"Wrote {len(sampled_ids)} IDs to {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
