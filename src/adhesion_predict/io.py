"""Input/output functions for handling FASTA files."""

import gzip
import multiprocessing as mp
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO


def clean_sequence(sequence):
    """Clean protein sequence by replacing invalid amino acids and removing stop codons.

    Args:
        sequence: Protein sequence string.

    Returns:
        Cleaned sequence with 'J' replaced by 'L' and '*' characters removed.
    """
    return sequence.replace("J", "L").replace("j", "L").replace("*", "")


def find_fasta_files(input_dir):
    """Finds all FASTA files in the input directory.

    Args:
        input_dir: Path to the input directory.

    Returns:
        List of Path objects for FASTA files found.
    """
    input_path = Path(input_dir)
    fasta_files = []
    extensions = [".aa", ".faa", ".pep", ".fa", ".fasta"]
    gzip_extensions = [".aa.gz", ".faa.gz", ".pep.gz", ".fa.gz", ".fasta.gz"]

    for ext in extensions + gzip_extensions:
        fasta_files.extend(input_path.glob(f"**/*{ext}"))
        fasta_files.extend(input_path.glob(f"**/*{ext.upper()}"))

    return fasta_files


def process_fasta_files_parallel(fasta_files, max_workers=None):
    """Process multiple FASTA files in parallel.

    Args:
        fasta_files: List of Path objects for FASTA files.
        max_workers: Maximum number of worker processes. If None, uses CPU count.

    Returns:
        List of all sequences from all files.
    """
    if not fasta_files:
        return []

    # Single file doesn't benefit from parallel processing
    if len(fasta_files) == 1:
        return process_fasta_file(fasta_files[0])

    if max_workers is None:
        max_workers = min(len(fasta_files), mp.cpu_count())

    print(f"Processing {len(fasta_files)} files using {max_workers} workers...")

    all_sequences = []
    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all file processing tasks
            future_to_file = {
                executor.submit(process_fasta_file, fasta_file): fasta_file
                for fasta_file in fasta_files
            }

            # Collect results as they complete
            for future in as_completed(future_to_file):
                fasta_file = future_to_file[future]
                try:
                    sequences = future.result()
                    all_sequences.extend(sequences)
                except Exception as e:
                    print(f"Error processing {fasta_file}: {e}", file=sys.stderr)

    except Exception as e:
        print(f"Error in parallel processing, falling back to sequential: {e}", file=sys.stderr)
        # Fallback to sequential processing
        for fasta_file in fasta_files:
            sequences = process_fasta_file(fasta_file)
            all_sequences.extend(sequences)

    print(f"Parallel processing completed. Total sequences: {len(all_sequences)}")
    return all_sequences


def process_fasta_file(file_path):
    """Reads and processes sequences from a single FASTA file.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        List of dictionaries with 'id' and 'sequence' keys.
    """
    print(f"Processing file: {file_path}")
    sequences = []
    try:
        if file_path.suffix == ".gz":
            print(f"Reading gzipped file: {file_path}")
            handle = gzip.open(file_path, "rt", encoding="utf-8")
        else:
            handle = open(file_path, encoding="utf-8")

        for seq_record in SeqIO.parse(handle, "fasta"):
            sequence = str(seq_record.seq)
            sequence = clean_sequence(sequence)
            sequences.append({"id": seq_record.id, "sequence": sequence})
        print(f"  Found {len(sequences)} sequences.")
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)

    return sequences


def load_sequences_from_dir(input_dir):
    """Load all sequences from FASTA files in a directory.

    Args:
        input_dir: Path to the input directory.

    Returns:
        List of dictionaries with 'id', 'sequence', and 'label' (optional) keys.
    """
    input_path = Path(input_dir)
    sequences = []

    extensions = [".aa", ".faa", ".pep", ".fa", ".fasta"]
    gzip_extensions = [".aa.gz", ".faa.gz", ".pep.gz", ".fa.gz", ".fasta.gz"]

    for ext in extensions + gzip_extensions:
        for fasta_file in input_path.glob(f"**/*{ext}"):
            fasta_file = fasta_file.resolve()
            if fasta_file.suffix == ".gz":
                print(f"Reading gzipped file: {fasta_file}")
                handle = gzip.open(fasta_file, "rt", encoding="utf-8")
            else:
                handle = open(fasta_file, encoding="utf-8")
            for seq_record in SeqIO.parse(handle, "fasta"):
                sequence = str(seq_record.seq)
                sequence = clean_sequence(sequence)
                sequences.append(
                    {
                        "id": seq_record.id,
                        "sequence": sequence,
                    }
                )

    return sequences
