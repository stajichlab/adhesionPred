"""Input/output functions for handling FASTA files."""

import gzip
import sys
from pathlib import Path

from Bio import SeqIO


def clean_sequence(sequence):
    """Clean protein sequence by replacing invalid amino acids.

    Args:
        sequence: Protein sequence string.

    Returns:
        Cleaned sequence with 'J' replaced by 'L'.
    """
    return sequence.replace("J", "L").replace("j", "L")


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
