import sys
from pathlib import Path
from Bio import SeqIO  # Assuming Biopython is installed


def find_fasta_files(input_dir):
    """Finds all FASTA files in the input directory."""
    input_path = Path(input_dir)
    fasta_files = []
    # Search for all specified extensions, case-insensitively
    extensions = [".faa", ".pep", ".pep.fa", ".fasta"]

    for ext in extensions:
        fasta_files.extend(input_path.glob(f"**/*{ext}"))
        fasta_files.extend(input_path.glob(f"**/*{ext.upper()}"))

    return fasta_files


def process_fasta_file(file_path):
    """Reads and processes sequences from a single FASTA file."""
    print(f"Processing file: {file_path}")
    sequences = []
    try:
        # Use Biopython's SeqIO for robust FASTA parsing
        for seq_record in SeqIO.parse(file_path, "fasta"):
            sequences.append({"id": seq_record.id, "sequence": str(seq_record.seq)})
        print(f"  Found {len(sequences)} sequences.")
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)

    return sequences


def main(input_dir="input"):
    """Main function to run the adhesion prediction framework."""
    print("Starting adhesion prediction framework...")

    # 1. Find input files
    fasta_files = find_fasta_files(input_dir)

    if not fasta_files:
        print(
            f"No FASTA files found in '{input_dir}' with extensions: "
            ".faa, .pep, .pep.fa, .fasta"
        )
        return

    print(f"Found {len(fasta_files)} FASTA files to process.")

    all_sequences = []
    # 2. Process each file
    for file_path in fasta_files:
        sequences = process_fasta_file(file_path)
        all_sequences.extend(sequences)

    # 3. Placeholder for prediction logic (future development)
    print(f"\nTotal sequences loaded: {len(all_sequences)}")
    print("--- Adhesion Prediction Logic Goes Here ---")


if __name__ == "__main__":
    # Allow input directory to be specified as a command-line argument
    if len(sys.argv) > 1:
        input_dir_path = sys.argv[1]
    else:
        input_dir_path = "input"

    main(input_dir_path)
