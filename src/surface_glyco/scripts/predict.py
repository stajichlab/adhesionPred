#!/usr/bin/env python
"""Prediction script for scoring fungal cell-surface glycoproteins."""

import argparse
import csv
import sys
from pathlib import Path

from surface_glyco.config import DEFAULT_MODEL, get_models_dir, model_filename
from surface_glyco.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from surface_glyco.io import find_fasta_files, process_fasta_file, process_fasta_files_parallel
from surface_glyco.model import load_model, load_model_card, predict, predict_proba

LABEL_POSITIVE = "surface_glycoprotein"
LABEL_NEGATIVE = "other"
SCORE_COLUMN = "surface_glycoprotein_score"


def default_output_path(input_path):
    """Output CSV in the current directory, named after the input."""
    input_path = Path(input_path)
    stem = input_path.name if input_path.is_dir() else input_path.stem
    return Path.cwd() / f"{stem}.surface_glyco.csv"


def write_results(results, output_file):
    """Write prediction rows as CSV (ids containing commas or quotes are quoted)."""
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["id", "prediction", SCORE_COLUMN])
        for result in results:
            writer.writerow([result["id"], result["prediction"], f"{result[SCORE_COLUMN]:.4f}"])


def build_results(seq_ids, predictions, probabilities):
    """Return one result row per sequence (every score, not only positive calls)."""
    return [
        {
            "id": seq_id,
            "prediction": LABEL_POSITIVE if predictions[i] == 1 else LABEL_NEGATIVE,
            SCORE_COLUMN: probabilities[i][1],
        }
        for i, seq_id in enumerate(seq_ids)
    ]


def main(
    input_path, model_path, output_file, model_name, silent=False, show_all=False, max_workers=None
):
    """Run prediction on input sequences.

    Args:
        input_path: Path to a FASTA file or directory containing FASTA files.
        model_path: Path to trained model.
        output_file: Optional output CSV file path.
        model_name: ESM-2 model variant to use.
        silent: Suppress per-sequence output to stdout.
        show_all: Print all predictions to stdout, not just positive calls. The output
            CSV always contains every sequence.
        max_workers: Maximum number of workers for parallel file processing.
    """
    print("=" * 50)
    print("Surface glycoprotein scoring")
    print("=" * 50)

    print(f"Loading model from {model_path}...")
    classifier = load_model(model_path)
    card = load_model_card(model_path)
    if card is None:
        print("  No model card found; assuming the model was trained with --model-name embeddings")
    elif card.get("esm_model") != model_name:
        print(
            f"Error: model was trained on {card.get('esm_model')} embeddings, "
            f"but --model-name is {model_name}"
        )
        sys.exit(1)

    # Handle both file and directory inputs
    input_path = Path(input_path)
    if input_path.is_file():
        print(f"Processing single file: {input_path}")
        fasta_files = [input_path]
    elif input_path.is_dir():
        print(f"Finding FASTA files in {input_path}...")
        fasta_files = find_fasta_files(input_path)
    else:
        print(f"Error: {input_path} is neither a file nor a directory")
        sys.exit(1)

    if not fasta_files:
        print(f"No FASTA files found in {input_path}")
        sys.exit(1)

    print(f"Found {len(fasta_files)} FASTA file(s)")

    # Use parallel processing for multiple files, sequential for single files
    if len(fasta_files) > 1:
        print("Using parallel processing for multiple files...")
        all_sequences = process_fasta_files_parallel(fasta_files, max_workers=max_workers)
    else:
        print("Processing single file...")
        all_sequences = []
        sequences = process_fasta_file(fasta_files[0])
        all_sequences.extend(sequences)

    if not all_sequences:
        print("No sequences found in input files")
        sys.exit(1)

    print(f"Total sequences to predict: {len(all_sequences)}")

    print(f"Extracting ESM-2 embeddings using {model_name}...")
    embeddings, seq_ids = get_esm_embeddings(all_sequences, model_name=model_name)

    if len(embeddings) == 0:
        print("Error: No embeddings extracted")
        sys.exit(1)

    print("Making predictions...")
    predictions = predict(classifier, embeddings)
    probabilities = predict_proba(classifier, embeddings)

    results = build_results(seq_ids, predictions, probabilities)

    if not silent:
        print("\nResults:")
        print("-" * 50)
        for result in results:
            # only print positive calls by default; the output file always has every score
            if show_all or result["prediction"] == LABEL_POSITIVE:
                print(f"{result['id']}: {result['prediction']} " f"(p={result[SCORE_COLUMN]:.3f})")

    if output_file is None:
        output_file = default_output_path(input_path)

    print(f"\nSaving results to {output_file}...")
    write_results(results, output_file)

    print("=" * 50)
    print("Prediction complete!")
    print("=" * 50)


def cli():
    """Command-line interface entry point."""
    parser = argparse.ArgumentParser(
        description="Score proteins as fungal cell-surface glycoproteins"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Input FASTA file or directory with FASTA files",
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=None,
        help="Path to trained model (defaults to a name based on --model-name)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output CSV file for predictions (defaults to a name based on input, "
            "written to the current working directory)"
        ),
    )
    parser.add_argument(
        "--model-name",
        default=DEFAULT_MODEL,
        choices=ESM2_MODEL_CHOICES,
        help="ESM-2 model variant (must match training)",
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Suppress per-sequence scores on stdout",
    )
    parser.add_argument(
        "--show-all",
        action="store_true",
        help="Print all predictions to stdout (the output CSV always has every score)",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Maximum number of workers for parallel file processing (defaults to CPU count)",
    )

    args = parser.parse_args()

    if args.input is None:
        parser.print_help()
        sys.exit(1)

    if args.model is None:
        args.model = get_models_dir() / model_filename(args.model_name)

    if not args.model.exists():
        print(
            f"Error: model file not found at {args.model}. "
            "Train one with surface_glyco_train, or pass --model.",
            file=sys.stderr,
        )
        sys.exit(1)

    main(
        args.input,
        args.model,
        args.output,
        args.model_name,
        silent=args.silent,
        show_all=args.show_all,
        max_workers=args.max_workers,
    )


if __name__ == "__main__":
    cli()
