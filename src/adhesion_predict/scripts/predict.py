#!/usr/bin/env python
"""Prediction script for classifying adhesion proteins."""

import argparse
import sys
from pathlib import Path

from adhesion_predict.config import DEFAULT_MODEL, MODELS_DIR
from adhesion_predict.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from adhesion_predict.io import find_fasta_files, process_fasta_file
from adhesion_predict.model import load_model, predict, predict_proba


def main(input_path, model_path, output_file, model_name, silent=False, show_all=False):
    """Run prediction on input sequences.

    Args:
        input_path: Path to a FASTA file or directory containing FASTA files.
        model_path: Path to trained model.
        output_file: Optional output CSV file path.
        model_name: ESM-2 model variant to use.
        silent: Suppress per-sequence output to stdout.
        show_all: Show all predictions, not just adhesion proteins.
    """
    print("=" * 50)
    print("Adhesion Protein Prediction")
    print("=" * 50)

    print(f"Loading model from {model_path}...")
    classifier = load_model(model_path)

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

    all_sequences = []
    for fasta_file in fasta_files:
        sequences = process_fasta_file(fasta_file)
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

    results = []
    for i, seq_id in enumerate(seq_ids):
        pred_label = "Adhesion" if predictions[i] == 1 else "Non-adhesion"
        if show_all or pred_label == "Adhesion":
            # only print adhesion predictions by default, but can show all if requested
            prob_adhesion = probabilities[i][1]
            results.append(
                {
                    "id": seq_id,
                    "prediction": pred_label,
                    "probability_adhesion": prob_adhesion,
                }
            )

    if not silent:
        print("\nResults:")
        print("-" * 50)
        for result in results:
            print(
                f"{result['id']}: {result['prediction']} (p={result['probability_adhesion']:.3f})"
            )

    if output_file is None:
        output_dir = Path.cwd()
        if input_path.is_dir():
            output_file = output_dir / f"{input_path.name}.adhesion_predict.csv"
        else:
            output_file = output_dir / f"{input_path.stem}.adhesion_predict.csv"

    print(f"\nSaving results to {output_file}...")
    with open(output_file, "w") as f:
        f.write("id,prediction,probability_adhesion\n")
        for result in results:
            f.write(f"{result['id']},{result['prediction']},{result['probability_adhesion']:.4f}\n")

    print("=" * 50)
    print("Prediction complete!")
    print("=" * 50)


def cli():
    """Command-line interface entry point."""
    parser = argparse.ArgumentParser(description="Predict adhesion proteins using trained model")
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
        help="Show all predictions (default: False)",
    )

    args = parser.parse_args()

    if args.model is None:
        args.model = MODELS_DIR / f"adhesion_model_{args.model_name}.pkl"

    if not args.model.exists():
        print(f"Error: Model file not found at {args.model}")
        print("Run training first: adhesion_train")
        sys.exit(1)

    main(
        args.input,
        args.model,
        args.output,
        args.model_name,
        silent=args.silent,
        show_all=args.show_all,
    )


if __name__ == "__main__":
    cli()
