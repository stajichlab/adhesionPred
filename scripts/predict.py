#!/usr/bin/env python
"""Prediction script for classifying adhesion proteins."""

import argparse
import sys
from pathlib import Path

from adhesion_predict.config import DEFAULT_MODEL, INPUT_DIR, MODELS_DIR
from adhesion_predict.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from adhesion_predict.io import find_fasta_files, process_fasta_file
from adhesion_predict.model import load_model, predict, predict_proba


def main(input_path, model_path, output_file, model_name):
    """Run prediction on input sequences.

    Args:
        input_path: Path to a FASTA file or directory containing FASTA files.
        model_path: Path to trained model.
        output_file: Optional output CSV file path.
        model_name: ESM-2 model variant to use.
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
        prob_adhesion = probabilities[i][1]
        results.append(
            {
                "id": seq_id,
                "prediction": pred_label,
                "probability_adhesion": prob_adhesion,
            }
        )

    print("\nResults:")
    print("-" * 50)
    for result in results:
        print(f"{result['id']}: {result['prediction']} (p={result['probability_adhesion']:.3f})")

    if output_file:
        print(f"\nSaving results to {output_file}...")
        with open(output_file, "w") as f:
            f.write("id,prediction,probability_adhesion\n")
            for result in results:
                f.write(
                    f"{result['id']},{result['prediction']},{result['probability_adhesion']:.4f}\n"
                )

    print("=" * 50)
    print("Prediction complete!")
    print("=" * 50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict adhesion proteins using trained model")
    parser.add_argument(
        "--input",
        type=Path,
        default=INPUT_DIR,
        help="Input FASTA file or directory with FASTA files",
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=MODELS_DIR / "adhesion_model.pkl",
        help="Path to trained model",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV file for predictions",
    )
    parser.add_argument(
        "--model-name",
        default=DEFAULT_MODEL,
        choices=ESM2_MODEL_CHOICES,
        help="ESM-2 model variant (must match training)",
    )

    args = parser.parse_args()

    if not args.model.exists():
        print(f"Error: Model file not found at {args.model}")
        print("Run training first: python scripts/train.py")
        sys.exit(1)

    main(args.input, args.model, args.output, args.model_name)
