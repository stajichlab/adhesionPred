#!/usr/bin/env python
"""Training script for adhesion protein classifier."""

import argparse
import sys
from pathlib import Path

import numpy as np
from adhesion_predict.config import (
    DEFAULT_MODEL,
    DEFAULT_TEST_SIZE,
    MODELS_DIR,
    NEGATIVE_DIR,
    POSITIVE_DIR,
)
from adhesion_predict.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from adhesion_predict.io import load_sequences_from_dir
from adhesion_predict.model import save_model, train_classifier


def prepare_data(positive_dir, negative_dir):
    """Load and label sequences from positive and negative directories."""
    print("Loading sequences...")
    positive_seqs = load_sequences_from_dir(positive_dir)
    negative_seqs = load_sequences_from_dir(negative_dir)

    print(f"  Positive samples: {len(positive_seqs)}")
    print(f"  Negative samples: {len(negative_seqs)}")

    if not positive_seqs:
        print(f"Error: No sequences found in {positive_dir}", file=sys.stderr)
        sys.exit(1)
    if not negative_seqs:
        print(f"Error: No sequences found in {negative_dir}", file=sys.stderr)
        sys.exit(1)

    all_sequences = []
    for seq in positive_seqs:
        seq["label"] = 1
        all_sequences.append(seq)
    for seq in negative_seqs:
        seq["label"] = 0
        all_sequences.append(seq)

    return all_sequences


def main(positive_dir, negative_dir, output_model, model_name, test_size):
    """Main training pipeline using ESM-2 embeddings."""
    print("=" * 50)
    print("Adhesion Protein Classification - LLM Training")
    print("=" * 50)

    sequences = prepare_data(positive_dir, negative_dir)

    embeddings, seq_ids = get_esm_embeddings(sequences, model_name=model_name)

    if len(embeddings) == 0:
        print("Error: No embeddings extracted", file=sys.stderr)
        sys.exit(1)

    print(f"  Embedding shape: {embeddings.shape}")

    labels = np.array([seq["label"] for seq in sequences[: len(embeddings)]])

    classifier, test_acc = train_classifier(embeddings, labels, test_size=test_size)

    save_model(classifier, output_model)

    print("=" * 50)
    print("Training complete!")
    print("=" * 50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train adhesion protein classifier using ESM-2 embeddings"
    )
    parser.add_argument(
        "--positive",
        type=Path,
        default=POSITIVE_DIR,
        help="Directory with positive training sequences",
    )
    parser.add_argument(
        "--negative",
        type=Path,
        default=NEGATIVE_DIR,
        help="Directory with negative training sequences",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=MODELS_DIR / "adhesion_model.pkl",
        help="Output path for trained model",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        choices=ESM2_MODEL_CHOICES,
        help="ESM-2 model variant",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=DEFAULT_TEST_SIZE,
        help="Proportion of data for test set",
    )

    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    main(args.positive, args.negative, args.output, args.model, args.test_size)
