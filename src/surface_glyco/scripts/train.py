#!/usr/bin/env python
"""Training script for the surface glycoprotein classifier."""

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import sklearn

from surface_glyco.config import (
    DEFAULT_MODEL,
    DEFAULT_TEST_SIZE,
    MODELS_DIR,
    NEGATIVE_DIR,
    POSITIVE_DIR,
)
from surface_glyco.embeddings import (
    DEFAULT_REPR_LAYER,
    ESM2_MODEL_CHOICES,
    MAX_RESIDUES,
    POOLING,
    get_esm_embeddings,
)
from surface_glyco.io import load_sequences_from_dir
from surface_glyco.model import save_model, save_model_card, train_classifier


def sequences_sha256(sequences):
    """Order-independent hash of the labelled training sequences."""
    h = hashlib.sha256()
    for line in sorted(f"{s['label']}\t{s['id']}\t{s['sequence']}" for s in sequences):
        h.update(line.encode())
        h.update(b"\n")
    return h.hexdigest()


def dedupe_sequences(sequences):
    """Drop records whose (label, sequence) was already seen; keep the first id.

    Exact duplicates inflate the positive class and leak across the train/test split.
    """
    seen = set()
    unique = []
    for seq in sequences:
        key = (seq["label"], seq["sequence"])
        if key not in seen:
            seen.add(key)
            unique.append(seq)
    return unique


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

    unique = dedupe_sequences(all_sequences)
    n_removed = len(all_sequences) - len(unique)
    if n_removed:
        print(f"  Removed {n_removed} exact-duplicate sequences ({len(unique)} remain)")
    return unique, n_removed


def main(positive_dir, negative_dir, output_model, model_name, test_size):
    """Main training pipeline using ESM-2 embeddings."""
    print("=" * 50)
    print("Surface glycoprotein classifier training")
    print("=" * 50)

    sequences, n_duplicates_removed = prepare_data(positive_dir, negative_dir)

    embeddings, seq_ids, kept = get_esm_embeddings(
        sequences, model_name=model_name, return_indices=True
    )

    if len(embeddings) == 0:
        print("Error: No embeddings extracted", file=sys.stderr)
        sys.exit(1)

    print(f"  Embedding shape: {embeddings.shape}")

    labels = np.array([sequences[k]["label"] for k in kept])

    classifier, test_acc = train_classifier(embeddings, labels, test_size=test_size)

    save_model(classifier, output_model)
    save_model_card(
        output_model,
        {
            "esm_model": model_name,
            "repr_layer": DEFAULT_REPR_LAYER,
            "pooling": POOLING,
            "max_residues": MAX_RESIDUES,
            "classifier": type(classifier).__name__,
            "sklearn_version": sklearn.__version__,
            "n_positive": int(labels.sum()),
            "n_negative": int(len(labels) - labels.sum()),
            "n_duplicates_removed": n_duplicates_removed,
            "n_not_embedded": len(sequences) - len(kept),
            "training_sequences_sha256": sequences_sha256(sequences),
            "positive_dir": str(positive_dir),
            "negative_dir": str(negative_dir),
            "holdout_accuracy": float(test_acc),
        },
    )

    print("=" * 50)
    print("Training complete!")
    print("=" * 50)


def cli():
    """Command-line interface entry point."""
    parser = argparse.ArgumentParser(
        description="Train a surface glycoprotein classifier using ESM-2 embeddings"
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
        default=None,
        help="Output path for trained model (defaults to a name based on --model)",
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

    if args.output is None:
        args.output = MODELS_DIR / f"adhesion_model_{args.model}.pkl"

    args.output.parent.mkdir(parents=True, exist_ok=True)

    main(args.positive, args.negative, args.output, args.model, args.test_size)


if __name__ == "__main__":
    cli()
