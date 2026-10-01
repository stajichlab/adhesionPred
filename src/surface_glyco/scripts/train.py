#!/usr/bin/env python
"""Training script for the surface glycoprotein classifier."""

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np

from surface_glyco.card import DEFAULT_REPR_LAYER, POOLING_RESIDUE_MEAN, new_card
from surface_glyco.config import (
    DEFAULT_MODEL,
    DEFAULT_TEST_SIZE,
    NEGATIVE_DIR,
    POSITIVE_DIR,
    model_filename,
)
from surface_glyco.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from surface_glyco.io import load_sequences_from_dir
from surface_glyco.model import save_model, save_model_card, train_classifier


def sequences_sha256(sequences):
    """Order-independent hash of the labelled training sequences."""
    h = hashlib.sha256()
    for line in sorted(f"{s['label']}\t{s['id']}\t{s['sequence']}" for s in sequences):
        h.update(line.encode())
        h.update(b"\n")
    return h.hexdigest()


def dedupe_with_counts(sequences):
    """Return (unique, counts).

    Empty sequences are counted first. A conflict is a non-empty sequence present with both
    labels; every copy is dropped, because it cannot be learned and leaks across any split.
    What remains are exact duplicates within a class; the first id is kept. The three counts
    and len(unique) sum to len(sequences).
    """
    non_empty = [s for s in sequences if s["sequence"]]
    labels_by_seq = {}
    for seq in non_empty:
        labels_by_seq.setdefault(seq["sequence"], set()).add(seq["label"])
    seen = set()
    unique = []
    n_conflicting = 0
    for seq in non_empty:
        if len(labels_by_seq[seq["sequence"]]) > 1:
            n_conflicting += 1
            continue
        key = (seq["label"], seq["sequence"])
        if key not in seen:
            seen.add(key)
            unique.append(seq)
    counts = {
        "n_empty_removed": len(sequences) - len(non_empty),
        "n_conflicting_removed": n_conflicting,
        "n_duplicates_removed": len(non_empty) - n_conflicting - len(unique),
    }
    return unique, counts


def dedupe_sequences(sequences):
    """Drop empty sequences, sequences present in both classes, and in-class duplicates."""
    return dedupe_with_counts(sequences)[0]


def _environment():
    """Versions of the libraries a model's scores depend on."""
    from importlib.metadata import PackageNotFoundError, version

    out = {}
    for key, dist in (
        ("numpy", "numpy"),
        ("scikit_learn", "scikit-learn"),
        ("torch", "torch"),
        ("fair_esm", "fair-esm"),
    ):
        try:
            out[key] = version(dist)
        except PackageNotFoundError:
            out[key] = "unknown"
    return out


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

    unique, counts = dedupe_with_counts(all_sequences)
    n_removed = len(all_sequences) - len(unique)
    if n_removed:
        print(
            f"  Removed {counts['n_empty_removed']} empty, "
            f"{counts['n_conflicting_removed']} in-both-classes and "
            f"{counts['n_duplicates_removed']} duplicate sequences ({len(unique)} remain)"
        )
    return unique, counts


def main(positive_dir, negative_dir, output_model, model_name, test_size):
    """Main training pipeline using ESM-2 embeddings."""
    print("=" * 50)
    print("Surface glycoprotein classifier training")
    print("=" * 50)

    sequences, counts = prepare_data(positive_dir, negative_dir)

    embeddings, seq_ids, kept = get_esm_embeddings(
        sequences,
        model_name=model_name,
        repr_layer=DEFAULT_REPR_LAYER,
        return_indices=True,
    )

    if len(embeddings) == 0:
        print("Error: No embeddings extracted", file=sys.stderr)
        sys.exit(1)

    print(f"  Embedding shape: {embeddings.shape}")

    labels = np.array([sequences[k]["label"] for k in kept])

    classifier, test_acc = train_classifier(embeddings, labels, test_size=test_size)

    save_model(classifier, output_model)
    n_positive = int(labels.sum())
    save_model_card(
        output_model,
        new_card(
            model_name,
            DEFAULT_REPR_LAYER,
            POOLING_RESIDUE_MEAN,
            n_positive,
            int(len(labels) - n_positive),
            classifier=type(classifier).__name__,
            **counts,
            n_not_embedded=len(sequences) - len(kept),
            training_sequences_sha256=sequences_sha256(sequences),
            positive_dir=str(positive_dir),
            negative_dir=str(negative_dir),
            holdout_accuracy=float(test_acc),
            environment=_environment(),
        ),
    )

    print("=" * 50)
    print("Training complete!")
    print("=" * 50)


def default_output_path(esm_model):
    """Default path for a newly trained model: ./models in the working directory.

    The installed package directory is not used because it is read-only in a normal install.
    """
    return Path.cwd() / "models" / model_filename(esm_model)


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
        args.output = default_output_path(args.model)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    main(args.positive, args.negative, args.output, args.model, args.test_size)


if __name__ == "__main__":
    cli()
