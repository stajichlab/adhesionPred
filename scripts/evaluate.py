#!/usr/bin/env python
"""Evaluation script for assessing model performance."""

import argparse
import sys
from pathlib import Path

import numpy as np
from adhesion_predict.config import DEFAULT_MODEL, MODELS_DIR, NEGATIVE_DIR, POSITIVE_DIR
from adhesion_predict.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from adhesion_predict.io import load_sequences_from_dir
from adhesion_predict.model import load_model, predict, predict_proba
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score


def main(positive_dir, negative_dir, model_path, model_name):
    """Evaluate model on test data."""
    print("=" * 50)
    print("Model Evaluation")
    print("=" * 50)

    print("Loading sequences...")
    positive_seqs = load_sequences_from_dir(positive_dir)
    negative_seqs = load_sequences_from_dir(negative_dir)

    for seq in positive_seqs:
        seq["label"] = 1
    for seq in negative_seqs:
        seq["label"] = 0

    all_sequences = positive_seqs + negative_seqs
    print(f"  Total sequences: {len(all_sequences)}")

    print(f"Loading model from {model_path}...")
    classifier = load_model(model_path)

    print(f"Extracting embeddings using {model_name}...")
    embeddings, seq_ids = get_esm_embeddings(all_sequences, model_name=model_name)

    if len(embeddings) == 0:
        print("Error: No embeddings extracted")
        sys.exit(1)

    labels = np.array([seq["label"] for seq in all_sequences[: len(embeddings)]])

    predictions = predict(classifier, embeddings)
    probabilities = predict_proba(classifier, embeddings)

    print("\n" + "=" * 50)
    print("Results")
    print("=" * 50)

    print("\nConfusion Matrix:")
    cm = confusion_matrix(labels, predictions)
    print(f"  TN: {cm[0][0]:4d}  FP: {cm[0][1]:4d}")
    print(f"  FN: {cm[1][0]:4d}  TP: {cm[1][1]:4d}")

    print("\nClassification Report:")
    target_names = ["Non-adhesion", "Adhesion"]
    print(classification_report(labels, predictions, target_names=target_names))

    try:
        auc = roc_auc_score(labels, probabilities[:, 1])
        print(f"ROC-AUC: {auc:.3f}")
    except ValueError:
        print("Could not calculate ROC-AUC (check class distribution)")

    accuracy = np.mean(predictions == labels)
    print(f"\nOverall Accuracy: {accuracy:.3f}")

    print("=" * 50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate trained adhesion protein model")
    parser.add_argument(
        "--positive",
        type=Path,
        default=POSITIVE_DIR,
        help="Directory with positive test sequences",
    )
    parser.add_argument(
        "--negative",
        type=Path,
        default=NEGATIVE_DIR,
        help="Directory with negative test sequences",
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=MODELS_DIR / "adhesion_model.pkl",
        help="Path to trained model",
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
        sys.exit(1)

    main(args.positive, args.negative, args.model, args.model_name)
