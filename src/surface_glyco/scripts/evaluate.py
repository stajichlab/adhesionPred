#!/usr/bin/env python
"""Evaluation script for assessing model performance."""

import argparse
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score

from surface_glyco.card import ModelCardError, resolve_embedding_settings
from surface_glyco.config import (
    DEFAULT_MODEL,
    NEGATIVE_DIR,
    POSITIVE_DIR,
    get_models_dir,
    model_filename,
)
from surface_glyco.embeddings import ESM2_MODEL_CHOICES, get_esm_embeddings
from surface_glyco.io import load_sequences_from_dir
from surface_glyco.model import (
    load_model,
    load_model_card,
    predict,
    predict_proba,
    require_model_file,
)


def main(positive_dir, negative_dir, model_path, model_name=None):
    """Evaluate model on test data. model_name None uses the model card's variant."""
    print("=" * 50)
    print("Model Evaluation")
    print("=" * 50)

    try:
        settings = resolve_embedding_settings(
            load_model_card(model_path), model_name, DEFAULT_MODEL
        )
    except ModelCardError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

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

    print(f"Extracting embeddings using {settings.esm_model}...")
    embeddings, seq_ids, kept = get_esm_embeddings(
        all_sequences,
        model_name=settings.esm_model,
        repr_layer=settings.repr_layer,
        return_indices=True,
    )

    if len(embeddings) == 0:
        print("Error: No embeddings extracted")
        sys.exit(1)

    labels = np.array([all_sequences[k]["label"] for k in kept])

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
    target_names = ["other", "surface_glycoprotein"]
    print(classification_report(labels, predictions, target_names=target_names))

    try:
        auc = roc_auc_score(labels, probabilities[:, 1])
        print(f"ROC-AUC: {auc:.3f}")
    except ValueError:
        print("Could not calculate ROC-AUC (check class distribution)")

    accuracy = np.mean(predictions == labels)
    print(f"\nOverall Accuracy: {accuracy:.3f}")

    print("=" * 50)


def cli():
    """Command-line interface entry point."""
    parser = argparse.ArgumentParser(description="Evaluate a surface glycoprotein model")
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
        default=None,
        help="Path to trained model (defaults to a name based on --model-name)",
    )
    parser.add_argument(
        "--model-name",
        default=None,
        choices=ESM2_MODEL_CHOICES,
        help="ESM-2 model variant; defaults to the model card's, and must match it if given",
    )

    args = parser.parse_args()

    if args.model is None:
        args.model = get_models_dir() / model_filename(args.model_name or DEFAULT_MODEL)

    require_model_file(args.model)

    main(args.positive, args.negative, args.model, args.model_name)


if __name__ == "__main__":
    cli()
