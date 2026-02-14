"""Feature extraction functions for protein sequences."""

import numpy as np


def extract_sequence_features(sequence):
    """Extract compositional features from a protein sequence.

    Args:
        sequence: Protein sequence string.

    Returns:
        List of features (20 amino acid frequencies + normalized length).
    """
    seq = sequence.upper()
    length = len(seq)
    seq.replace("J", "L")  # Treat J as L for simplicity
    if length == 0:
        return [0] * 21

    aa_counts = {
        "A": 0,
        "R": 0,
        "N": 0,
        "D": 0,
        "C": 0,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 0,
        "I": 0,
        "L": 0,
        "K": 0,
        "M": 0,
        "F": 0,
        "P": 0,
        "S": 0,
        "T": 0,
        "W": 0,
        "Y": 0,
        "V": 0,
    }

    for aa in seq:
        if aa in aa_counts:
            aa_counts[aa] += 1

    features = [aa_counts[aa] / length for aa in sorted(aa_counts.keys())]
    features.append(length / 1000.0)

    return features


def prepare_features(sequences):
    """Convert sequences to feature matrix.

    Args:
        sequences: List of sequence dictionaries.

    Returns:
        numpy array of features.
    """
    x = []
    for seq in sequences:
        features = extract_sequence_features(seq["sequence"])
        x.append(features)
    return np.array(x)
