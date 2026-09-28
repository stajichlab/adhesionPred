"""Tests for analysis/embedding_clustering/embed.py — exercises only the
pure bookkeeping logic (chunk save/load, protein-id reading), never
real model inference (no GPU/network needed for these tests)."""

import numpy as np
import pandas as pd
from embed import (
    load_all_embedding_chunks,
    read_protein_universe_adhesion_ids,
    save_embeddings_chunk,
)


def test_read_protein_universe_adhesion_ids(tmp_path):
    universe = tmp_path / "protein_universe.csv"
    pd.DataFrame(
        {
            "protein_id": ["P1", "P2", "P3"],
            "group": ["adhesion", "background", "adhesion"],
        }
    ).to_csv(universe, index=False)
    ids = read_protein_universe_adhesion_ids(universe)
    assert ids == ["P1", "P3"]


def test_save_and_load_embedding_chunks(tmp_path):
    chunk0 = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
    chunk1 = np.array([[5.0, 6.0]], dtype=np.float32)
    save_embeddings_chunk(chunk0, ["A", "B"], tmp_path, 0)
    save_embeddings_chunk(chunk1, ["C"], tmp_path, 1)

    embeddings, ids = load_all_embedding_chunks(tmp_path)
    assert embeddings.shape == (3, 2)
    assert ids == ["A", "B", "C"]
    np.testing.assert_array_equal(embeddings, np.vstack([chunk0, chunk1]))
