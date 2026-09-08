"""Tests for analysis/embedding_clustering/cluster.py — uses small
synthetic embeddings with two obviously-separated Gaussian blobs, never
the real 749,697x480/960 data."""

import numpy as np
from cluster import compute_umap_projection, reduce_and_cluster


def _two_blob_embeddings(n_per_blob=100, dim=20, seed=1):
    rng = np.random.default_rng(seed)
    blob_a = rng.normal(loc=0.0, scale=0.5, size=(n_per_blob, dim))
    blob_b = rng.normal(loc=20.0, scale=0.5, size=(n_per_blob, dim))
    return np.vstack([blob_a, blob_b]).astype(np.float32)


def test_reduce_and_cluster_finds_two_separated_blobs():
    embeddings = _two_blob_embeddings()
    result = reduce_and_cluster(embeddings, n_pca_components=10, min_cluster_size=20)
    assert result["n_clusters"] >= 2
    labels = result["cluster_labels"]
    assert len(labels) == len(embeddings)
    # the two blobs should mostly not share a cluster label (allowing some noise, label -1)
    blob_a_labels = set(labels[:100]) - {-1}
    blob_b_labels = set(labels[100:]) - {-1}
    assert not (
        blob_a_labels & blob_b_labels
    ), "the two separated blobs were merged into one cluster"


def test_reduce_and_cluster_reports_noise_as_its_own_category():
    embeddings = _two_blob_embeddings()
    result = reduce_and_cluster(embeddings, n_pca_components=10, min_cluster_size=20)
    assert "n_noise" in result
    assert result["n_noise"] >= 0


def test_compute_umap_projection_shape():
    embeddings = _two_blob_embeddings(n_per_blob=30)
    projection = compute_umap_projection(embeddings, n_neighbors=5)
    assert projection.shape == (60, 2)
