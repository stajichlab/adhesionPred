"""PCA + HDBSCAN clustering, and UMAP projection, for embedding-space
protein clustering. Uses sklearn.cluster.HDBSCAN (built into scikit-learn
>=1.3, already a project dependency) rather than the separate hdbscan
PyPI package."""

import numpy as np
import umap
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA


def reduce_and_cluster(
    embeddings: np.ndarray, n_pca_components: int = 50, min_cluster_size: int = 50
) -> dict:
    """PCA-reduce embeddings, then HDBSCAN-cluster. Returns cluster
    labels (HDBSCAN's -1 = noise, reported as its own category, never
    discarded), the PCA-reduced embeddings, and summary counts."""
    n_components = min(n_pca_components, embeddings.shape[0] - 1, embeddings.shape[1])
    pca = PCA(n_components=n_components, random_state=42)
    pca_embeddings = pca.fit_transform(embeddings)

    clusterer = HDBSCAN(min_cluster_size=min_cluster_size)
    cluster_labels = clusterer.fit_predict(pca_embeddings)

    unique_labels = set(cluster_labels)
    n_clusters = len(unique_labels - {-1})
    n_noise = int(np.sum(cluster_labels == -1))

    return {
        "cluster_labels": cluster_labels,
        "pca_embeddings": pca_embeddings,
        "n_clusters": n_clusters,
        "n_noise": n_noise,
    }


def compute_umap_projection(
    embeddings: np.ndarray, n_neighbors: int = 15, random_state: int = 42
) -> np.ndarray:
    """2D UMAP projection of embeddings, for visualization only (not used
    for clustering itself — that happens on the PCA-reduced space)."""
    reducer = umap.UMAP(n_neighbors=n_neighbors, n_components=2, random_state=random_state)
    return reducer.fit_transform(embeddings)
