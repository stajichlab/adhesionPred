#!/usr/bin/env python
"""Cluster ESM-C 300M-space embeddings and compute a UMAP projection for
visualization."""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cluster import compute_umap_projection, reduce_and_cluster  # noqa: E402
from embed import load_all_embedding_chunks  # noqa: E402

EMBEDDINGS_DIR = Path(__file__).resolve().parent / "tables" / "esmc300m"
TABLES_DIR = Path(__file__).resolve().parent / "tables"


def main() -> None:
    embeddings, ids = load_all_embedding_chunks(EMBEDDINGS_DIR)
    print(f"Loaded {embeddings.shape[0]} embeddings, dim={embeddings.shape[1]}")

    # n_pca_components=15 (not the reduce_and_cluster default of 50): see
    # 03_cluster_esm2_classifier.py and task-6-report.md for the benchmark
    # showing HDBSCAN's near-O(n^~2) scaling at 50 PCA dims on ~750K points,
    # and the ~2.4x speedup measured at 15 dims -- applied identically here
    # for consistency, still clustering the FULL dataset (no subsampling).
    print("Starting PCA + HDBSCAN (reduce_and_cluster)...", flush=True)
    t0 = time.time()
    result = reduce_and_cluster(embeddings, n_pca_components=15)
    print(
        f"reduce_and_cluster done in {time.time() - t0:.1f}s: "
        f"{result['n_clusters']} clusters, {result['n_noise']} noise points",
        flush=True,
    )

    pd.DataFrame({"protein_id": ids, "cluster_label": result["cluster_labels"]}).to_csv(
        TABLES_DIR / "cluster_labels_esmc300m.csv", index=False
    )

    print("Starting UMAP projection...", flush=True)
    t0 = time.time()
    umap_projection = compute_umap_projection(embeddings)
    print(f"UMAP done in {time.time() - t0:.1f}s", flush=True)
    np.save(TABLES_DIR / "umap_esmc300m.npy", umap_projection)
    print("Wrote cluster labels and UMAP projection")


if __name__ == "__main__":
    main()
