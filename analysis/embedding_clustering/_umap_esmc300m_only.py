#!/usr/bin/env python
"""Recovery script: compute ONLY the UMAP projection for the ESM-C 300M
embedding space and write tables/umap_esmc300m.npy.

Used after 04_cluster_esmc300m.py's HDBSCAN step completed and was written
to tables/cluster_labels_esmc300m.csv, but the UMAP step was OOM-killed by
the interactive SLURM job's 16G memory cgroup (ESM-C's raw embeddings are
960-dim vs. ESM2-classifier's 480-dim, and compute_umap_projection operates
on the raw embeddings, not the PCA-reduced ones -- roughly double the
nearest-neighbor-graph memory footprint). Re-running only this step (not the
~2.5h HDBSCAN step, whose output is already on disk and intact) under a
dedicated higher-memory SLURM allocation.

Relies on load_all_embedding_chunks' deterministic (sorted glob) ordering
matching the ids order 04_cluster_esmc300m.py used when it wrote
cluster_labels_esmc300m.csv, so the row order here lines up with that CSV.
Verifies this explicitly before saving.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cluster import compute_umap_projection  # noqa: E402
from embed import load_all_embedding_chunks  # noqa: E402

EMBEDDINGS_DIR = Path(__file__).resolve().parent / "tables" / "esmc300m"
TABLES_DIR = Path(__file__).resolve().parent / "tables"


def main() -> None:
    embeddings, ids = load_all_embedding_chunks(EMBEDDINGS_DIR)
    print(f"Loaded {embeddings.shape[0]} embeddings, dim={embeddings.shape[1]}", flush=True)

    existing = pd.read_csv(TABLES_DIR / "cluster_labels_esmc300m.csv")
    assert list(existing["protein_id"]) == list(ids), (
        "id order mismatch between existing cluster_labels_esmc300m.csv and "
        "a fresh load_all_embedding_chunks call -- UMAP rows would not line "
        "up with the HDBSCAN cluster labels"
    )
    print("id order verified consistent with cluster_labels_esmc300m.csv", flush=True)

    umap_projection = compute_umap_projection(embeddings)
    print(f"UMAP shape: {umap_projection.shape}", flush=True)
    np.save(TABLES_DIR / "umap_esmc300m.npy", umap_projection)
    print("Wrote umap_esmc300m.npy", flush=True)


if __name__ == "__main__":
    main()
