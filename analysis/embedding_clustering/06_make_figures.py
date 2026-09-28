#!/usr/bin/env python
"""Generate UMAP figures for both embedding spaces: colored by cluster,
and colored by key features (pct_ser_thr_pro, has_signal_peptide,
has_cazy) for visual alignment checking."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures_ext import umap_scatter_by_cluster, umap_scatter_by_feature  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"
ADHESION_PROPERTIES_TABLES = Path(__file__).resolve().parents[1] / "adhesion_properties" / "tables"


def _make_space_figures(space_name: str) -> None:
    umap_coords = np.load(TABLES_DIR / f"umap_{space_name}.npy")
    clusters = pd.read_csv(TABLES_DIR / f"cluster_labels_{space_name}.csv")
    props = pd.read_csv(ADHESION_PROPERTIES_TABLES / "protein_sequence_properties.csv")
    flags = pd.read_csv(ADHESION_PROPERTIES_TABLES / "protein_domain_flags.csv")
    merged = clusters.merge(props, on="protein_id", how="left").merge(
        flags, on="protein_id", how="left"
    )

    umap_scatter_by_cluster(
        umap_coords,
        merged["cluster_label"].to_numpy(),
        FIGURES_DIR / f"umap_cluster_{space_name}.png",
        title=f"Clusters ({space_name})",
    )
    umap_scatter_by_feature(
        umap_coords,
        merged["pct_ser_thr_pro"].to_numpy(),
        FIGURES_DIR / f"umap_pct_ser_thr_pro_{space_name}.png",
        title=f"Ser/Thr/Pro % ({space_name})",
        is_categorical=False,
    )
    umap_scatter_by_feature(
        umap_coords,
        merged["has_signal_peptide"].to_numpy(),
        FIGURES_DIR / f"umap_has_signal_peptide_{space_name}.png",
        title=f"Signal peptide ({space_name})",
        is_categorical=True,
    )
    umap_scatter_by_feature(
        umap_coords,
        merged["has_cazy"].to_numpy(),
        FIGURES_DIR / f"umap_has_cazy_{space_name}.png",
        title=f"CAZy presence ({space_name})",
        is_categorical=True,
    )


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for space_name in ["esm2_classifier", "esmc300m"]:
        _make_space_figures(space_name)
    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
