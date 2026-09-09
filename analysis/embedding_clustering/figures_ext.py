"""Figure generation for the embedding-space clustering analysis.

Colors follow the dataviz skill's guidance and reuse
`analysis/kingdom_survey/figures.py`'s validated categorical palette
where a categorical encoding is appropriate; a perceptually-uniform
continuous colormap (viridis) is used for a continuous feature like
`pct_ser_thr_pro`.

`umap_scatter_by_cluster` must work for both embedding spaces this
project produces, which have very different cluster counts:
ESM2-classifier space has only 2 real clusters + noise, while ESM-C
300M space has ~512 real clusters + noise. Giving every one of 512
clusters its own distinct color and legend entry (the naive approach,
e.g. cycling a 20-color tab20 colormap) would not be useful: colors
would repeat every 20 clusters (visually implying false identity
between unrelated clusters), the legend would run to hundreds of
entries, and the figure would not communicate any real structure. So
this module caps the number of individually-colored/legended clusters
at `MAX_DISTINCT_CLUSTERS` (ranked by cluster size, largest first);
any additional real clusters are folded into a single muted-gray
"other clusters (N)" bucket. HDBSCAN noise (-1) is always kept as its
own separate, distinct light-gray entry -- it is not "just another
cluster" and must not be confused with the "other clusters" bucket.
This same logic runs unconditionally for both embedding spaces: when
the real cluster count is <= MAX_DISTINCT_CLUSTERS (as in the
ESM2-classifier space), every cluster still gets its own distinct
color and legend entry and the "other clusters" bucket is simply
never populated -- there is no dataset-specific fork.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Cap on how many real clusters get their own distinct color + legend
# entry in umap_scatter_by_cluster. Chosen to keep a tab20-based legend
# legible (roughly one full pass through tab20's 20 colors) while still
# highlighting real structure for datasets with hundreds of clusters.
MAX_DISTINCT_CLUSTERS = 15

# dataviz-skill-validated categorical/gray tones, reused from
# analysis/kingdom_survey/figures.py so this project's figures share a
# consistent palette rather than inventing a new one.
_NOISE_COLOR = "lightgray"
_OTHER_COLOR = "#898781"  # muted gray, validated "Other"/de-emphasis tone
_CATEGORICAL_ORDER = [
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#4a3aa7",  # 7 violet
    "#e34948",  # 8 red
]
# tab20 indices 14/15 are matplotlib's built-in gray pair -- visually
# indistinguishable from _OTHER_COLOR/_NOISE_COLOR above. Skipped so a
# tab20-fallback cluster color never collides with the "other clusters"
# bucket or the noise color (a real collision found in whole-branch review).
_TAB20_GRAY_INDICES = {14, 15}


def _tab20_fallback_color(slot: int):
    """0-indexed fallback slot -> a tab20 color, skipping the gray pair."""
    tab20 = plt.get_cmap("tab20")
    idx = 0
    seen = 0
    while True:
        if (idx % 20) not in _TAB20_GRAY_INDICES:
            if seen == slot:
                return tab20(idx % 20)
            seen += 1
        idx += 1


def umap_scatter_by_cluster(
    umap_coords: np.ndarray, cluster_labels: np.ndarray, out_path: Path, title: str
) -> None:
    """2D UMAP scatter colored by cluster label.

    Noise (-1) is rendered in light gray as its own distinct legend
    entry. Real clusters are ranked by size (largest first); only the
    top `MAX_DISTINCT_CLUSTERS` get a distinct color + legend entry
    (drawn from a fixed categorical order, falling back to `tab20` for
    the tail of that top-N when there are more than 8), and any
    remaining smaller clusters are folded into a single muted-gray
    "other clusters (N)" bucket. When there are <= MAX_DISTINCT_CLUSTERS
    real clusters, every cluster keeps its own entry and the "other"
    bucket is never drawn -- this logic runs unconditionally, so it
    degrades gracefully for small cluster counts and folds correctly
    for large ones without any dataset-specific branching.
    """
    cluster_labels = np.asarray(cluster_labels)
    fig, ax = plt.subplots(figsize=(8, 7))

    real_labels = [label for label in set(cluster_labels.tolist()) if label != -1]
    sizes = {label: int((cluster_labels == label).sum()) for label in real_labels}
    ranked = sorted(real_labels, key=lambda label: sizes[label], reverse=True)
    top_labels = ranked[:MAX_DISTINCT_CLUSTERS]
    overflow_labels = ranked[MAX_DISTINCT_CLUSTERS:]

    for i, label in enumerate(sorted(top_labels, key=lambda label: sizes[label], reverse=True)):
        color = (
            _CATEGORICAL_ORDER[i]
            if i < len(_CATEGORICAL_ORDER)
            else _tab20_fallback_color(i - len(_CATEGORICAL_ORDER))
        )
        mask = cluster_labels == label
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            s=3,
            alpha=0.5,
            color=color,
            label=f"cluster {label} (n={sizes[label]})",
        )

    if overflow_labels:
        overflow_mask = np.isin(cluster_labels, overflow_labels)
        ax.scatter(
            umap_coords[overflow_mask, 0],
            umap_coords[overflow_mask, 1],
            s=3,
            alpha=0.3,
            color=_OTHER_COLOR,
            # Explicitly distinguish cluster count from protein count here --
            # every other legend entry's "n=" is a protein count, so a bare
            # "(n={len(overflow_labels)})" (a cluster count) reads as a
            # protein count by comparison (found in whole-branch review).
            label=f"other clusters ({len(overflow_labels)} clusters, n={int(overflow_mask.sum())} proteins)",
            zorder=1,
        )

    noise_mask = cluster_labels == -1
    if noise_mask.any():
        ax.scatter(
            umap_coords[noise_mask, 0],
            umap_coords[noise_mask, 1],
            s=3,
            alpha=0.3,
            color=_NOISE_COLOR,
            label=f"noise (n={int(noise_mask.sum())})",
            zorder=0,
        )

    ax.set_title(title)
    ax.legend(markerscale=4, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def umap_scatter_by_feature(
    umap_coords: np.ndarray,
    feature_values: np.ndarray,
    out_path: Path,
    title: str,
    is_categorical: bool = False,
) -> None:
    """2D UMAP scatter colored by a continuous (viridis) or boolean/
    categorical (fixed categorical palette) feature, for visually
    checking whether clusters align with known properties. Feature
    cardinality is expected to stay low (booleans/small categoricals
    like has_signal_peptide, has_cazy), so no top-N-plus-other folding
    is needed here."""
    fig, ax = plt.subplots(figsize=(8, 7))
    if is_categorical:
        unique_values = sorted(set(feature_values.tolist()), key=str)
        for i, value in enumerate(unique_values):
            color = (
                _CATEGORICAL_ORDER[i]
                if i < len(_CATEGORICAL_ORDER)
                else _tab20_fallback_color(i - len(_CATEGORICAL_ORDER))
            )
            mask = feature_values == value
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                s=3,
                alpha=0.5,
                color=color,
                label=str(value),
            )
        ax.legend(markerscale=4, fontsize=8)
    else:
        scatter = ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            s=3,
            alpha=0.5,
            c=feature_values,
            cmap="viridis",
        )
        fig.colorbar(scatter, ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
