"""Smoke tests for analysis/embedding_clustering/figures_ext.py."""

import numpy as np
from figures_ext import (
    MAX_DISTINCT_CLUSTERS,
    umap_scatter_by_cluster,
    umap_scatter_by_feature,
)


def _sample_umap(n=50):
    rng = np.random.default_rng(1)
    return rng.normal(size=(n, 2))


def test_umap_scatter_by_cluster_creates_file(tmp_path):
    coords = _sample_umap()
    labels = np.array([0, 1, -1] * 16 + [0, 1])
    out = tmp_path / "umap_cluster.png"
    umap_scatter_by_cluster(coords, labels, out, title="Test")
    assert out.exists()
    assert out.stat().st_size > 0


def test_umap_scatter_by_feature_continuous(tmp_path):
    coords = _sample_umap()
    values = np.linspace(0, 100, 50)
    out = tmp_path / "umap_feature.png"
    umap_scatter_by_feature(coords, values, out, title="Test", is_categorical=False)
    assert out.exists()
    assert out.stat().st_size > 0


def test_umap_scatter_by_feature_categorical(tmp_path):
    coords = _sample_umap()
    values = np.array([True, False] * 25)
    out = tmp_path / "umap_bool.png"
    umap_scatter_by_feature(coords, values, out, title="Test", is_categorical=True)
    assert out.exists()
    assert out.stat().st_size > 0


def test_umap_scatter_by_cluster_many_clusters_folds_to_other(tmp_path, monkeypatch):
    """Real ESM-C 300M data has ~512 real clusters + noise. Plotting one
    legend entry per cluster would be unreadable (tab20 cycling every 20
    colors, legend hundreds of entries long) and wouldn't communicate real
    structure. umap_scatter_by_cluster must instead give distinct legend
    entries only to the largest MAX_DISTINCT_CLUSTERS clusters and fold
    everything else into a single "other clusters" bucket, with noise (-1)
    kept as its own separate gray entry -- so the legend stays bounded
    regardless of how many real clusters exist."""
    n_clusters = 120
    rng = np.random.default_rng(2)
    n = 5000
    coords = rng.normal(size=(n, 2))
    # Give clusters very unequal sizes so "largest N" is well-defined:
    # cluster i gets size proportional to (n_clusters - i), plus some noise.
    weights = np.arange(n_clusters, 0, -1, dtype=float)
    weights = weights / weights.sum()
    cluster_choices = rng.choice(n_clusters, size=int(n * 0.9), p=weights)
    noise_choices = np.full(n - len(cluster_choices), -1)
    labels = np.concatenate([cluster_choices, noise_choices])
    rng.shuffle(labels)

    out = tmp_path / "umap_many_clusters.png"

    recorded_labels = []
    import figures_ext as fx

    real_scatter = fx.plt.Axes.scatter

    def _spy_scatter(self, *args, **kwargs):
        if "label" in kwargs:
            recorded_labels.append(kwargs["label"])
        return real_scatter(self, *args, **kwargs)

    monkeypatch.setattr(fx.plt.Axes, "scatter", _spy_scatter)

    umap_scatter_by_cluster(coords, labels, out, title="Many clusters")

    assert out.exists()
    assert out.stat().st_size > 0

    # Legend entries must stay bounded: at most MAX_DISTINCT_CLUSTERS named
    # clusters + 1 "other clusters" bucket + 1 "noise" entry -- never one
    # entry per real cluster (which would be ~120 here, or ~512 for the
    # real ESM-C 300M data).
    assert len(recorded_labels) <= MAX_DISTINCT_CLUSTERS + 2
    assert any("other" in label.lower() for label in recorded_labels)
    assert any("noise" in label.lower() for label in recorded_labels)


def test_umap_scatter_by_cluster_small_cluster_count_stays_distinct(tmp_path, monkeypatch):
    """When the number of real clusters is small (e.g. ESM2-classifier
    space's 2 clusters), every cluster should still get its own distinct
    legend entry -- the top-N-plus-other logic must degrade gracefully to
    "one entry per cluster" rather than always folding into "other"."""
    coords = _sample_umap()
    labels = np.array([0, 1, -1] * 16 + [0, 1])
    out = tmp_path / "umap_small.png"

    recorded_labels = []
    import figures_ext as fx

    real_scatter = fx.plt.Axes.scatter

    def _spy_scatter(self, *args, **kwargs):
        if "label" in kwargs:
            recorded_labels.append(kwargs["label"])
        return real_scatter(self, *args, **kwargs)

    monkeypatch.setattr(fx.plt.Axes, "scatter", _spy_scatter)

    umap_scatter_by_cluster(coords, labels, out, title="Small clusters")

    assert out.exists()
    # 2 real clusters + noise, no "other" bucket needed.
    assert not any("other" in label.lower() for label in recorded_labels)
    assert any("noise" in label.lower() for label in recorded_labels)
    assert len(recorded_labels) == 3
