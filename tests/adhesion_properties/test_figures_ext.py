"""Smoke tests for analysis/adhesion_properties/figures_ext.py — verifies
each function runs without error and writes a non-empty file."""

import pandas as pd
from figures_ext import (
    domain_fraction_barplot,
    paired_boxplot_by_rank,
    scatter_feature_vs_probability,
    scatter_species_mean_length_vs_fraction,
    top_domain_enrichment_barplot,
)


def _master_df():
    rows = []
    for i in range(12):
        rows.append(
            {
                "phylum": "PhylumX" if i < 6 else "PhylumY",
                "group": "adhesion" if i % 2 == 0 else "background",
                "length": 200 + 10 * i,
                "probability_adhesion": 0.5 + 0.02 * i,
                "locustag": f"LOC{i % 3}",
            }
        )
    return pd.DataFrame(rows)


def test_paired_boxplot_by_rank_creates_file(tmp_path):
    df = _master_df()
    out = tmp_path / "box.png"
    paired_boxplot_by_rank(df, "phylum", "length", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_domain_fraction_barplot_creates_file(tmp_path):
    fractions = pd.DataFrame(
        [
            {"phylum": "PhylumX", "group": "adhesion", "n": 6, "has_pfam": 0.5},
            {"phylum": "PhylumX", "group": "background", "n": 6, "has_pfam": 0.3},
            {"phylum": "PhylumY", "group": "adhesion", "n": 6, "has_pfam": 0.6},
            {"phylum": "PhylumY", "group": "background", "n": 6, "has_pfam": 0.4},
        ]
    )
    out = tmp_path / "bar.png"
    domain_fraction_barplot(fractions, "phylum", "has_pfam", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_top_domain_enrichment_barplot_creates_file(tmp_path):
    top = pd.DataFrame(
        [
            {"domain_id": "DomainA", "enrichment_ratio": 3.5},
            {"domain_id": "DomainB", "enrichment_ratio": float("inf")},
        ]
    )
    out = tmp_path / "top.png"
    top_domain_enrichment_barplot(top, out, title="Test")
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_feature_vs_probability_creates_file(tmp_path):
    df = _master_df()
    adhesion = df[df["group"] == "adhesion"]
    out = tmp_path / "scatter.png"
    scatter_feature_vs_probability(adhesion, "length", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_species_mean_length_vs_fraction_creates_file(tmp_path):
    df = _master_df()
    adhesion = df[df["group"] == "adhesion"]
    kingdom_summary = pd.DataFrame(
        {"locustag": ["LOC0", "LOC1", "LOC2"], "adhesion_fraction": [0.01, 0.02, 0.03]}
    )
    out = tmp_path / "species_scatter.png"
    scatter_species_mean_length_vs_fraction(adhesion, kingdom_summary, out)
    assert out.exists()
    assert out.stat().st_size > 0
