"""Tests for analysis/kingdom_survey/stats.py."""

import numpy as np
import pandas as pd
from stats import (
    aggregate_summary,
    dunn_posthoc,
    genus_average,
    kruskal_wallis_by_rank,
    mixedlm_by_rank,
)


def _sample_df():
    rows = []
    rng = np.random.default_rng(42)
    for _ in range(12):
        rows.append(
            {
                "genus": "GenusA",
                "phylum": "PhylumX",
                "class": "ClassX",
                "order": "OrderX",
                "family": "FamilyX",
                "adhesion_fraction": 0.10 + rng.normal(0, 0.01),
                "adhesion_count": 10,
                "mean_adhesion_prob": 0.7,
            }
        )
    for _ in range(12):
        rows.append(
            {
                "genus": "GenusB",
                "phylum": "PhylumY",
                "class": "ClassY",
                "order": "OrderY",
                "family": "FamilyY",
                "adhesion_fraction": 0.30 + rng.normal(0, 0.01),
                "adhesion_count": 30,
                "mean_adhesion_prob": 0.8,
            }
        )
    return pd.DataFrame(rows)


def test_aggregate_summary():
    df = _sample_df()
    summary = aggregate_summary(df, "phylum")
    assert set(summary["phylum"]) == {"PhylumX", "PhylumY"}
    y_frac = summary.loc[summary["phylum"] == "PhylumY", "median_fraction"].iloc[0]
    x_frac = summary.loc[summary["phylum"] == "PhylumX", "median_fraction"].iloc[0]
    assert y_frac > x_frac


def test_kruskal_wallis_by_rank_detects_difference():
    df = _sample_df()
    result = kruskal_wallis_by_rank(df, "phylum", "adhesion_fraction")
    assert result["n_groups"] == 2
    assert result["p_value"] < 0.05


def test_kruskal_wallis_excludes_small_groups():
    df = _sample_df()
    small = pd.DataFrame(
        [
            {
                "genus": "GenusC",
                "phylum": "PhylumZ",
                "class": "ClassZ",
                "order": "OrderZ",
                "family": "FamilyZ",
                "adhesion_fraction": 0.99,
                "adhesion_count": 1,
                "mean_adhesion_prob": 0.99,
            }
        ]
    )
    df = pd.concat([df, small], ignore_index=True)
    result = kruskal_wallis_by_rank(df, "phylum", "adhesion_fraction")
    assert result["n_groups"] == 2  # PhylumZ (n=1) excluded, not counted


def test_dunn_posthoc_shape():
    df = _sample_df()
    posthoc = dunn_posthoc(df, "phylum", "adhesion_fraction")
    assert "PhylumX" in posthoc.index
    assert "PhylumY" in posthoc.columns


def test_genus_average_collapses_species():
    df = _sample_df()
    averaged = genus_average(df)
    assert len(averaged) == 2
    assert set(averaged["n_species"]) == {12}


def test_mixedlm_by_rank_runs():
    df = _sample_df()
    result = mixedlm_by_rank(df, "phylum", "adhesion_fraction")
    assert result["converged"] is True
    assert "genus_variance" in result
