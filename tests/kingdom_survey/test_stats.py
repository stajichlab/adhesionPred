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


def test_mixedlm_by_rank_handles_class_reserved_keyword():
    """Regression test: rank_col="class" is a Python/patsy reserved keyword.

    mixedlm_by_rank interpolates rank_col into a patsy formula string
    (C(rank_col)); without Q(...) quoting, patsy tries to parse the bare
    token "class" as Python code and raises a SyntaxError. This must run
    without raising for the real "class" taxonomic rank.
    """
    df = _sample_df()
    result = mixedlm_by_rank(df, "class", "adhesion_fraction")
    assert "converged" in result
    assert result["converged"] in (True, False)


def test_dunn_posthoc_returns_empty_when_fewer_than_two_groups_eligible():
    df = _sample_df()
    single_group = df[df["phylum"] == "PhylumX"]
    posthoc = dunn_posthoc(single_group, "phylum", "adhesion_fraction")
    assert isinstance(posthoc, pd.DataFrame)
    assert posthoc.empty


def test_mixedlm_by_rank_degenerate_input_returns_uniform_keys():
    df = _sample_df()
    single_rank_value = df[df["phylum"] == "PhylumX"]
    result = mixedlm_by_rank(single_rank_value, "phylum", "adhesion_fraction")
    expected_keys = {
        "rank",
        "value_col",
        "converged",
        "genus_variance",
        "residual_variance",
        "llf",
        "summary",
    }
    assert set(result) == expected_keys
    assert result["converged"] is False
    assert result["genus_variance"] is None
    assert result["residual_variance"] is None
    assert result["llf"] is None
    assert isinstance(result["summary"], str)
