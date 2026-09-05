"""Tests for analysis/adhesion_properties/stats_ext.py."""

import numpy as np
import pandas as pd
from stats_ext import MIN_GROUP_N, aggregate_by_clade, mannwhitney_within_clade


def _clade_df(clade_values, feature_values):
    return pd.DataFrame({"phylum": clade_values, "length": feature_values})


def test_aggregate_by_clade_basic():
    df = _clade_df(
        ["PhylumX"] * 6 + ["PhylumY"] * 6,
        [100, 110, 120, 130, 140, 150, 200, 210, 220, 230, 240, 250],
    )
    summary = aggregate_by_clade(df, "phylum", "length")
    assert set(summary["phylum"]) == {"PhylumX", "PhylumY"}
    y_median = summary.loc[summary["phylum"] == "PhylumY", "median"].iloc[0]
    x_median = summary.loc[summary["phylum"] == "PhylumX", "median"].iloc[0]
    assert y_median > x_median
    assert not summary["small_n"].any()  # both groups have 6 >= MIN_GROUP_N


def test_aggregate_by_clade_flags_small_n():
    df = _clade_df(["PhylumZ"] * 2, [100, 110])
    summary = aggregate_by_clade(df, "phylum", "length")
    assert summary.iloc[0]["small_n"] is np.True_ or summary.iloc[0]["small_n"] is True
    assert 2 < MIN_GROUP_N


def test_mannwhitney_within_clade_detects_difference():
    rng = np.random.default_rng(1)
    adhesion = _clade_df(["PhylumX"] * 10, list(rng.normal(200, 5, 10)))
    background = _clade_df(["PhylumX"] * 10, list(rng.normal(100, 5, 10)))
    result = mannwhitney_within_clade(adhesion, background, "phylum", "length")
    assert len(result) == 1
    assert result.iloc[0]["p_value"] < 0.05
    assert "p_value_bh" in result.columns


def test_mannwhitney_within_clade_excludes_small_n_clades():
    adhesion = _clade_df(["PhylumSmall"] * 2, [100, 110])
    background = _clade_df(["PhylumSmall"] * 10, list(range(100, 110)))
    result = mannwhitney_within_clade(adhesion, background, "phylum", "length")
    assert result.empty


def test_mannwhitney_within_clade_empty_result_has_expected_columns():
    """No clade is eligible for testing (both below MIN_GROUP_N) -> result
    must be an empty-but-well-formed DataFrame, not a bare pd.DataFrame(),
    so downstream to_csv()/read_csv() round-trips without an EmptyDataError."""
    adhesion = _clade_df(["PhylumSmall"] * 2, [100, 110])
    background = _clade_df(["PhylumSmall"] * 2, [90, 95])
    result = mannwhitney_within_clade(adhesion, background, "phylum", "length")
    assert result.empty
    assert list(result.columns) == [
        "phylum",
        "n_adhesion",
        "n_background",
        "median_adhesion",
        "median_background",
        "u_stat",
        "p_value",
        "p_value_bh",
    ]
