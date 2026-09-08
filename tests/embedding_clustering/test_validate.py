"""Tests for analysis/embedding_clustering/validate.py."""

import pandas as pd
from validate import (
    check_feature_concentration,
    cluster_contingency_table,
    summarize_clusters,
)


def _master_df():
    return pd.DataFrame(
        {
            "protein_id": [f"P{i}" for i in range(8)],
            "cluster_label": [0, 0, 0, 0, 1, 1, 1, 1],
            "length": [100, 110, 120, 130, 400, 410, 420, 430],
            "pct_ser_thr_pro": [40, 42, 38, 44, 10, 12, 9, 11],
            "has_pfam": [False, False, True, False, True, True, True, True],
        }
    )


def test_summarize_clusters():
    summary = summarize_clusters(_master_df(), "cluster_label")
    assert set(summary["cluster_label"]) == {0, 1}
    cluster0 = summary[summary["cluster_label"] == 0].iloc[0]
    assert cluster0["n"] == 4
    assert cluster0["median_length"] == 115.0


def test_cluster_contingency_table():
    labels_a = pd.Series([0, 0, 1, 1])
    labels_b = pd.Series(["X", "X", "Y", "Y"])
    table = cluster_contingency_table(labels_a, labels_b)
    assert table.loc[0, "X"] == 2
    assert table.loc[1, "Y"] == 2


def test_check_feature_concentration():
    df = _master_df()
    df["has_aa1"] = [False, False, False, False, True, True, True, False]
    result = check_feature_concentration(df, "cluster_label", "has_aa1")
    cluster1_frac = result.loc[result["cluster_label"] == 1, "fraction"].iloc[0]
    assert cluster1_frac == 0.75
