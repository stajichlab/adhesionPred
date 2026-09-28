"""Per-cluster property/domain validation functions."""

import pandas as pd


def summarize_clusters(master_df: pd.DataFrame, cluster_col: str) -> pd.DataFrame:
    """Per-cluster N, median length/pct_ser_thr_pro, mean hydrophobicity,
    and domain-flag fractions (for whichever has_* columns are present in
    master_df)."""
    domain_flags = [c for c in master_df.columns if c.startswith("has_")]
    rows = []
    for cluster, gdf in master_df.groupby(cluster_col):
        row = {
            cluster_col: cluster,
            "n": len(gdf),
            "median_length": gdf["length"].median() if "length" in gdf else None,
            "median_pct_ser_thr_pro": (
                gdf["pct_ser_thr_pro"].median() if "pct_ser_thr_pro" in gdf else None
            ),
            "mean_hydrophobicity": (
                gdf["mean_hydrophobicity"].mean() if "mean_hydrophobicity" in gdf else None
            ),
        }
        for flag in domain_flags:
            row[flag] = gdf[flag].mean()
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n", ascending=False).reset_index(drop=True)


def cluster_contingency_table(labels_a: pd.Series, labels_b: pd.Series) -> pd.DataFrame:
    """Cross-tabulation between two clusterings of the same proteins
    (same row order/index correspondence assumed)."""
    return pd.crosstab(labels_a, labels_b)


def check_feature_concentration(
    master_df: pd.DataFrame, cluster_col: str, flag_col: str
) -> pd.DataFrame:
    """Per-cluster fraction with flag_col set, sorted descending -- used to
    check whether a flagged feature (e.g. has_aa1_family) concentrates in
    one cluster or is scattered across all of them."""
    result = master_df.groupby(cluster_col)[flag_col].agg(["mean", "size"]).reset_index()
    result = result.rename(columns={"mean": "fraction", "size": "n"})
    return result.sort_values("fraction", ascending=False).reset_index(drop=True)
