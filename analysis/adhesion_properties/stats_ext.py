"""Per-clade statistics for the adhesion-protein properties analysis.

Reuses kingdom_survey's MIN_GROUP_N/kruskal_wallis_by_rank directly
(already generic over an arbitrary value_col) rather than duplicating
them; adds a generic per-clade summary (kingdom_survey's aggregate_summary
is hardcoded to its own adhesion_fraction/adhesion_count column names, so
is not reusable here) and a within-clade adhesion-vs-background test,
which kingdom_survey has no equivalent of.
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

_KINGDOM_SURVEY_DIR = Path(__file__).resolve().parent.parent / "kingdom_survey"
sys.path.insert(0, str(_KINGDOM_SURVEY_DIR))
from stats import MIN_GROUP_N, kruskal_wallis_by_rank  # noqa: E402,F401


def aggregate_by_clade(df: pd.DataFrame, rank_col: str, value_col: str) -> pd.DataFrame:
    """Per-clade N, median, IQR of an arbitrary value_col."""
    sub = df[df[rank_col].notna()]
    rows = []
    for group, gdf in sub.groupby(rank_col):
        rows.append(
            {
                rank_col: group,
                "n": len(gdf),
                "median": gdf[value_col].median(),
                "iqr": gdf[value_col].quantile(0.75) - gdf[value_col].quantile(0.25),
                "small_n": len(gdf) < MIN_GROUP_N,
            }
        )
    return pd.DataFrame(rows).sort_values("median", ascending=False).reset_index(drop=True)


def mannwhitney_within_clade(
    adhesion_df: pd.DataFrame, background_df: pd.DataFrame, rank_col: str, value_col: str
) -> pd.DataFrame:
    """For each clade with both groups meeting MIN_GROUP_N, Mann-Whitney U
    comparing adhesion vs. background value_col; BH-corrected across clades."""
    a_counts = adhesion_df[rank_col].value_counts()
    b_counts = background_df[rank_col].value_counts()
    eligible = sorted(
        set(a_counts[a_counts >= MIN_GROUP_N].index) & set(b_counts[b_counts >= MIN_GROUP_N].index)
    )

    rows = []
    for group in eligible:
        a_vals = adhesion_df.loc[adhesion_df[rank_col] == group, value_col].dropna()
        b_vals = background_df.loc[background_df[rank_col] == group, value_col].dropna()
        if len(a_vals) < 2 or len(b_vals) < 2:
            continue
        u_stat, p_value = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        rows.append(
            {
                rank_col: group,
                "n_adhesion": len(a_vals),
                "n_background": len(b_vals),
                "median_adhesion": a_vals.median(),
                "median_background": b_vals.median(),
                "u_stat": u_stat,
                "p_value": p_value,
            }
        )

    result = pd.DataFrame(rows)
    if not result.empty:
        _, corrected, _, _ = multipletests(result["p_value"], method="fdr_bh")
        result["p_value_bh"] = corrected
    return result
