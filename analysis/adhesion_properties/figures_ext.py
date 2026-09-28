"""Figure generation for the adhesion-protein properties analysis.

Reuses kingdom_survey/figures.py's validated color constants (dataviz
skill palette) for consistency across both reports: "adhesion" = full
color, "background" = the same muted gray kingdom_survey uses for its
de-emphasized/reference states.
"""

import sys
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

_KINGDOM_SURVEY_DIR = Path(__file__).resolve().parent.parent / "kingdom_survey"
sys.path.insert(0, str(_KINGDOM_SURVEY_DIR))
from figures import _FULL_COLOR, _OTHER_COLOR  # noqa: E402

_GROUP_PALETTE = {"adhesion": _FULL_COLOR, "background": _OTHER_COLOR}


def paired_boxplot_by_rank(
    master_df: pd.DataFrame,
    rank_col: str,
    value_col: str,
    out_path: Path,
    top_n: Optional[int] = None,
) -> None:
    """Boxplot of value_col by rank_col, adhesion vs. background as paired
    boxes per clade (hue='group'), ordered by adhesion-group median."""
    sub = master_df[master_df[rank_col].notna()].copy()
    adhesion_medians = (
        sub[sub["group"] == "adhesion"]
        .groupby(rank_col)[value_col]
        .median()
        .sort_values(ascending=False)
    )
    if top_n is not None:
        counts = sub[sub["group"] == "adhesion"][rank_col].value_counts()
        top_groups = counts.sort_values(ascending=False).head(top_n).index
        adhesion_medians = adhesion_medians[adhesion_medians.index.isin(top_groups)]
    order = list(adhesion_medians.index)
    sub = sub[sub[rank_col].isin(order)]

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.6), 6))
    sns.boxplot(
        data=sub,
        x=rank_col,
        y=value_col,
        hue="group",
        order=order,
        palette=_GROUP_PALETTE,
        ax=ax,
    )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90 if len(order) > 8 else 45, ha="right")
    ax.set_xlabel(rank_col.capitalize())
    ax.set_ylabel(value_col.replace("_", " "))
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def domain_fraction_barplot(
    fractions_df: pd.DataFrame, rank_col: str, flag_col: str, out_path: Path
) -> None:
    """Grouped bar chart of domain-presence fraction by clade, adhesion vs
    background."""
    order = (
        fractions_df[fractions_df["group"] == "adhesion"]
        .sort_values(flag_col, ascending=False)[rank_col]
        .tolist()
    )
    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.6), 6))
    sns.barplot(
        data=fractions_df,
        x=rank_col,
        y=flag_col,
        hue="group",
        order=order,
        palette=_GROUP_PALETTE,
        ax=ax,
    )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90 if len(order) > 8 else 45, ha="right")
    ax.set_xlabel(rank_col.capitalize())
    ax.set_ylabel(f"Fraction with {flag_col.replace('has_', '').replace('_', ' ')}")
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def top_domain_enrichment_barplot(
    top_df: pd.DataFrame, out_path: Path, title: str, ratio_cap: float = 50.0
) -> None:
    """Horizontal bar chart of enrichment_ratio for the top domains,
    capping infinite/very large ratios at ratio_cap for display (noted in
    the axis label, not silently)."""
    sub = top_df.copy()
    sub["display_ratio"] = sub["enrichment_ratio"].clip(upper=ratio_cap)
    sub = sub.sort_values("display_ratio")

    fig, ax = plt.subplots(figsize=(8, max(4, len(sub) * 0.35)))
    ax.barh(sub["domain_id"], sub["display_ratio"], color=_FULL_COLOR)
    ax.set_xlabel(f"Adhesion enrichment ratio (capped at {ratio_cap:.0f})")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def scatter_feature_vs_probability(
    adhesion_df: pd.DataFrame, feature_col: str, out_path: Path
) -> None:
    """Scatter of feature_col vs. probability_adhesion, adhesion proteins only."""
    sub = adhesion_df.dropna(subset=[feature_col, "probability_adhesion"])
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(sub[feature_col], sub["probability_adhesion"], alpha=0.3, s=10, color=_FULL_COLOR)
    ax.set_xlabel(feature_col.replace("_", " "))
    ax.set_ylabel("Adhesion probability")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def scatter_species_mean_length_vs_fraction(
    adhesion_df: pd.DataFrame, kingdom_summary_df: pd.DataFrame, out_path: Path
) -> None:
    """Scatter of per-species mean adhesion-protein length vs. that
    species' adhesion_fraction (from kingdom_survey's summary table)."""
    mean_length = adhesion_df.groupby("locustag")["length"].mean().rename("mean_length")
    merged = kingdom_summary_df.set_index("locustag").join(mean_length, how="inner")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        merged["mean_length"], merged["adhesion_fraction"], alpha=0.4, s=15, color=_FULL_COLOR
    )
    ax.set_xlabel("Mean adhesion-protein length (species)")
    ax.set_ylabel("Species adhesion fraction")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
