"""Figure generation for the kingdom-wide adhesion survey report.

Every group figure annotates each group with its N and visually
desaturates groups with N < MIN_GROUP_N, since sample-size imbalance
across fungal phyla/orders is the main confound in this analysis and
needs to be visible in the figure itself, not only in a table.

Colors follow the dataviz skill's validated, colorblind-safe default
palette (see .superpowers dataviz skill, references/palette.md):
- Full-N groups use categorical slot 1 (blue); small-N groups use the
  muted axis/label gray, so "insufficient sample size" reads as a
  de-emphasized/warning-like state rather than just another hue.
- The scatter plot's phylum coloring uses the fixed 8-slot categorical
  order (assigned deterministically, never cycled/re-ordered per call);
  a scatter plot is an "all-pairs" context under the skill's validator,
  so any phylum beyond the first 8 folds into the muted gray "Other"
  bucket rather than reusing/cycling hues.
"""

from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

MIN_GROUP_N = 5

# dataviz skill validated palette (references/palette.md)
_FULL_COLOR = "#2a78d6"  # categorical slot 1 (blue)
_SMALL_N_COLOR = "#898781"  # muted axis/label gray
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
_OTHER_COLOR = "#898781"  # muted gray, for overflow past the 8 validated slots


def _order_groups(df: pd.DataFrame, rank_col: str, value_col: str, top_n: Optional[int]) -> list:
    medians = df.groupby(rank_col)[value_col].median().sort_values(ascending=False)
    if top_n is not None:
        counts = df[rank_col].value_counts()
        top_groups = counts.sort_values(ascending=False).head(top_n).index
        medians = medians[medians.index.isin(top_groups)]
    return list(medians.index)


def boxplot_by_rank(
    df: pd.DataFrame,
    rank_col: str,
    value_col: str,
    out_path: Path,
    top_n: Optional[int] = None,
    title: str = "",
) -> None:
    """Ordered box plot of value_col across rank_col groups, N-annotated, small-N greyed."""
    sub = df[df[rank_col].notna()].copy()
    order = _order_groups(sub, rank_col, value_col, top_n)
    sub = sub[sub[rank_col].isin(order)]
    counts = sub[rank_col].value_counts()

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.5), 6))
    palette = {g: (_FULL_COLOR if counts[g] >= MIN_GROUP_N else _SMALL_N_COLOR) for g in order}
    sns.boxplot(
        data=sub,
        x=rank_col,
        y=value_col,
        order=order,
        hue=rank_col,
        palette=palette,
        legend=False,
        ax=ax,
    )
    labels = [f"{g} (n={counts[g]})" for g in order]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, rotation=90 if len(order) > 8 else 45, ha="right")
    ax.set_title(title)
    ax.set_xlabel(rank_col.capitalize())
    ax.set_ylabel(value_col.replace("_", " "))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _phylum_color_map(phyla: list) -> dict:
    """Deterministic, fixed-order categorical color assignment.

    A scatter plot is an "all-pairs" context (every point can sit next to
    every other), so only the palette's first 8 validated categorical
    slots are used; any phylum beyond that folds into a shared muted-gray
    "Other" bucket instead of cycling colors.
    """
    ordered_phyla = sorted(phyla)
    color_map = {}
    for i, p in enumerate(ordered_phyla):
        color_map[p] = _CATEGORICAL_ORDER[i] if i < len(_CATEGORICAL_ORDER) else _OTHER_COLOR
    return color_map


def scatter_proteome_vs_adhesion(df: pd.DataFrame, out_path: Path) -> None:
    """Scatter of total_proteins vs adhesion_count, colored by phylum, with a
    reference line for 'expected count under the kingdom-wide median
    fraction' — separates genome-size effects from real enrichment."""
    sub = df[df["phylum"].notna()].copy()
    kingdom_median_fraction = sub["adhesion_fraction"].median()
    color_map = _phylum_color_map(sub["phylum"].unique())

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.scatterplot(
        data=sub,
        x="total_proteins",
        y="adhesion_count",
        hue="phylum",
        palette=color_map,
        alpha=0.6,
        s=20,
        ax=ax,
        legend="brief",
    )
    x_range = [sub["total_proteins"].min(), sub["total_proteins"].max()]
    ax.plot(
        x_range,
        [x * kingdom_median_fraction for x in x_range],
        color="black",
        linestyle="--",
        linewidth=1,
        label=f"expected @ kingdom median fraction ({kingdom_median_fraction:.4f})",
    )
    ax.set_xlabel("Total proteins in proteome")
    ax.set_ylabel("Adhesion-called proteins")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def probability_boxplot_by_rank(
    df: pd.DataFrame, rank_col: str, out_path: Path, top_n: Optional[int] = None
) -> None:
    """Distribution of mean_adhesion_prob by rank_col — secondary signal of
    'how confidently adhesive', same N-annotation/greying as boxplot_by_rank."""
    boxplot_by_rank(
        df,
        rank_col,
        "mean_adhesion_prob",
        out_path,
        top_n=top_n,
        title=f"Mean adhesion probability by {rank_col}",
    )


def ranked_summary_table(summary_df: pd.DataFrame, rank_col: str) -> str:
    """Render a sorted markdown table (replaces a min-max-normalized heatmap,
    which would visually inflate trivial differences over ~9 phyla).

    Small-N groups (small_n is True, i.e. n < MIN_GROUP_N) are marked with
    a dagger on the group name so the table itself flags "descriptive
    only, excluded from formal tests" without needing an extra column.
    """
    cols = [rank_col, "n_species", "median_fraction", "median_count", "median_prob"]
    sorted_df = summary_df.sort_values("median_fraction", ascending=False)[cols + ["small_n"]]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for _, row in sorted_df.iterrows():
        values = []
        for c in cols:
            v = row[c]
            if c == rank_col and row["small_n"]:
                values.append(f"{v} †")
            else:
                values.append(f"{v:.4f}" if isinstance(v, float) else str(v))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)
