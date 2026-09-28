"""Taxonomic aggregation and statistics for the kingdom-wide adhesion survey.

Runs three parallel views at each rank, per the design's taxonomy-based
pseudo-phylogenetic correction: naive per-species tests (kruskal_wallis_by_rank/
dunn_posthoc directly on the input df), genus-averaged tests (same functions
called on genus_average(df)'s output), and a genus-random-intercept mixed
model (mixedlm_by_rank). None of this uses real branch lengths.
"""

import pandas as pd
import scikit_posthocs as sp
import statsmodels.formula.api as smf
from scipy.stats import kruskal

MIN_GROUP_N = 5


def aggregate_summary(df: pd.DataFrame, rank_col: str) -> pd.DataFrame:
    """Per-group N, median/IQR of adhesion_fraction and adhesion_count."""
    sub = df[df[rank_col].notna()]
    rows = []
    for group, gdf in sub.groupby(rank_col):
        rows.append(
            {
                rank_col: group,
                "n_species": len(gdf),
                "median_fraction": gdf["adhesion_fraction"].median(),
                "iqr_fraction": (
                    gdf["adhesion_fraction"].quantile(0.75)
                    - gdf["adhesion_fraction"].quantile(0.25)
                ),
                "median_count": gdf["adhesion_count"].median(),
                "iqr_count": (
                    gdf["adhesion_count"].quantile(0.75) - gdf["adhesion_count"].quantile(0.25)
                ),
                "median_prob": gdf["mean_adhesion_prob"].median(),
                "small_n": len(gdf) < MIN_GROUP_N,
            }
        )
    return pd.DataFrame(rows).sort_values("median_fraction", ascending=False).reset_index(drop=True)


def kruskal_wallis_by_rank(df: pd.DataFrame, rank_col: str, value_col: str) -> dict:
    """Kruskal-Wallis omnibus test across groups with N >= MIN_GROUP_N.

    Reports epsilon-squared effect size alongside p-value always — with
    ~5,800 species and wildly unbalanced phylum sample sizes, p-value alone
    is not a reliable signal of a biologically meaningful difference.
    """
    sub = df[df[rank_col].notna()]
    counts = sub[rank_col].value_counts()
    eligible_groups = counts[counts >= MIN_GROUP_N].index
    samples = [sub.loc[sub[rank_col] == g, value_col].to_numpy() for g in eligible_groups]

    if len(samples) < 2:
        return {
            "rank": rank_col,
            "value_col": value_col,
            "h_stat": float("nan"),
            "p_value": float("nan"),
            "epsilon_squared": float("nan"),
            "n_groups": len(samples),
            "n_total": sum(len(s) for s in samples),
        }

    h_stat, p_value = kruskal(*samples)
    n_total = sum(len(s) for s in samples)
    k = len(samples)
    # True epsilon-squared (Kruskal-Wallis effect size): H / (n_total - 1).
    # (Note: (H - k + 1) / (n_total - k) is eta-squared-from-H, a related
    # but distinct statistic — not epsilon-squared.)
    epsilon_sq = h_stat / (n_total - 1) if n_total > 1 else float("nan")
    return {
        "rank": rank_col,
        "value_col": value_col,
        "h_stat": h_stat,
        "p_value": p_value,
        "epsilon_squared": epsilon_sq,
        "n_groups": k,
        "n_total": n_total,
    }


def dunn_posthoc(df: pd.DataFrame, rank_col: str, value_col: str) -> pd.DataFrame:
    """Dunn's post-hoc pairwise comparisons (BH-corrected) among N>=MIN_GROUP_N groups.

    Returns an empty DataFrame when fewer than 2 groups meet MIN_GROUP_N —
    no post-hoc comparison is possible in that case.
    """
    sub = df[df[rank_col].notna()].copy()
    counts = sub[rank_col].value_counts()
    eligible = counts[counts >= MIN_GROUP_N].index
    sub = sub[sub[rank_col].isin(eligible)]
    if len(eligible) < 2:
        return pd.DataFrame()
    return sp.posthoc_dunn(sub, val_col=value_col, group_col=rank_col, p_adjust="fdr_bh")


def genus_average(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per genus: mean of numeric metrics, first taxonomy seen.

    This is the "genus-averaged" leg of the taxonomy-based pseudo-
    phylogenetic correction: it removes the dominant pseudoreplication
    source (a genus with many sequenced species dominating a coarser
    rank's median).
    """
    sub = df[df["genus"].notna()]
    agg = (
        sub.groupby("genus")
        .agg(
            phylum=("phylum", "first"),
            class_=("class", "first"),
            order=("order", "first"),
            family=("family", "first"),
            n_species=("adhesion_fraction", "size"),
            adhesion_fraction=("adhesion_fraction", "mean"),
            adhesion_count=("adhesion_count", "mean"),
            mean_adhesion_prob=("mean_adhesion_prob", "mean"),
        )
        .reset_index()
    )
    return agg.rename(columns={"class_": "class"})


def mixedlm_by_rank(df: pd.DataFrame, rank_col: str, value_col: str = "adhesion_fraction") -> dict:
    """Fit value_col ~ C(rank_col) with genus as a random intercept.

    Genus (not a fully nested order/family/genus structure) is used as
    the single random-effects level: statsmodels' MixedLM does not
    cleanly support 3-4 levels of nesting without variance-component
    machinery, and genus is the finest per-species grouping, capturing
    the largest share of non-independence. Still an approximation of a
    true phylogenetic comparative method (equal "branch lengths").
    """
    sub = df[df[rank_col].notna() & df["genus"].notna()].copy()
    counts = sub[rank_col].value_counts()
    eligible_groups = counts[counts >= MIN_GROUP_N].index
    sub = sub[sub[rank_col].isin(eligible_groups)]
    if sub[rank_col].nunique() < 2 or sub["genus"].nunique() < 2:
        return {
            "rank": rank_col,
            "value_col": value_col,
            "converged": False,
            "genus_variance": None,
            "residual_variance": None,
            "llf": None,
            "summary": "did not converge: fewer than 2 groups at this rank or fewer than 2 genera",
        }

    # Q(...) quoting is required: rank_col may be "class", a Python reserved
    # word, which breaks patsy's formula parser when interpolated bare into
    # C(...) (patsy parses the C(...) argument as a Python expression).
    model = smf.mixedlm(f"{value_col} ~ C(Q('{rank_col}'))", sub, groups=sub["genus"])
    fit = model.fit(reml=False)
    return {
        "rank": rank_col,
        "value_col": value_col,
        "converged": bool(fit.converged),
        "genus_variance": float(fit.cov_re.iloc[0, 0]),
        "residual_variance": float(fit.scale),
        "llf": float(fit.llf),
        "summary": fit.summary().as_text(),
    }
