"""Smoke tests for analysis/kingdom_survey/figures.py — verifies each
function runs without error and writes a non-empty file; does not
inspect pixel content."""

import pandas as pd
from figures import (
    boxplot_by_rank,
    probability_boxplot_by_rank,
    ranked_summary_table,
    scatter_proteome_vs_adhesion,
)
from stats import aggregate_summary


def _sample_df():
    rows = []
    for i in range(8):
        rows.append(
            {
                "phylum": "PhylumX" if i < 4 else "PhylumY",
                "adhesion_fraction": 0.1 + 0.01 * i,
                "adhesion_count": 10 + i,
                "total_proteins": 1000 + 50 * i,
                "mean_adhesion_prob": 0.6 + 0.01 * i,
            }
        )
    return pd.DataFrame(rows)


def test_boxplot_by_rank_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "box.png"
    boxplot_by_rank(df, "phylum", "adhesion_fraction", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "scatter.png"
    scatter_proteome_vs_adhesion(df, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_probability_boxplot_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "prob.png"
    probability_boxplot_by_rank(df, "phylum", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_ranked_summary_table_is_markdown():
    df = _sample_df()
    summary = aggregate_summary(df, "phylum")
    md = ranked_summary_table(summary, "phylum")
    assert "PhylumX" in md
    assert "|" in md
