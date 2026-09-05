"""Smoke tests for analysis/kingdom_survey/figures.py — verifies each
function runs without error and writes a non-empty file; does not
inspect pixel content."""

import pandas as pd
from figures import (
    boxplot_by_rank,
    probability_boxplot_by_rank,
    ranked_summary_table,
    scatter_proteome_vs_adhesion,
    scatter_proteome_vs_adhesion_by_phylum,
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


def _many_phyla_df(n_phyla=10):
    rows = []
    for p in range(n_phyla):
        # descending species counts so phylum ranking is unambiguous
        for i in range(n_phyla - p):
            rows.append(
                {
                    "phylum": f"Phylum{p}",
                    "adhesion_fraction": 0.1 + 0.01 * i,
                    "adhesion_count": 10 + i,
                    "total_proteins": 1000 + 50 * i,
                    "mean_adhesion_prob": 0.6 + 0.01 * i,
                }
            )
    return pd.DataFrame(rows)


def test_scatter_by_phylum_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "facet.png"
    scatter_proteome_vs_adhesion_by_phylum(df, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_by_phylum_handles_more_phyla_than_grid(tmp_path):
    """10 phyla into a 3x3 (9-panel) grid should truncate gracefully,
    not raise — only the 9 largest-by-count phyla get a panel."""
    df = _many_phyla_df(n_phyla=10)
    out = tmp_path / "facet_overflow.png"
    scatter_proteome_vs_adhesion_by_phylum(df, out, ncols=3, nrows=3)
    assert out.exists()
    assert out.stat().st_size > 0


def test_phylum_color_map_caps_at_three_highlighted_slots():
    from figures import _CATEGORICAL_ORDER, _OTHER_COLOR, _phylum_color_map

    counts = pd.Series({"A": 100, "B": 80, "C": 60, "D": 40, "E": 20})
    color_map = _phylum_color_map(counts)

    highlighted = {"A", "B", "C"}
    other = {"D", "E"}
    for p in highlighted:
        assert color_map[p] in _CATEGORICAL_ORDER[:3]
    for p in other:
        assert color_map[p] == _OTHER_COLOR
    # highlighted phyla assigned distinct colors, in count-descending order
    assert color_map["A"] != color_map["B"] != color_map["C"]
