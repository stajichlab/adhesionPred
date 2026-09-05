#!/usr/bin/env python
"""Generate all figures for the kingdom-wide adhesion survey report."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures import (  # noqa: E402
    boxplot_by_rank,
    probability_boxplot_by_rank,
    scatter_proteome_vs_adhesion,
    scatter_proteome_vs_adhesion_by_phylum,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(TABLES_DIR / "species_adhesion_summary.csv")

    boxplot_by_rank(
        df,
        "phylum",
        "adhesion_fraction",
        FIGURES_DIR / "box_fraction_by_phylum.png",
        title="Adhesion fraction by phylum",
    )
    boxplot_by_rank(
        df,
        "order",
        "adhesion_fraction",
        FIGURES_DIR / "box_fraction_by_order.png",
        top_n=20,
        title="Adhesion fraction by order (top 20 by species count)",
    )
    scatter_proteome_vs_adhesion(df, FIGURES_DIR / "scatter_proteome_vs_adhesion.png")
    scatter_proteome_vs_adhesion_by_phylum(
        df, FIGURES_DIR / "scatter_proteome_vs_adhesion_by_phylum.png"
    )
    probability_boxplot_by_rank(df, "phylum", FIGURES_DIR / "box_prob_by_phylum.png")
    probability_boxplot_by_rank(
        df,
        "order",
        FIGURES_DIR / "box_prob_by_order.png",
        top_n=20,
    )
    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
