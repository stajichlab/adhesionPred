#!/usr/bin/env python
"""Generate all figures for the adhesion-protein properties report."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures_ext import (  # noqa: E402
    domain_fraction_barplot,
    paired_boxplot_by_rank,
    scatter_feature_vs_probability,
    scatter_species_mean_length_vs_fraction,
    top_domain_enrichment_barplot,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"
KINGDOM_SUMMARY_CSV = (
    Path(__file__).resolve().parents[1]
    / "kingdom_survey"
    / "tables"
    / "species_adhesion_summary.csv"
)

DOMAIN_FLAGS = ["has_pfam", "has_cazy", "has_merops", "has_signal_peptide", "has_tm_helix"]


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(TABLES_DIR / "adhesion_protein_properties.csv")
    adhesion = master[master["group"] == "adhesion"]

    paired_boxplot_by_rank(master, "phylum", "length", FIGURES_DIR / "box_length_by_phylum.png")
    paired_boxplot_by_rank(
        master, "phylum", "pct_ser_thr_pro", FIGURES_DIR / "box_pct_ser_thr_pro_by_phylum.png"
    )

    fractions = pd.read_csv(TABLES_DIR / "clade_phylum_domain_fractions.csv")
    for flag in DOMAIN_FLAGS:
        domain_fraction_barplot(
            fractions, "phylum", flag, FIGURES_DIR / f"bar_{flag}_by_phylum.png"
        )

    for domain_type in ["pfam", "cazy", "merops"]:
        top = pd.read_csv(TABLES_DIR / f"top_domains_{domain_type}.csv")
        if top.empty:
            print(f"Skipping top_domains_{domain_type}.png — no enriched domains found")
            continue
        top_domain_enrichment_barplot(
            top,
            FIGURES_DIR / f"top_domains_{domain_type}.png",
            title=f"Top enriched {domain_type.upper()} domains among adhesion proteins",
        )

    scatter_feature_vs_probability(
        adhesion, "length", FIGURES_DIR / "scatter_length_vs_probability.png"
    )
    scatter_feature_vs_probability(
        adhesion, "pct_ser_thr_pro", FIGURES_DIR / "scatter_pct_ser_thr_pro_vs_probability.png"
    )

    kingdom_summary = pd.read_csv(KINGDOM_SUMMARY_CSV)
    scatter_species_mean_length_vs_fraction(
        adhesion, kingdom_summary, FIGURES_DIR / "scatter_species_mean_length_vs_fraction.png"
    )

    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
