#!/usr/bin/env python
"""Merge sequence properties + domain flags into one master table, and
run per-clade aggregation/statistics on adhesion vs. background proteins."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats_ext import (  # noqa: E402
    aggregate_by_clade,
    kruskal_wallis_by_rank,
    mannwhitney_within_clade,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
RANKS = ["phylum", "class", "order"]
CONTINUOUS_FEATURES = [
    "length",
    "pct_ser_thr_pro",
    "aromaticity",
    "mean_hydrophobicity",
    "net_charge_ph7",
]
DOMAIN_FLAGS = ["has_pfam", "has_cazy", "has_merops", "has_signal_peptide", "has_tm_helix"]


def main() -> None:
    props = pd.read_csv(TABLES_DIR / "protein_sequence_properties.csv")
    flags = pd.read_csv(TABLES_DIR / "protein_domain_flags.csv")
    master = props.merge(flags, on="protein_id", how="inner")
    master.to_csv(TABLES_DIR / "adhesion_protein_properties.csv", index=False)
    print(f"Master table: {len(master)} proteins")

    adhesion = master[master["group"] == "adhesion"]
    background = master[master["group"] == "background"]

    for rank in RANKS:
        for feature in CONTINUOUS_FEATURES:
            aggregate_by_clade(adhesion, rank, feature).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_adhesion.csv", index=False
            )
            aggregate_by_clade(background, rank, feature).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_background.csv", index=False
            )
            omnibus_adhesion = kruskal_wallis_by_rank(adhesion, rank, feature)
            omnibus_background = kruskal_wallis_by_rank(background, rank, feature)
            pd.DataFrame([omnibus_adhesion]).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_omnibus_adhesion.csv", index=False
            )
            pd.DataFrame([omnibus_background]).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_omnibus_background.csv", index=False
            )
            mw = mannwhitney_within_clade(adhesion, background, rank, feature)
            mw.to_csv(TABLES_DIR / f"clade_{rank}_{feature}_mannwhitney.csv", index=False)

        rows = []
        for group_name, gdf in [("adhesion", adhesion), ("background", background)]:
            sub = gdf[gdf[rank].notna()]
            for clade, cdf in sub.groupby(rank):
                row = {rank: clade, "group": group_name, "n": len(cdf)}
                for flag in DOMAIN_FLAGS:
                    row[flag] = cdf[flag].mean()
                rows.append(row)
        pd.DataFrame(rows).to_csv(TABLES_DIR / f"clade_{rank}_domain_fractions.csv", index=False)

    print("Aggregation complete.")


if __name__ == "__main__":
    main()
