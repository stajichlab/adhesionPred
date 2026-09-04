#!/usr/bin/env python
"""Run taxonomic-rank statistics (naive, genus-averaged, mixed model) on the
master species table produced by 01_build_species_table.py."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats import (  # noqa: E402
    aggregate_summary,
    dunn_posthoc,
    genus_average,
    kruskal_wallis_by_rank,
    mixedlm_by_rank,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
RANKS = ["phylum", "class", "order"]


def _run_rank(df: pd.DataFrame, rank: str, suffix: str) -> None:
    summary = aggregate_summary(df, rank)
    omnibus = kruskal_wallis_by_rank(df, rank, "adhesion_fraction")
    summary.to_csv(TABLES_DIR / f"stats_by_{rank}{suffix}.csv", index=False)
    pd.DataFrame([omnibus]).to_csv(TABLES_DIR / f"stats_by_{rank}{suffix}_omnibus.csv", index=False)
    if omnibus["p_value"] == omnibus["p_value"] and omnibus["p_value"] < 0.05:  # not NaN
        posthoc = dunn_posthoc(df, rank, "adhesion_fraction")
        posthoc.to_csv(TABLES_DIR / f"stats_by_{rank}{suffix}_posthoc.csv")
    print(
        f"{rank}{suffix}: n_groups={omnibus['n_groups']} p={omnibus['p_value']:.4g} "
        f"eps2={omnibus['epsilon_squared']:.4g}"
    )


def main() -> None:
    df = pd.read_csv(TABLES_DIR / "species_adhesion_summary.csv")
    genus_df = genus_average(df)

    for rank in RANKS:
        _run_rank(df, rank, "")
        _run_rank(genus_df, rank, "_genus_avg")

        mixed = mixedlm_by_rank(df, rank)
        pd.DataFrame([{k: v for k, v in mixed.items() if k != "summary"}]).to_csv(
            TABLES_DIR / f"mixedlm_by_{rank}.csv", index=False
        )
        (TABLES_DIR / f"mixedlm_by_{rank}_summary.txt").write_text(
            mixed.get("summary", "did not converge")
        )
        print(f"{rank} mixedlm: converged={mixed['converged']}")


if __name__ == "__main__":
    main()
