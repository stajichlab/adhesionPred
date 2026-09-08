#!/usr/bin/env python
"""Join cluster labels against property/domain data, summarize per
cluster, compare the two clusterings, and check AA1/laccase (CAZy)
concentration by cluster."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate import (  # noqa: E402
    check_feature_concentration,
    cluster_contingency_table,
    summarize_clusters,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
ADHESION_PROPERTIES_TABLES = Path(__file__).resolve().parents[1] / "adhesion_properties" / "tables"


def _build_master(cluster_labels_csv: Path) -> pd.DataFrame:
    props = pd.read_csv(ADHESION_PROPERTIES_TABLES / "protein_sequence_properties.csv")
    flags = pd.read_csv(ADHESION_PROPERTIES_TABLES / "protein_domain_flags.csv")
    clusters = pd.read_csv(cluster_labels_csv)
    merged = props.merge(flags, on="protein_id", how="inner").merge(
        clusters, on="protein_id", how="inner"
    )
    return merged[merged["group"] == "adhesion"] if "group" in merged.columns else merged


def main() -> None:
    esm2_master = _build_master(TABLES_DIR / "cluster_labels_esm2_classifier.csv")
    esmc_master = _build_master(TABLES_DIR / "cluster_labels_esmc300m.csv")

    summarize_clusters(esm2_master, "cluster_label").to_csv(
        TABLES_DIR / "cluster_summary_esm2_classifier.csv", index=False
    )
    summarize_clusters(esmc_master, "cluster_label").to_csv(
        TABLES_DIR / "cluster_summary_esmc300m.csv", index=False
    )
    print("Wrote per-cluster summaries for both embedding spaces")

    merged_both = esm2_master[["protein_id", "cluster_label"]].merge(
        esmc_master[["protein_id", "cluster_label"]], on="protein_id", suffixes=("_esm2", "_esmc")
    )
    contingency = cluster_contingency_table(
        merged_both["cluster_label_esm2"], merged_both["cluster_label_esmc"]
    )
    contingency.to_csv(TABLES_DIR / "cluster_contingency_table.csv")
    print(f"Contingency table: {contingency.shape}")

    # AA1/laccase concentration check: does this protein carry ANY
    # AA1-family CAZy hit? protein_domain_flags.csv only has a broad
    # boolean has_cazy (any CAZy hit at all), not a specific per-protein
    # AA1-family flag. A precise per-protein AA1-only flag would require
    # re-querying cazy_overview.cazyme_fam for this protein set, filtered
    # to AA1-family codes -- out of scope for this task. has_cazy is used
    # here as an approximate proxy and this limitation is stated plainly
    # in the task report.
    for name, master in [("esm2_classifier", esm2_master), ("esmc300m", esmc_master)]:
        result = check_feature_concentration(master, "cluster_label", "has_cazy")
        result.to_csv(TABLES_DIR / f"aa1_concentration_by_cluster_{name}.csv", index=False)

    print("Done.")


if __name__ == "__main__":
    main()
