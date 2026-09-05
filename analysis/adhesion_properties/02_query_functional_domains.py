#!/usr/bin/env python
"""Query functional domain annotations for the adhesion/background protein
universe, writing per-protein domain flags and top-domain enrichment tables."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import connect  # noqa: E402
from domains import compute_domain_flags, top_domain_table  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"

DOMAIN_TABLES = [
    ("pfam", "pfam_id", "top_domains_pfam.csv"),
    ("cazy_overview", "cazyme_fam", "top_domains_cazy.csv"),
    ("merops", "merops_id", "top_domains_merops.csv"),
]


def main() -> None:
    universe = pd.read_csv(TABLES_DIR / "protein_universe.csv")
    adhesion_ids = universe.loc[universe["group"] == "adhesion", "protein_id"].tolist()
    background_ids = universe.loc[universe["group"] == "background", "protein_id"].tolist()

    con = connect()

    flags = compute_domain_flags(con, universe["protein_id"].tolist())
    flags.to_csv(TABLES_DIR / "protein_domain_flags.csv", index=False)
    print(f"Domain flags computed for {len(flags)} proteins")

    for table, id_field, out_name in DOMAIN_TABLES:
        top = top_domain_table(con, table, id_field, adhesion_ids, background_ids)
        top.to_csv(TABLES_DIR / out_name, index=False)
        print(f"{out_name}: {len(top)} rows")

    con.close()


if __name__ == "__main__":
    main()
