#!/usr/bin/env python
"""Compute length and amino-acid-composition properties for every protein
in the adhesion/background universe."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import connect, fetch_lengths  # noqa: E402
from properties import sequence_properties  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"

# The full universe is ~1.5M proteins; pulling all peptide sequences into
# memory at once (fetchdf() text + per-row property dicts) exceeded an 8GB
# job memory cap and was OOM-killed by the cgroup. Processing in chunks and
# writing incrementally keeps peak memory to one chunk's worth of sequence
# text instead of the full universe's.
CHUNK_SIZE = 50_000


def main() -> None:
    universe = pd.read_csv(TABLES_DIR / "protein_universe.csv")
    out_path = TABLES_DIR / "protein_sequence_properties.csv"
    ids = universe["protein_id"].tolist()

    con = connect()
    total_matched = 0
    first_chunk = True
    for start in range(0, len(ids), CHUNK_SIZE):
        chunk_ids = ids[start : start + CHUNK_SIZE]
        seq_data = fetch_lengths(con, chunk_ids)

        prop_rows = []
        for _, row in seq_data.iterrows():
            props = sequence_properties(row["peptide"] or "")
            prop_rows.append({"protein_id": row["protein_id"], "length": row["length"], **props})

        props_df = pd.DataFrame(prop_rows)
        universe_chunk = universe.iloc[start : start + CHUNK_SIZE]
        merged = universe_chunk.merge(props_df, on="protein_id", how="inner")
        merged.to_csv(out_path, mode="w" if first_chunk else "a", header=first_chunk, index=False)
        first_chunk = False
        total_matched += len(merged)
    con.close()

    n_missing = len(universe) - total_matched
    print(
        f"Computed properties for {total_matched} proteins ({n_missing} not found in gene_proteins)"
    )


if __name__ == "__main__":
    main()
