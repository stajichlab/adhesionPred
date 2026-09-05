#!/usr/bin/env python
"""Extract classifier-space (esm2_t12_35M_UR50D, layer 6) embeddings for
every adhesion-predicted protein. Run as a GPU SLURM job -- see
run_extract_esm2.sbatch. Runs entirely in the project's main Python 3.9
environment (fair-esm) -- do not activate .venv_esmc for this script,
that environment is reserved for the separate ESM-C extraction task."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from embed import (  # noqa: E402
    extract_esm2_classifier_embeddings,
    read_protein_universe_adhesion_ids,
    save_embeddings_chunk,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "adhesion_properties"))
from db import connect, fetch_lengths  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
UNIVERSE_CSV = REPO_ROOT / "analysis" / "adhesion_properties" / "tables" / "protein_universe.csv"
OUT_DIR = Path(__file__).resolve().parent / "tables" / "esm2_classifier"
CHUNK_SIZE = 5_000


def main() -> None:
    protein_ids = read_protein_universe_adhesion_ids(UNIVERSE_CSV)
    print(f"Extracting classifier-space embeddings for {len(protein_ids)} proteins")

    con = connect()
    total_done = 0
    for start in range(0, len(protein_ids), CHUNK_SIZE):
        chunk_idx = start // CHUNK_SIZE
        existing_ids = OUT_DIR / f"chunk_{chunk_idx:05d}_ids.txt"
        existing_emb = OUT_DIR / f"chunk_{chunk_idx:05d}_embeddings.npy"
        if existing_ids.exists() and existing_emb.exists():
            n_existing = sum(1 for _ in existing_ids.open())
            total_done += n_existing
            print(f"  chunk {chunk_idx}: already done, skipping (total {total_done})")
            continue
        chunk_ids = protein_ids[start : start + CHUNK_SIZE]
        seq_data = fetch_lengths(con, chunk_ids)
        sequences = [
            {"id": row["protein_id"], "sequence": row["peptide"] or ""}
            for _, row in seq_data.iterrows()
            if row["peptide"]
        ]
        embeddings, ids = extract_esm2_classifier_embeddings(sequences)
        save_embeddings_chunk(embeddings, ids, OUT_DIR, start // CHUNK_SIZE)
        total_done += len(ids)
        print(f"  chunk {start // CHUNK_SIZE}: {len(ids)} embedded (total {total_done})")
    con.close()

    n_missing = len(protein_ids) - total_done
    print(f"Done: {total_done} embedded, {n_missing} missing/skipped")


if __name__ == "__main__":
    main()
