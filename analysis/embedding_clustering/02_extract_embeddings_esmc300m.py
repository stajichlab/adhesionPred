#!/usr/bin/env python
"""Extract ESM Cambrian 300M (esmc_300m) embeddings for every
adhesion-predicted protein. Must be run with .venv_esmc/bin/python (a
separate Python 3.11 environment) -- see run_extract_esmc300m.sbatch.

Do NOT import embed.py's extract_esm2_classifier_embeddings from this
script: that function pulls in fair-esm, which requires the project's
main Python 3.9 environment and is a different meaning of `import esm`
than the ESM-C SDK package installed in .venv_esmc (see Task 1's
report and embed.py's own SCOPE NOTE docstring for why these two are
deliberately kept apart). This script only reuses the version-agnostic
bookkeeping helpers from embed.py (read_protein_universe_adhesion_ids,
save_embeddings_chunk) and db.py (connect, fetch_lengths), and contains
its own inline ESM-C model-loading + special-token-aware mean-pooling
code, based exactly on Task 1's verified snippet:

    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig
    model = ESMC.from_pretrained('esmc_300m').eval()
    tok = model.tokenizer
    cls_id, eos_id, pad_id = tok.cls_token_id, tok.eos_token_id, tok.pad_token_id
    protein_tensor = model.encode(ESMProtein(sequence=seq))
    out = model.logits(protein_tensor, LogitsConfig(sequence=True, return_embeddings=True))
    emb = out.embeddings[0]           # (L, 960)
    tokens = protein_tensor.sequence   # (L,)
    mask = (tokens != cls_id) & (tokens != eos_id) & (tokens != pad_id)
    pooled = emb[mask].mean(dim=0)     # (960,)
"""

import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from embed import (  # noqa: E402
    read_protein_universe_adhesion_ids,
    save_embeddings_chunk,
)
from extract_esmc300m_helpers import chunk_is_complete  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "adhesion_properties"))
from db import connect, fetch_lengths  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
UNIVERSE_CSV = REPO_ROOT / "analysis" / "adhesion_properties" / "tables" / "protein_universe.csv"
OUT_DIR = Path(__file__).resolve().parent / "tables" / "esmc300m"
CHUNK_SIZE = 5_000
EMBED_DIM = 960


def load_esmc_model():
    """Load ESM Cambrian 300M and return (model, cls_id, eos_id, pad_id).

    Uses the exact approach verified in Task 1's feasibility check
    (approach A: EvolutionaryScale's `esm` PyPI package, only installable
    under Python >=3.10 -- this script must run under .venv_esmc)."""
    from esm.models.esmc import ESMC

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Loading ESM C 300M on {device}...")
    model = ESMC.from_pretrained("esmc_300m").eval().to(device)
    tok = model.tokenizer
    return model, tok.cls_token_id, tok.eos_token_id, tok.pad_token_id, device


def extract_esmc300m_embeddings_batch(
    sequences: list, model, cls_id: int, eos_id: int, pad_id: int, device: str
) -> tuple:
    """Extract ESM C 300M embeddings for a list of {"id", "sequence"}
    dicts, mean-pooled over residue tokens only (special tokens masked
    out per Task 1's verified pooling logic). One sequence at a time
    through model.encode/model.logits (the ESM-C SDK's own API shape,
    which does not expose a batched-tensor path the way fair-esm does)."""
    from esm.sdk.api import ESMProtein, LogitsConfig

    embeddings = np.zeros((len(sequences), EMBED_DIM), dtype=np.float32)
    ids = []
    logits_config = LogitsConfig(sequence=True, return_embeddings=True)
    with torch.no_grad():
        for i, item in enumerate(sequences):
            protein = ESMProtein(sequence=item["sequence"])
            protein_tensor = model.encode(protein)
            out = model.logits(protein_tensor, logits_config)
            emb = out.embeddings[0]  # (L, 960)
            tokens = protein_tensor.sequence  # (L,)
            mask = (tokens != cls_id) & (tokens != eos_id) & (tokens != pad_id)
            pooled = emb[mask].mean(dim=0)  # (960,)
            embeddings[i] = pooled.detach().cpu().numpy().astype(np.float32)
            ids.append(item["id"])
    return embeddings, ids


def main() -> None:
    protein_ids = read_protein_universe_adhesion_ids(UNIVERSE_CSV)
    print(f"Extracting ESM C 300M embeddings for {len(protein_ids)} proteins")

    model, cls_id, eos_id, pad_id, device = load_esmc_model()

    con = connect()
    total_done = 0
    for start in range(0, len(protein_ids), CHUNK_SIZE):
        chunk_idx = start // CHUNK_SIZE
        ids_path = OUT_DIR / f"chunk_{chunk_idx:05d}_ids.txt"
        emb_path = OUT_DIR / f"chunk_{chunk_idx:05d}_embeddings.npy"

        is_complete, n_existing = chunk_is_complete(ids_path, emb_path)
        if is_complete:
            total_done += n_existing
            print(f"  chunk {chunk_idx}: already done, skipping (total {total_done})")
            continue
        if ids_path.exists() or emb_path.exists():
            print(f"  chunk {chunk_idx}: existing files incomplete/inconsistent, recomputing")

        chunk_ids = protein_ids[start : start + CHUNK_SIZE]
        seq_data = fetch_lengths(con, chunk_ids)
        sequences = [
            {"id": row["protein_id"], "sequence": row["peptide"] or ""}
            for _, row in seq_data.iterrows()
            if row["peptide"]
        ]
        embeddings, ids = extract_esmc300m_embeddings_batch(
            sequences, model, cls_id, eos_id, pad_id, device
        )
        save_embeddings_chunk(embeddings, ids, OUT_DIR, chunk_idx)
        total_done += len(ids)
        print(f"  chunk {chunk_idx}: {len(ids)} embedded (total {total_done})")
    con.close()

    n_missing = len(protein_ids) - total_done
    print(f"Done: {total_done} embedded, {n_missing} missing/skipped")


if __name__ == "__main__":
    main()
