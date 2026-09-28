"""Shared embedding-extraction helpers for the ESM2-classifier embedding
space, plus chunked save/load bookkeeping used across the embedding-space
clustering pipeline.

Reuses src/adhesion_predict/embeddings.py's exact classifier recipe
(esm2_t12_35M_UR50D, layer 6, mean-pooled) rather than reimplementing it,
so "the classifier's own embedding space" is genuinely the same
representation the classifier was trained/predicts on.

SCOPE NOTE: this module intentionally contains ONLY the ESM2-classifier-space
extraction and shared bookkeeping functions. ESM Cambrian (ESM C) requires a
separate Python >=3.10 virtual environment (`.venv_esmc/`) and a different
package meaning for `import esm` (EvolutionaryScale's ESM-C SDK, vs.
`fair-esm` used here) -- mixing both in one file risks a confusing
AttributeError if this module is ever imported in the wrong interpreter.
ESM-C extraction lives entirely in its own self-contained CLI script
(a separate task), not here.

IMPORT NOTE: `adhesion_predict.embeddings` (which pulls in `fair-esm`) is
imported lazily, inside extract_esm2_classifier_embeddings, rather than at
module level. This lets the ESM-C extraction script import this module's
version-agnostic bookkeeping functions (read_protein_universe_adhesion_ids,
save_embeddings_chunk, load_all_embedding_chunks) under `.venv_esmc`
(Python 3.11), where `adhesion_predict` is not installed and `import esm`
resolves to the unrelated ESM-C SDK package -- without ever needing
adhesion_predict.embeddings to actually load there.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_protein_universe_adhesion_ids(universe_csv_path: Path) -> list[str]:
    """Return every protein_id in the adhesion group of protein_universe.csv."""
    df = pd.read_csv(universe_csv_path)
    return df.loc[df["group"] == "adhesion", "protein_id"].tolist()


def extract_esm2_classifier_embeddings(sequences: list) -> tuple[np.ndarray, list[str]]:
    """Extract embeddings using the exact recipe adhesion_predict's
    classifier was trained/predicts on: esm2_t12_35M_UR50D, layer 6,
    mean-pooled. `sequences` is a list of {"id", "sequence"} dicts,
    matching get_esm_embeddings' expected input shape."""
    _repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(_repo_root / "src"))
    from adhesion_predict.embeddings import get_esm_embeddings

    return get_esm_embeddings(sequences, model_name="esm2_t12_35M_UR50D")


def save_embeddings_chunk(
    embeddings: np.ndarray, ids: list[str], out_dir: Path, chunk_idx: int
) -> None:
    """Save one chunk of embeddings + ids to disk (never accumulate all
    chunks in memory at once -- mirrors the memory-safety lesson from the
    adhesion_properties project's Task 7 OOM incident)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"chunk_{chunk_idx:05d}_embeddings.npy", embeddings)
    (out_dir / f"chunk_{chunk_idx:05d}_ids.txt").write_text("\n".join(ids) + "\n")


def load_all_embedding_chunks(out_dir: Path) -> tuple[np.ndarray, list[str]]:
    """Load and concatenate all chunks written by save_embeddings_chunk,
    in chunk-index order."""
    chunk_files = sorted(out_dir.glob("chunk_*_embeddings.npy"))
    all_embeddings = [np.load(f) for f in chunk_files]
    all_ids: list[str] = []
    for f in chunk_files:
        ids_path = out_dir / f.name.replace("_embeddings.npy", "_ids.txt")
        all_ids.extend(ids_path.read_text().strip().split("\n"))
    return np.vstack(all_embeddings), all_ids
