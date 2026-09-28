"""Pure bookkeeping helpers for 02_extract_embeddings_esmc300m.py, split
out so they can be unit-tested (and imported by tests) without needing
`torch`/`esm` or a real GPU -- this module has zero model-related
imports and works unchanged under both the project's main Python 3.9
environment and the ESM-C-specific .venv_esmc Python 3.11 environment.
"""

from pathlib import Path

import numpy as np


def chunk_is_complete(ids_path: Path, emb_path: Path) -> tuple:
    """Check whether a previously-written chunk pair is complete and
    trustworthy, returning (is_complete, n_ids).

    Hardened over Task 3's resume-skip logic (which only checked file
    *existence*): also verifies the ids.txt line count matches the
    embeddings.npy row count before trusting the chunk as done. A kill
    mid-write of ids.txt itself (leaving it truncated-but-present, with
    the .npy already written in full, or vice versa) would pass an
    existence-only check but fail this shape/line-count check -- this
    closes exactly the gap Task 3's reviewer flagged as a latent risk.
    """
    if not (ids_path.exists() and emb_path.exists()):
        return False, 0
    n_ids = sum(1 for _ in ids_path.open())
    # mmap_mode='r' avoids loading the full array into memory just to
    # check its shape.
    try:
        arr = np.load(emb_path, mmap_mode="r")
    except (OSError, ValueError):
        # Truncated/corrupt .npy file -- treat as incomplete.
        return False, 0
    n_rows = arr.shape[0]
    if n_ids != n_rows:
        return False, 0
    return True, n_ids
