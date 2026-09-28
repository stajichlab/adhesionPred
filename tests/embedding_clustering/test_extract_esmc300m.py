"""Tests for the resume-skip bookkeeping in
analysis/embedding_clustering/02_extract_embeddings_esmc300m.py --
exercises only chunk_is_complete's pure file/shape-consistency logic,
never real ESM-C model inference (no GPU/network needed for these
tests, and this module cannot even import cleanly outside the
.venv_esmc Python >=3.10 environment -- see the module's own docstring
-- so these tests import only the one pure-Python helper function via
a private module-loading trick that avoids importing `torch`/`esm` at
collection time is unnecessary here since torch is already a project
dependency; only `esm` (ESM-C SDK) is venv-specific, and
chunk_is_complete itself never touches it)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analysis" / "embedding_clustering"))

from extract_esmc300m_helpers import chunk_is_complete  # noqa: E402


def _write_chunk(tmp_path, ids, embeddings):
    ids_path = tmp_path / "chunk_00000_ids.txt"
    emb_path = tmp_path / "chunk_00000_embeddings.npy"
    ids_path.write_text("\n".join(ids) + "\n")
    np.save(emb_path, embeddings)
    return ids_path, emb_path


def test_chunk_is_complete_true_when_consistent(tmp_path):
    ids_path, emb_path = _write_chunk(tmp_path, ["A", "B"], np.zeros((2, 4), dtype=np.float32))
    is_complete, n = chunk_is_complete(ids_path, emb_path)
    assert is_complete
    assert n == 2


def test_chunk_is_complete_false_when_missing(tmp_path):
    ids_path = tmp_path / "chunk_00000_ids.txt"
    emb_path = tmp_path / "chunk_00000_embeddings.npy"
    is_complete, n = chunk_is_complete(ids_path, emb_path)
    assert not is_complete
    assert n == 0


def test_chunk_is_complete_false_when_truncated_ids(tmp_path):
    # embeddings.npy claims 2 rows but ids.txt has only 1 line -- a
    # truncated-but-present ids file, the exact gap Task 3's reviewer
    # flagged as a latent risk in file-existence-only checking.
    ids_path, emb_path = _write_chunk(tmp_path, ["A"], np.zeros((2, 4), dtype=np.float32))
    is_complete, n = chunk_is_complete(ids_path, emb_path)
    assert not is_complete
    assert n == 0


def test_chunk_is_complete_false_when_corrupt_npy(tmp_path):
    ids_path = tmp_path / "chunk_00000_ids.txt"
    emb_path = tmp_path / "chunk_00000_embeddings.npy"
    ids_path.write_text("A\nB\n")
    emb_path.write_bytes(b"not a real npy file")
    is_complete, n = chunk_is_complete(ids_path, emb_path)
    assert not is_complete
    assert n == 0
