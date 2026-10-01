"""Chunk files of the Phase B embeddings: save, verify, load. Needs numpy.

Each chunk is two files in <publish>/<model>/: `<chunk_id>.npy` (float32, one row per chunk
member, in chunk order) and `<chunk_id>.json` (the done marker). The JSON is written last, so
a chunk without JSON is not done. Both files are first written to a temporary name and then
moved with os.replace, so a killed job never leaves a half-written chunk under its final name.
"""

import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
from embed_constants import MODEL_DIM, REPR_LAYER

DTYPE = np.float32


def array_sha256(arr: np.ndarray) -> str:
    """SHA-256 of the raw array bytes (C order). Shape and dtype are stored beside it."""
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def chunk_paths(publish: Path, model: str, chunk_id: str) -> tuple[Path, Path]:
    base = Path(publish) / model
    return base / f"{chunk_id}.npy", base / f"{chunk_id}.json"


def _atomic_copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".tmp.{dest.name}")
    try:
        shutil.copyfile(src, tmp)
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)


def save_chunk(publish: Path, scratch: Path, model: str, chunk_id: str, arr, meta: dict) -> dict:
    """Write the chunk to `scratch`, then copy it to `publish`. Return the stored metadata."""
    arr = np.ascontiguousarray(arr, dtype=DTYPE)
    if arr.ndim != 2:
        raise ValueError(f"{chunk_id}: expected a 2-D array, got shape {arr.shape}")
    if not np.isfinite(arr).all():
        raise ValueError(f"{chunk_id}: array has NaN or inf")
    meta = {
        **meta,
        "chunk_id": chunk_id,
        "model": model,
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "array_sha256": array_sha256(arr),
    }
    s_npy, s_json = chunk_paths(scratch, model, chunk_id)
    s_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(s_npy, arr)
    s_json.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    p_npy, p_json = chunk_paths(publish, model, chunk_id)
    _atomic_copy(s_npy, p_npy)
    _atomic_copy(s_json, p_json)
    return meta


def load_chunk(publish: Path, model: str, chunk_id: str):
    """Return (array, metadata) of a finished chunk. Raise ValueError if it fails a check."""
    npy, js = chunk_paths(publish, model, chunk_id)
    meta = json.loads(js.read_text())
    arr = np.load(npy, allow_pickle=False)
    if list(arr.shape) != meta["shape"] or str(arr.dtype) != meta["dtype"]:
        raise ValueError(f"{model}/{chunk_id}: shape or dtype differs from its JSON")
    if array_sha256(arr) != meta["array_sha256"]:
        raise ValueError(f"{model}/{chunk_id}: array SHA-256 differs from its JSON")
    return arr, meta


def check_chunk_meta(model: str, arr: np.ndarray, meta: dict) -> str | None:
    """Reason why a loaded chunk is not an ESM-2 layer-6 chunk of `model`, else None."""
    if meta.get("model") != model:
        return f"the file records model {meta.get('model')!r}, not {model!r}"
    if meta.get("repr_layer") != REPR_LAYER:
        return f"the file records repr_layer {meta.get('repr_layer')!r}, not {REPR_LAYER}"
    dim = MODEL_DIM.get(model)
    if dim is None:
        return f"model {model!r} has no known embedding width (embed_constants.MODEL_DIM)"
    if arr.ndim != 2 or arr.shape[1] != dim:
        return f"the array has shape {list(arr.shape)}, the width of {model} is {dim}"
    return None


def chunk_is_done(
    publish: Path,
    model: str,
    chunk_id: str,
    members_sha256: str,
    n: int,
    window: str | None = None,
) -> bool:
    """True if the chunk exists, matches the plan (members, row count, window), the model
    (name, layer, width) and its hash."""
    npy, js = chunk_paths(publish, model, chunk_id)
    if not (npy.exists() and js.exists()):
        return False
    try:
        arr, meta = load_chunk(publish, model, chunk_id)
    except (OSError, ValueError, KeyError):
        return False
    if check_chunk_meta(model, arr, meta) is not None:
        return False
    if window is not None and meta.get("window") != window:
        return False
    return meta.get("members_sha256") == members_sha256 and arr.shape[0] == n
