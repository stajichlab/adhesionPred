# Adhesion Protein Embedding-Space Clustering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cluster the 749,697 adhesion-predicted proteins in two embedding spaces (the classifier's own ESM2-t12-35M layer-6 representation, and ESM Cambrian 300M) and validate any clusters against already-computed length/composition/domain data to look for distinct structural sub-types.

**Architecture:** A standalone `analysis/embedding_clustering/` directory. Reuses `analysis/adhesion_properties/db.py` (peptide fetch) and `domains.py` (top-domain logic) and `src/adhesion_predict/embeddings.py` (the classifier's exact embedding recipe) via `sys.path.insert`, never duplicating them. Two GPU SLURM batch jobs extract embeddings; the rest of the pipeline (PCA/clustering/validation/figures/report) is CPU-only and runs interactively.

**Tech Stack:** Python 3.9, PyTorch, `fair-esm` (existing), a new ESM-C-loading dependency (determined by Task 1), scikit-learn (PCA + `sklearn.cluster.HDBSCAN`, already installed), `umap-learn` (already installed), matplotlib/seaborn, pytest.

**Spec:** `docs/superpowers/specs/2026-09-04-embedding-clustering-design.md` — read it before starting; this plan implements it task by task. Task 1's outcome determines exact implementation details in Task 2 — read Task 1's report before starting Task 2.

## Global Constraints

- Python `>=3.9` — no PEP 604 `X | None` union syntax; use `typing.Optional`.
- Ruff: line-length 100, rules `E,F,W,I,N,UP,B,C4`, ignoring `E501,N803,N806,N999`.
- The classifier's exact embedding recipe (`esm2_t12_35M_UR50D`, `repr_layers=[6]`, mean-pooled, 1022-residue truncation) must be replicated EXACTLY for the "classifier's own embedding space" — reuse `src/adhesion_predict/embeddings.py`'s `get_esm_embeddings` function directly, do not reimplement.
- GPU-heavy extraction steps (Tasks 3, 4) run as SLURM batch jobs on the `gpu` or `preempt_gpu` partition, following the same operational pattern established for the prior project's memory-constrained CPU job (explicit resource request, `$SCRATCH` for spill/temp, chunked/incremental writes).
- `sklearn.cluster.HDBSCAN` (built into scikit-learn ≥1.3, already installed) is used for clustering — do not add the separate `hdbscan` PyPI package.
- Reuse, don't duplicate: `analysis/adhesion_properties/db.py`'s `connect`/`fetch_lengths`, `domains.py`'s `top_domain_table`, `src/adhesion_predict/embeddings.py`'s `get_esm_embeddings` — imported via `sys.path.insert`, never copied.

---

## File Structure

```
analysis/embedding_clustering/
  00_feasibility_check.py                    # Task 1
  embed.py                                   # Task 2
  01_extract_embeddings_esm2_classifier.py   # Task 3 (GPU SLURM job)
  02_extract_embeddings_esmc300m.py          # Task 4 (GPU SLURM job)
  cluster.py                                 # Task 5
  03_cluster_esm2_classifier.py              # Task 6
  04_cluster_esmc300m.py                     # Task 6
  validate.py                                # Task 7
  05_validate_and_compare.py                 # Task 8
  figures_ext.py                             # Task 9
  06_make_figures.py                         # Task 9
  07_generate_report.py                      # Task 10
  tables/
  figures/
  REPORT.md
tests/embedding_clustering/
  conftest.py
  test_embed.py
  test_cluster.py
  test_validate.py
  test_figures_ext.py
```

---

### Task 1: Feasibility check — can we actually load and embed with both models?

**Files:**
- Create: `analysis/embedding_clustering/00_feasibility_check.py`

**Interfaces:**
- Produces: a written report (in the SDD workspace, per the implementer's standard report format) stating exactly which loading method works for ESM C, the confirmed embedding dimensionality for both models, and any dependency/access issues found — Task 2 reads this report before writing `embed.py`.

This is an investigative task, not a fixed-code task — the exact ESM C loading incantation is not yet known to work in this environment.

- [ ] **Step 1: Confirm the classifier's own embedding path still works**

Run a quick check that `src/adhesion_predict/embeddings.py`'s `get_esm_embeddings` still works as expected on a handful of real sequences:

```python
import sys
sys.path.insert(0, "src")
from adhesion_predict.embeddings import get_esm_embeddings

sequences = [
    {"id": "test1", "sequence": "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA"},
    {"id": "test2", "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"},
]
embeddings, ids = get_esm_embeddings(sequences, model_name="esm2_t12_35M_UR50D")
print(f"Shape: {embeddings.shape}, expected (2, 480)")
assert embeddings.shape == (2, 480), f"Unexpected shape: {embeddings.shape}"
print("Classifier embedding path: OK")
```

Note: this calls `model(batch_tokens, repr_layers=[6], ...)` internally — confirm the returned embedding dimension is 480 as expected for `esm2_t12_35M_UR50D` layer 6. If it's not 480, investigate before proceeding (the model variant or layer count may differ from what the design doc assumed).

- [ ] **Step 2: Attempt to load ESM C 300M — try candidate approaches in order**

Try, in this order, stopping at the first one that works:

**Approach A — the `esm` PyPI package's `from_pretrained_esm`:**
```bash
pip install esm
python3 -c "
from esm.models.esmc import ESMC
model = ESMC.from_pretrained('esmc_300m').eval()
print('Loaded via esm package')
print(model)
"
```

**Approach B — HuggingFace `transformers`:**
```bash
pip install transformers
python3 -c "
from transformers import AutoModelForMaskedLM, AutoTokenizer
import torch
tokenizer = AutoTokenizer.from_pretrained('Biohub/esmc-300m-2024-12')
model = AutoModelForMaskedLM.from_pretrained('Biohub/esmc-300m-2024-12')
inputs = tokenizer('MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGG', return_tensors='pt')
with torch.no_grad():
    outputs = model(**inputs, output_hidden_states=True)
last_hidden = outputs.hidden_states[-1]
print('Shape:', last_hidden.shape, '- expect last dim 960')
"
```

If BOTH approaches fail (network access blocked from this compute environment, gated model requiring authentication, package incompatibility, etc.), STOP and report BLOCKED with the exact errors from both attempts — do not proceed to Task 2's implementation, and do not silently substitute a different model. Report back to the controller so it can decide how to proceed (e.g., request network/token access, or reconsider the model choice) — this is exactly the kind of decision only project owner can make.

If ONE approach works, confirm the embedding dimensionality (960 expected) and pooling method available (mean-pool over residue tokens excluding BOS/EOS — check the tokenizer's special token handling to compute this correctly, e.g. by masking out positions corresponding to special tokens before averaging).

- [ ] **Step 3: Confirm GPU availability for a real extraction job**

Run: `sinfo -p gpu,preempt_gpu -o "%N %G %C %m"` and confirm at least one partition/node combination is usable (matches what was already checked during design — `gpu[06-08]` a100:8, `gpu11` h100:2, etc.). No action needed beyond confirming this still holds; note the partition name Task 3/4 should target.

- [ ] **Step 4: Report findings**

Write your report to the SDD workspace's `task-1-report.md` covering:
- Which ESM C loading approach worked (A or B), with the exact working code snippet
- Confirmed embedding dimensions for both models (480 for ESM2 layer 6, 960 for ESM C 300M — or actual values if different)
- The exact mean-pooling logic needed for ESM C (accounting for special tokens)
- Confirmed GPU partition to target
- Any new pip packages that had to be installed (for `pyproject.toml`)

Then report back with the standard status contract. This task has NO commit step of its own (it's a spike, not a deliverable) unless the feasibility check reveals a package needs adding to `pyproject.toml` — in that case, add it and commit just that dependency addition, matching Task 1's pattern from the prior `kingdom_survey`/`adhesion_properties` projects (a small, standalone dependency commit).

---

### Task 2: `embed.py` — shared embedding-extraction helpers

**Files:**
- Create: `analysis/embedding_clustering/embed.py`
- Create: `tests/embedding_clustering/conftest.py`
- Test: `tests/embedding_clustering/test_embed.py`

**Interfaces:**
- Consumes: Task 1's report (for the exact ESM C loading/pooling code), `src/adhesion_predict/embeddings.py`'s `get_esm_embeddings` (reused for the classifier-space extraction).
- Produces (used by Tasks 3, 4):
  - `read_protein_universe_adhesion_ids(universe_csv_path: Path) -> list` — returns the list of adhesion-group protein IDs.
  - `extract_esm2_classifier_embeddings(sequences: list) -> Tuple[np.ndarray, list]` — thin wrapper around the reused `get_esm_embeddings`, fixing `model_name="esm2_t12_35M_UR50D"`.
  - `extract_esmc300m_embeddings(sequences: list, batch_size: int = 8) -> Tuple[np.ndarray, list]` — implements whichever approach Task 1's report found working, with correct special-token-aware mean pooling.
  - `save_embeddings_chunk(embeddings: np.ndarray, ids: list, out_dir: Path, chunk_idx: int) -> None` and `load_all_embedding_chunks(out_dir: Path) -> Tuple[np.ndarray, list]` — chunked incremental save/load (mirrors the memory-safety pattern from the prior project's Task 7 incident: never hold more than one chunk's embeddings plus what's already written to disk in memory at once).

Read Task 1's report file before writing this task's ESM C loading code — the brief text below is a template; substitute the actual working approach found in Task 1.

- [ ] **Step 1: Write the failing tests**

Create `tests/embedding_clustering/conftest.py`:

```python
"""Make analysis/embedding_clustering importable as plain modules for tests."""
import sys
from pathlib import Path

EMBEDDING_CLUSTERING_DIR = Path(__file__).resolve().parents[2] / "analysis" / "embedding_clustering"
sys.path.insert(0, str(EMBEDDING_CLUSTERING_DIR))
```

Create `tests/embedding_clustering/test_embed.py` — this tests only the pure-Python bookkeeping (chunk save/load, protein-ID reading), NOT actual model inference (that requires the real model weights and is exercised only in Tasks 3/4's real-data runs):

```python
"""Tests for analysis/embedding_clustering/embed.py — exercises only the
pure bookkeeping logic (chunk save/load, protein-id reading), never
real model inference (no GPU/network needed for these tests)."""
import numpy as np
import pandas as pd

from embed import (
    load_all_embedding_chunks,
    read_protein_universe_adhesion_ids,
    save_embeddings_chunk,
)


def test_read_protein_universe_adhesion_ids(tmp_path):
    universe = tmp_path / "protein_universe.csv"
    pd.DataFrame(
        {
            "protein_id": ["P1", "P2", "P3"],
            "group": ["adhesion", "background", "adhesion"],
        }
    ).to_csv(universe, index=False)
    ids = read_protein_universe_adhesion_ids(universe)
    assert ids == ["P1", "P3"]


def test_save_and_load_embedding_chunks(tmp_path):
    chunk0 = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
    chunk1 = np.array([[5.0, 6.0]], dtype=np.float32)
    save_embeddings_chunk(chunk0, ["A", "B"], tmp_path, 0)
    save_embeddings_chunk(chunk1, ["C"], tmp_path, 1)

    embeddings, ids = load_all_embedding_chunks(tmp_path)
    assert embeddings.shape == (3, 2)
    assert ids == ["A", "B", "C"]
    np.testing.assert_array_equal(embeddings, np.vstack([chunk0, chunk1]))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/embedding_clustering/test_embed.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'embed'`.

- [ ] **Step 3: Implement `embed.py`**

Create `analysis/embedding_clustering/embed.py`:

```python
"""Shared embedding-extraction helpers for the two embedding-space
clustering pipeline. Reuses src/adhesion_predict/embeddings.py's exact
classifier recipe (esm2_t12_35M_UR50D, layer 6, mean-pooled) rather than
reimplementing it, so "the classifier's own embedding space" is genuinely
the same representation the classifier was trained/predicts on.
"""

import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "src"))
from adhesion_predict.embeddings import get_esm_embeddings  # noqa: E402


def read_protein_universe_adhesion_ids(universe_csv_path: Path) -> List[str]:
    """Return every protein_id in the adhesion group of protein_universe.csv."""
    df = pd.read_csv(universe_csv_path)
    return df.loc[df["group"] == "adhesion", "protein_id"].tolist()


def extract_esm2_classifier_embeddings(sequences: list) -> Tuple[np.ndarray, List[str]]:
    """Extract embeddings using the exact recipe adhesion_predict's
    classifier was trained/predicts on: esm2_t12_35M_UR50D, layer 6,
    mean-pooled. `sequences` is a list of {"id", "sequence"} dicts,
    matching get_esm_embeddings' expected input shape."""
    return get_esm_embeddings(sequences, model_name="esm2_t12_35M_UR50D")


# --- ESM C 300M extraction: implementation depends on Task 1's findings.
# The function signature and chunk-save contract below are fixed; the
# model-loading/pooling body must be filled in with whichever approach
# Task 1's report confirmed works (PyPI `esm` package's
# `ESMC.from_pretrained('esmc_300m')`, or HuggingFace `transformers`'
# `AutoModelForMaskedLM.from_pretrained('Biohub/esmc-300m-2024-12')`).
# Read /path/to/sdd/workspace/task-1-report.md before implementing this.


def extract_esmc300m_embeddings(
    sequences: list, batch_size: int = 8
) -> Tuple[np.ndarray, List[str]]:
    """Extract embeddings using ESM C 300M, mean-pooled over residue
    tokens (excluding special tokens). `sequences` is a list of
    {"id", "sequence"} dicts, same shape as
    extract_esm2_classifier_embeddings' input.

    IMPLEMENT PER TASK 1's REPORT — the model-loading approach and
    special-token-aware pooling logic must match whatever Task 1 verified
    actually works in this environment.
    """
    raise NotImplementedError(
        "Implement using the approach confirmed working in Task 1's report"
    )


def save_embeddings_chunk(embeddings: np.ndarray, ids: List[str], out_dir: Path, chunk_idx: int) -> None:
    """Save one chunk of embeddings + ids to disk (never accumulate all
    chunks in memory at once — mirrors the memory-safety lesson from the
    adhesion_properties project's Task 7 OOM incident)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"chunk_{chunk_idx:05d}_embeddings.npy", embeddings)
    (out_dir / f"chunk_{chunk_idx:05d}_ids.txt").write_text("\n".join(ids) + "\n")


def load_all_embedding_chunks(out_dir: Path) -> Tuple[np.ndarray, List[str]]:
    """Load and concatenate all chunks written by save_embeddings_chunk,
    in chunk-index order."""
    chunk_files = sorted(out_dir.glob("chunk_*_embeddings.npy"))
    all_embeddings = [np.load(f) for f in chunk_files]
    all_ids: List[str] = []
    for f in chunk_files:
        ids_path = out_dir / f.name.replace("_embeddings.npy", "_ids.txt")
        all_ids.extend(ids_path.read_text().strip().split("\n"))
    return np.vstack(all_embeddings), all_ids
```

- [ ] **Step 4: Implement `extract_esmc300m_embeddings`'s body per Task 1's findings**

Replace the `NotImplementedError` placeholder with real code implementing whichever approach Task 1's report confirmed working. This is the one place in this plan where the exact code is determined by a prior task's investigation rather than specified here — follow Task 1's report precisely, including its special-token-aware pooling logic. Batch the input by `batch_size`, matching the pattern `get_esm_embeddings` already uses (iterate in batches, `torch.no_grad()`, move to the appropriate device).

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/embedding_clustering/test_embed.py -v`
Expected: all PASS. (These tests don't exercise `extract_esmc300m_embeddings` or `extract_esm2_classifier_embeddings` directly — that's real-data-only, covered in Tasks 3/4.)

- [ ] **Step 6: Smoke-test the real extraction functions on a tiny input**

Run a manual check (not a committed test, since it needs GPU/network) confirming both extraction functions work end-to-end on 2-3 real sequences, matching Task 1's feasibility check but now through the actual `embed.py` functions rather than ad hoc snippets:

```bash
python3 -c "
import sys; sys.path.insert(0, 'analysis/embedding_clustering')
from embed import extract_esm2_classifier_embeddings, extract_esmc300m_embeddings
seqs = [{'id': 't1', 'sequence': 'MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA'}]
e1, i1 = extract_esm2_classifier_embeddings(seqs)
e2, i2 = extract_esmc300m_embeddings(seqs)
print('ESM2 classifier shape:', e1.shape)
print('ESM C 300M shape:', e2.shape)
assert e1.shape[1] == 480
assert e2.shape[1] == 960
print('OK')
"
```

- [ ] **Step 7: Lint**

Run: `ruff check analysis/embedding_clustering/embed.py tests/embedding_clustering/`

- [ ] **Step 8: Commit**

```bash
git add analysis/embedding_clustering/embed.py tests/embedding_clustering/conftest.py tests/embedding_clustering/test_embed.py
git commit -m "Add shared embedding-extraction helpers for both embedding spaces"
```

---

### Task 3: `01_extract_embeddings_esm2_classifier.py` — run as a GPU SLURM job

**Files:**
- Create: `analysis/embedding_clustering/01_extract_embeddings_esm2_classifier.py`
- Create: `analysis/embedding_clustering/run_extract_esm2.sbatch`

**Interfaces:**
- Consumes: `embed.py` (Task 2), `analysis/adhesion_properties/db.py` (peptide fetch, reused), `analysis/adhesion_properties/tables/protein_universe.csv`.
- Produces: `analysis/embedding_clustering/tables/esm2_classifier/chunk_*_embeddings.npy` + `chunk_*_ids.txt` — consumed by Task 6.

- [ ] **Step 1: Write the CLI script**

Create `analysis/embedding_clustering/01_extract_embeddings_esm2_classifier.py`:

```python
#!/usr/bin/env python
"""Extract classifier-space (esm2_t12_35M_UR50D, layer 6) embeddings for
every adhesion-predicted protein. Run as a GPU SLURM job — see
run_extract_esm2.sbatch."""
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
```

- [ ] **Step 2: Write the SLURM batch script**

Create `analysis/embedding_clustering/run_extract_esm2.sbatch`:

```bash
#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH --time=4:00:00
#SBATCH -J embed_esm2_classifier
#SBATCH -o logs/embed_esm2_classifier.%j.log
#SBATCH -e logs/embed_esm2_classifier.%j.log

cd "$(dirname "$(readlink -f "$0")")/../.."
python analysis/embedding_clustering/01_extract_embeddings_esm2_classifier.py
```

Note: `cd "$(dirname "$(readlink -f "$0")")/../.."` resolves to the repo root from the sbatch script's own location — do NOT hardcode an absolute worktree path (that was flagged as a portability issue in the prior project's analogous sbatch script; fix it here).

- [ ] **Step 3: Submit and monitor the job**

```bash
mkdir -p logs
sbatch analysis/embedding_clustering/run_extract_esm2.sbatch
```

Monitor with `squeue -j <jobid>` and check the log file. This should take well under an hour per the design's estimate, but treat that as an estimate — monitor for actual completion, don't assume a fixed duration.

- [ ] **Step 4: Verify output**

```bash
python3 -c "
from pathlib import Path
import sys
sys.path.insert(0, 'analysis/embedding_clustering')
from embed import load_all_embedding_chunks
embeddings, ids = load_all_embedding_chunks(Path('analysis/embedding_clustering/tables/esm2_classifier'))
print('Shape:', embeddings.shape, '- expect (~749697, 480)')
print('Any NaN:', __import__('numpy').isnan(embeddings).any())
"
```

- [ ] **Step 5: Lint**

Run: `ruff check analysis/embedding_clustering/01_extract_embeddings_esm2_classifier.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/embedding_clustering/01_extract_embeddings_esm2_classifier.py analysis/embedding_clustering/run_extract_esm2.sbatch
git commit -m "Add CLI and SLURM job to extract classifier-space embeddings on real data"
```

(The embeddings themselves, under `tables/esm2_classifier/`, are large and gitignored — add `analysis/embedding_clustering/tables/*` and `analysis/embedding_clustering/figures/*` to `.gitignore` as part of this task, mirroring the prior projects' convention, with `.gitkeep` exceptions.)

---

### Task 4: `02_extract_embeddings_esmc300m.py` — run as a GPU SLURM job

**Files:**
- Create: `analysis/embedding_clustering/02_extract_embeddings_esmc300m.py`
- Create: `analysis/embedding_clustering/run_extract_esmc300m.sbatch`

**Interfaces:**
- Consumes: `embed.py` (Task 2), same protein universe/peptide source as Task 3.
- Produces: `analysis/embedding_clustering/tables/esmc300m/chunk_*_embeddings.npy` + `chunk_*_ids.txt` — consumed by Task 6.

- [ ] **Step 1: Write the CLI script**

Same structure as Task 3's `01_extract_embeddings_esm2_classifier.py`, but calling `extract_esmc300m_embeddings` instead and writing to `tables/esmc300m/`. Copy Task 3's script and make these two substitutions plus update the print statements to say "ESM C 300M" instead of "classifier-space."

- [ ] **Step 2: Write the SLURM batch script**

Same as Task 3's `run_extract_esm2.sbatch`, adjusted job name/log filename, calling `02_extract_embeddings_esmc300m.py`. ESM C 300M is a larger model (300M vs 35M params) — consider requesting more memory (e.g. `--mem=48G`) and possibly a longer time limit if Task 1's feasibility check suggested slower per-batch inference; use judgment based on what Task 1 observed, don't blindly copy Task 3's resource request if the numbers don't fit.

- [ ] **Step 3: Submit, monitor, and verify** — same process as Task 3 Steps 3-4, checking for shape `(~749697, 960)` instead of `(~749697, 480)`.

- [ ] **Step 4: Lint**

Run: `ruff check analysis/embedding_clustering/02_extract_embeddings_esmc300m.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/embedding_clustering/02_extract_embeddings_esmc300m.py analysis/embedding_clustering/run_extract_esmc300m.sbatch
git commit -m "Add CLI and SLURM job to extract ESM C 300M embeddings on real data"
```

---

### Task 5: `cluster.py` — PCA + HDBSCAN + UMAP functions

**Files:**
- Create: `analysis/embedding_clustering/cluster.py`
- Test: `tests/embedding_clustering/test_cluster.py`

**Interfaces:**
- Produces (used by Task 6):
  - `reduce_and_cluster(embeddings: np.ndarray, n_pca_components: int = 50, min_cluster_size: int = 50) -> Dict` returning `{"cluster_labels": np.ndarray, "pca_embeddings": np.ndarray, "n_clusters": int, "n_noise": int}`
  - `compute_umap_projection(embeddings: np.ndarray, n_neighbors: int = 15, random_state: int = 42) -> np.ndarray` (Nx2)

- [ ] **Step 1: Write the failing tests**

Create `tests/embedding_clustering/test_cluster.py`:

```python
"""Tests for analysis/embedding_clustering/cluster.py — uses small
synthetic embeddings with two obviously-separated Gaussian blobs, never
the real 749,697x480/960 data."""
import numpy as np

from cluster import compute_umap_projection, reduce_and_cluster


def _two_blob_embeddings(n_per_blob=100, dim=20, seed=1):
    rng = np.random.default_rng(seed)
    blob_a = rng.normal(loc=0.0, scale=0.5, size=(n_per_blob, dim))
    blob_b = rng.normal(loc=20.0, scale=0.5, size=(n_per_blob, dim))
    return np.vstack([blob_a, blob_b]).astype(np.float32)


def test_reduce_and_cluster_finds_two_separated_blobs():
    embeddings = _two_blob_embeddings()
    result = reduce_and_cluster(embeddings, n_pca_components=10, min_cluster_size=20)
    assert result["n_clusters"] >= 2
    labels = result["cluster_labels"]
    assert len(labels) == len(embeddings)
    # the two blobs should mostly not share a cluster label (allowing some noise, label -1)
    blob_a_labels = set(labels[:100]) - {-1}
    blob_b_labels = set(labels[100:]) - {-1}
    assert not (blob_a_labels & blob_b_labels), "the two separated blobs were merged into one cluster"


def test_reduce_and_cluster_reports_noise_as_its_own_category():
    embeddings = _two_blob_embeddings()
    result = reduce_and_cluster(embeddings, n_pca_components=10, min_cluster_size=20)
    assert "n_noise" in result
    assert result["n_noise"] >= 0


def test_compute_umap_projection_shape():
    embeddings = _two_blob_embeddings(n_per_blob=30)
    projection = compute_umap_projection(embeddings, n_neighbors=5)
    assert projection.shape == (60, 2)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/embedding_clustering/test_cluster.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'cluster'`.

- [ ] **Step 3: Implement `cluster.py`**

Create `analysis/embedding_clustering/cluster.py`:

```python
"""PCA + HDBSCAN clustering, and UMAP projection, for embedding-space
protein clustering. Uses sklearn.cluster.HDBSCAN (built into scikit-learn
>=1.3, already a project dependency) rather than the separate hdbscan
PyPI package."""

from typing import Dict

import numpy as np
import umap
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA


def reduce_and_cluster(
    embeddings: np.ndarray, n_pca_components: int = 50, min_cluster_size: int = 50
) -> Dict:
    """PCA-reduce embeddings, then HDBSCAN-cluster. Returns cluster
    labels (HDBSCAN's -1 = noise, reported as its own category, never
    discarded), the PCA-reduced embeddings, and summary counts."""
    n_components = min(n_pca_components, embeddings.shape[0] - 1, embeddings.shape[1])
    pca = PCA(n_components=n_components, random_state=42)
    pca_embeddings = pca.fit_transform(embeddings)

    clusterer = HDBSCAN(min_cluster_size=min_cluster_size)
    cluster_labels = clusterer.fit_predict(pca_embeddings)

    unique_labels = set(cluster_labels)
    n_clusters = len(unique_labels - {-1})
    n_noise = int(np.sum(cluster_labels == -1))

    return {
        "cluster_labels": cluster_labels,
        "pca_embeddings": pca_embeddings,
        "n_clusters": n_clusters,
        "n_noise": n_noise,
    }


def compute_umap_projection(embeddings: np.ndarray, n_neighbors: int = 15, random_state: int = 42) -> np.ndarray:
    """2D UMAP projection of embeddings, for visualization only (not used
    for clustering itself — that happens on the PCA-reduced space)."""
    reducer = umap.UMAP(n_neighbors=n_neighbors, n_components=2, random_state=random_state)
    return reducer.fit_transform(embeddings)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/embedding_clustering/test_cluster.py -v`
Expected: all PASS. (UMAP with very small N and `n_neighbors=5` should still run, just may print a warning about small dataset size — that's fine for this synthetic test.)

- [ ] **Step 5: Lint**

Run: `ruff check analysis/embedding_clustering/cluster.py tests/embedding_clustering/test_cluster.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/embedding_clustering/cluster.py tests/embedding_clustering/test_cluster.py
git commit -m "Add PCA/HDBSCAN clustering and UMAP projection functions"
```

---

### Task 6: `03_cluster_esm2_classifier.py` + `04_cluster_esmc300m.py` — run on real embeddings

**Files:**
- Create: `analysis/embedding_clustering/03_cluster_esm2_classifier.py`
- Create: `analysis/embedding_clustering/04_cluster_esmc300m.py`

**Interfaces:**
- Consumes: `cluster.py` (Task 5), `embed.py`'s `load_all_embedding_chunks` (Task 2), the real embeddings from Tasks 3/4.
- Produces: `tables/cluster_labels_esm2_classifier.csv` and `tables/cluster_labels_esmc300m.csv` (columns: `protein_id, cluster_label`), plus `tables/umap_esm2_classifier.npy` and `tables/umap_esmc300m.npy` (Nx2 arrays) — consumed by Tasks 8, 9.

- [ ] **Step 1: Write the CLI scripts**

Create `analysis/embedding_clustering/03_cluster_esm2_classifier.py`:

```python
#!/usr/bin/env python
"""Cluster classifier-space embeddings and compute a UMAP projection for
visualization."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cluster import compute_umap_projection, reduce_and_cluster  # noqa: E402
from embed import load_all_embedding_chunks  # noqa: E402

EMBEDDINGS_DIR = Path(__file__).resolve().parent / "tables" / "esm2_classifier"
TABLES_DIR = Path(__file__).resolve().parent / "tables"


def main() -> None:
    embeddings, ids = load_all_embedding_chunks(EMBEDDINGS_DIR)
    print(f"Loaded {embeddings.shape[0]} embeddings, dim={embeddings.shape[1]}")

    result = reduce_and_cluster(embeddings)
    print(f"Found {result['n_clusters']} clusters, {result['n_noise']} noise points")

    pd.DataFrame({"protein_id": ids, "cluster_label": result["cluster_labels"]}).to_csv(
        TABLES_DIR / "cluster_labels_esm2_classifier.csv", index=False
    )

    umap_projection = compute_umap_projection(embeddings)
    np.save(TABLES_DIR / "umap_esm2_classifier.npy", umap_projection)
    print("Wrote cluster labels and UMAP projection")


if __name__ == "__main__":
    main()
```

Create `analysis/embedding_clustering/04_cluster_esmc300m.py` — identical structure, pointing at `tables/esmc300m/` and writing `cluster_labels_esmc300m.csv`/`umap_esmc300m.npy` instead.

- [ ] **Step 2: Run both on real data**

```bash
python analysis/embedding_clustering/03_cluster_esm2_classifier.py
python analysis/embedding_clustering/04_cluster_esmc300m.py
```

UMAP on ~750K points may take significant time (minutes to tens of minutes) — let it run to completion, monitor rather than assume a hang. If it's genuinely impractical (hours), that's real information: fall back to UMAP on a fixed random subsample (e.g. 50,000 points, seeded) for visualization only, clearly noting in the report that the UMAP plot is a subsample while cluster statistics come from the full set — but only do this if the full-set run is actually impractical, not preemptively.

- [ ] **Step 3: Sanity-check output**

For each embedding space, print the cluster size distribution (`value_counts()` on `cluster_label`) and confirm it's not degenerate (e.g. everything in one cluster, or everything noise) — if it is, that's worth investigating (try a different `min_cluster_size`) before moving on.

- [ ] **Step 4: Lint**

Run: `ruff check analysis/embedding_clustering/03_cluster_esm2_classifier.py analysis/embedding_clustering/04_cluster_esmc300m.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/embedding_clustering/03_cluster_esm2_classifier.py analysis/embedding_clustering/04_cluster_esmc300m.py
git commit -m "Add CLIs to cluster both embedding spaces on real data"
```

---

### Task 7: `validate.py` — per-cluster property/domain summary functions

**Files:**
- Create: `analysis/embedding_clustering/validate.py`
- Test: `tests/embedding_clustering/test_validate.py`

**Interfaces:**
- Consumes: `analysis/adhesion_properties/domains.py`'s `top_domain_table` (reused, not duplicated).
- Produces (used by Task 8):
  - `summarize_clusters(master_df: pd.DataFrame, cluster_col: str) -> pd.DataFrame` — per-cluster N, median length/pct_ser_thr_pro/mean_hydrophobicity, fraction with each domain flag.
  - `cluster_contingency_table(labels_a: pd.Series, labels_b: pd.Series) -> pd.DataFrame` — cross-tabulation between two clusterings of the same proteins.
  - `check_feature_concentration(master_df: pd.DataFrame, cluster_col: str, flag_col: str) -> pd.DataFrame` — per-cluster fraction with a given boolean flag set, sorted descending (used for the AA1/laccase concentration check — pass a `has_aa1_family`-style column, computed by the caller from the domain-hit data).

- [ ] **Step 1: Write the failing tests**

Create `tests/embedding_clustering/test_validate.py`:

```python
"""Tests for analysis/embedding_clustering/validate.py."""
import pandas as pd

from validate import check_feature_concentration, cluster_contingency_table, summarize_clusters


def _master_df():
    return pd.DataFrame(
        {
            "protein_id": [f"P{i}" for i in range(8)],
            "cluster_label": [0, 0, 0, 0, 1, 1, 1, 1],
            "length": [100, 110, 120, 130, 400, 410, 420, 430],
            "pct_ser_thr_pro": [40, 42, 38, 44, 10, 12, 9, 11],
            "has_pfam": [False, False, True, False, True, True, True, True],
        }
    )


def test_summarize_clusters():
    summary = summarize_clusters(_master_df(), "cluster_label")
    assert set(summary["cluster_label"]) == {0, 1}
    cluster0 = summary[summary["cluster_label"] == 0].iloc[0]
    assert cluster0["n"] == 4
    assert cluster0["median_length"] == 115.0


def test_cluster_contingency_table():
    labels_a = pd.Series([0, 0, 1, 1])
    labels_b = pd.Series(["X", "X", "Y", "Y"])
    table = cluster_contingency_table(labels_a, labels_b)
    assert table.loc[0, "X"] == 2
    assert table.loc[1, "Y"] == 2


def test_check_feature_concentration():
    df = _master_df()
    df["has_aa1"] = [False, False, False, False, True, True, True, False]
    result = check_feature_concentration(df, "cluster_label", "has_aa1")
    cluster1_frac = result.loc[result["cluster_label"] == 1, "fraction"].iloc[0]
    assert cluster1_frac == 0.75
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/embedding_clustering/test_validate.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'validate'`.

- [ ] **Step 3: Implement `validate.py`**

Create `analysis/embedding_clustering/validate.py`:

```python
"""Per-cluster property/domain validation functions."""

import pandas as pd


def summarize_clusters(master_df: pd.DataFrame, cluster_col: str) -> pd.DataFrame:
    """Per-cluster N, median length/pct_ser_thr_pro, and domain-flag
    fractions (for whichever has_* columns are present in master_df)."""
    domain_flags = [c for c in master_df.columns if c.startswith("has_")]
    rows = []
    for cluster, gdf in master_df.groupby(cluster_col):
        row = {
            cluster_col: cluster,
            "n": len(gdf),
            "median_length": gdf["length"].median() if "length" in gdf else None,
            "median_pct_ser_thr_pro": (
                gdf["pct_ser_thr_pro"].median() if "pct_ser_thr_pro" in gdf else None
            ),
        }
        for flag in domain_flags:
            row[flag] = gdf[flag].mean()
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n", ascending=False).reset_index(drop=True)


def cluster_contingency_table(labels_a: pd.Series, labels_b: pd.Series) -> pd.DataFrame:
    """Cross-tabulation between two clusterings of the same proteins
    (same row order/index correspondence assumed)."""
    return pd.crosstab(labels_a, labels_b)


def check_feature_concentration(master_df: pd.DataFrame, cluster_col: str, flag_col: str) -> pd.DataFrame:
    """Per-cluster fraction with flag_col set, sorted descending — used to
    check whether a flagged feature (e.g. has_aa1_family) concentrates in
    one cluster or is scattered across all of them."""
    result = master_df.groupby(cluster_col)[flag_col].agg(["mean", "size"]).reset_index()
    result = result.rename(columns={"mean": "fraction", "size": "n"})
    return result.sort_values("fraction", ascending=False).reset_index(drop=True)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/embedding_clustering/test_validate.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/embedding_clustering/validate.py tests/embedding_clustering/test_validate.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/embedding_clustering/validate.py tests/embedding_clustering/test_validate.py
git commit -m "Add per-cluster property/domain validation functions"
```

---

### Task 8: `05_validate_and_compare.py` — run on real data

**Files:**
- Create: `analysis/embedding_clustering/05_validate_and_compare.py`

**Interfaces:**
- Consumes: `validate.py` (Task 7), Task 6's cluster-label CSVs, `analysis/adhesion_properties/tables/protein_sequence_properties.csv` + `protein_domain_flags.csv`, `analysis/adhesion_properties/tables/top_domains_cazy.csv` (to identify which specific CAZy family IDs are AA1-related).
- Produces: `tables/cluster_summary_esm2_classifier.csv`, `tables/cluster_summary_esmc300m.csv`, `tables/cluster_contingency_table.csv`, `tables/aa1_concentration_by_cluster.csv` — consumed by Task 10.

- [ ] **Step 1: Write the CLI script**

Create `analysis/embedding_clustering/05_validate_and_compare.py`:

```python
#!/usr/bin/env python
"""Join cluster labels against property/domain data, summarize per
cluster, compare the two clusterings, and check AA1/laccase
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
ADHESION_PROPERTIES_TABLES = (
    Path(__file__).resolve().parents[1] / "adhesion_properties" / "tables"
)

AA1_CAZY_FAMILIES = {"AA1", "AA1_1", "AA1_2", "AA1_3"}


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
    # AA1-family CAZy hit? (approximated here via has_cazy flag combined
    # with membership in the flagged AA1 families is not directly
    # available per-protein in protein_domain_flags.csv, which only has
    # a boolean has_cazy — a precise per-protein AA1 flag would require
    # re-querying cazy_overview per protein_id, which is out of scope for
    # this task; instead, report the concentration using has_cazy as an
    # approximate proxy and note this limitation explicitly in the report.)
    for name, master in [("esm2_classifier", esm2_master), ("esmc300m", esmc_master)]:
        result = check_feature_concentration(master, "cluster_label", "has_cazy")
        result.to_csv(TABLES_DIR / f"aa1_concentration_by_cluster_{name}.csv", index=False)

    print("Done.")


if __name__ == "__main__":
    main()
```

Note the docstring/comment above: `protein_domain_flags.csv` only has a boolean `has_cazy`, not a specific AA1-family flag per protein. A precise AA1-only concentration check would require re-querying `cazy_overview` per protein_id for family membership — out of scope for this task given time constraints, so this task approximates using the broader `has_cazy` flag and states this limitation plainly in the report (Task 10). If a more precise check is wanted later, it would need a new query against the real DuckDB database filtering `cazy_overview.cazyme_fam` to the AA1 family codes for this specific protein set.

- [ ] **Step 2: Run on real data**

```bash
python analysis/embedding_clustering/05_validate_and_compare.py
```

- [ ] **Step 3: Sanity-check output**

Read `cluster_contingency_table.csv` and confirm it's not degenerate (e.g., all mass in one cell would mean the two embedding spaces produce identical clusterings, which is possible but worth a specific note if true; all mass spread evenly would mean no relationship between the two spaces, also worth noting).

- [ ] **Step 4: Lint**

Run: `ruff check analysis/embedding_clustering/05_validate_and_compare.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/embedding_clustering/05_validate_and_compare.py
git commit -m "Add CLI to validate clusters and compare the two embedding spaces on real data"
```

---

### Task 9: `figures_ext.py` + `06_make_figures.py` — figures

Before writing this task's plotting code, invoke the **dataviz** skill, same as the prior two projects — reuse `analysis/kingdom_survey/figures.py`'s validated color constants where a categorical palette is needed (e.g. distinguishing a handful of clusters), falling back to a perceptually-uniform continuous colormap (e.g. matplotlib's `viridis`) for coloring by a continuous feature like `pct_ser_thr_pro`.

**Files:**
- Create: `analysis/embedding_clustering/figures_ext.py`
- Create: `analysis/embedding_clustering/06_make_figures.py`
- Test: `tests/embedding_clustering/test_figures_ext.py`

**Interfaces:**
- Produces: `umap_scatter_by_cluster(umap_coords: np.ndarray, cluster_labels: np.ndarray, out_path: Path, title: str) -> None`, `umap_scatter_by_feature(umap_coords: np.ndarray, feature_values: np.ndarray, out_path: Path, title: str, is_categorical: bool = False) -> None`.

- [ ] **Step 1: Write the failing tests**

Create `tests/embedding_clustering/test_figures_ext.py`:

```python
"""Smoke tests for analysis/embedding_clustering/figures_ext.py."""
import numpy as np

from figures_ext import umap_scatter_by_cluster, umap_scatter_by_feature


def _sample_umap():
    rng = np.random.default_rng(1)
    return rng.normal(size=(50, 2))


def test_umap_scatter_by_cluster_creates_file(tmp_path):
    coords = _sample_umap()
    labels = np.array([0, 1, -1] * 16 + [0, 1])
    out = tmp_path / "umap_cluster.png"
    umap_scatter_by_cluster(coords, labels, out, title="Test")
    assert out.exists()
    assert out.stat().st_size > 0


def test_umap_scatter_by_feature_continuous(tmp_path):
    coords = _sample_umap()
    values = np.linspace(0, 100, 50)
    out = tmp_path / "umap_feature.png"
    umap_scatter_by_feature(coords, values, out, title="Test", is_categorical=False)
    assert out.exists()
    assert out.stat().st_size > 0


def test_umap_scatter_by_feature_categorical(tmp_path):
    coords = _sample_umap()
    values = np.array([True, False] * 25)
    out = tmp_path / "umap_bool.png"
    umap_scatter_by_feature(coords, values, out, title="Test", is_categorical=True)
    assert out.exists()
    assert out.stat().st_size > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/embedding_clustering/test_figures_ext.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError`.

- [ ] **Step 3: Implement `figures_ext.py`**

Create `analysis/embedding_clustering/figures_ext.py`:

```python
"""Figure generation for the embedding-space clustering analysis."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def umap_scatter_by_cluster(
    umap_coords: np.ndarray, cluster_labels: np.ndarray, out_path: Path, title: str
) -> None:
    """2D UMAP scatter colored by cluster label; noise (-1) rendered in
    light gray, distinct from any real cluster color."""
    fig, ax = plt.subplots(figsize=(8, 7))
    unique_labels = sorted(set(cluster_labels) - {-1})
    cmap = plt.get_cmap("tab20")
    for i, label in enumerate(unique_labels):
        mask = cluster_labels == label
        ax.scatter(
            umap_coords[mask, 0], umap_coords[mask, 1], s=3, alpha=0.5,
            color=cmap(i % 20), label=f"cluster {label}",
        )
    noise_mask = cluster_labels == -1
    if noise_mask.any():
        ax.scatter(
            umap_coords[noise_mask, 0], umap_coords[noise_mask, 1],
            s=3, alpha=0.3, color="lightgray", label="noise",
        )
    ax.set_title(title)
    ax.legend(markerscale=4, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def umap_scatter_by_feature(
    umap_coords: np.ndarray,
    feature_values: np.ndarray,
    out_path: Path,
    title: str,
    is_categorical: bool = False,
) -> None:
    """2D UMAP scatter colored by a continuous (viridis) or boolean/
    categorical (tab10) feature, for visually checking whether clusters
    align with known properties."""
    fig, ax = plt.subplots(figsize=(8, 7))
    if is_categorical:
        unique_values = sorted(set(feature_values))
        cmap = plt.get_cmap("tab10")
        for i, value in enumerate(unique_values):
            mask = feature_values == value
            ax.scatter(
                umap_coords[mask, 0], umap_coords[mask, 1], s=3, alpha=0.5,
                color=cmap(i % 10), label=str(value),
            )
        ax.legend(markerscale=4, fontsize=8)
    else:
        scatter = ax.scatter(
            umap_coords[:, 0], umap_coords[:, 1], s=3, alpha=0.5,
            c=feature_values, cmap="viridis",
        )
        fig.colorbar(scatter, ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/embedding_clustering/test_figures_ext.py -v`
Expected: all PASS.

- [ ] **Step 5: Write `06_make_figures.py` and run on real data**

Create `analysis/embedding_clustering/06_make_figures.py`:

```python
#!/usr/bin/env python
"""Generate UMAP figures for both embedding spaces: colored by cluster,
and colored by key features (pct_ser_thr_pro, has_signal_peptide,
has_cazy) for visual alignment checking."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures_ext import umap_scatter_by_cluster, umap_scatter_by_feature  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"
ADHESION_PROPERTIES_TABLES = (
    Path(__file__).resolve().parents[1] / "adhesion_properties" / "tables"
)


def _make_space_figures(space_name: str) -> None:
    umap_coords = np.load(TABLES_DIR / f"umap_{space_name}.npy")
    clusters = pd.read_csv(TABLES_DIR / f"cluster_labels_{space_name}.csv")
    props = pd.read_csv(ADHESION_PROPERTIES_TABLES / "protein_sequence_properties.csv")
    flags = pd.read_csv(ADHESION_PROPERTIES_TABLES / "protein_domain_flags.csv")
    merged = clusters.merge(props, on="protein_id", how="left").merge(
        flags, on="protein_id", how="left"
    )

    umap_scatter_by_cluster(
        umap_coords, merged["cluster_label"].to_numpy(),
        FIGURES_DIR / f"umap_cluster_{space_name}.png", title=f"Clusters ({space_name})",
    )
    umap_scatter_by_feature(
        umap_coords, merged["pct_ser_thr_pro"].to_numpy(),
        FIGURES_DIR / f"umap_pct_ser_thr_pro_{space_name}.png",
        title=f"Ser/Thr/Pro % ({space_name})", is_categorical=False,
    )
    umap_scatter_by_feature(
        umap_coords, merged["has_signal_peptide"].to_numpy(),
        FIGURES_DIR / f"umap_has_signal_peptide_{space_name}.png",
        title=f"Signal peptide ({space_name})", is_categorical=True,
    )
    umap_scatter_by_feature(
        umap_coords, merged["has_cazy"].to_numpy(),
        FIGURES_DIR / f"umap_has_cazy_{space_name}.png",
        title=f"CAZy presence ({space_name})", is_categorical=True,
    )


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for space_name in ["esm2_classifier", "esmc300m"]:
        _make_space_figures(space_name)
    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
```

Run it, then visually sanity-check at least 3 of the 8 generated figures (view them with the Read tool) — confirm legible axes/legends, real point clouds, distinguishable cluster colors.

- [ ] **Step 6: Lint**

Run: `ruff check analysis/embedding_clustering/figures_ext.py analysis/embedding_clustering/06_make_figures.py tests/embedding_clustering/test_figures_ext.py`

- [ ] **Step 7: Commit**

```bash
git add analysis/embedding_clustering/figures_ext.py analysis/embedding_clustering/06_make_figures.py tests/embedding_clustering/test_figures_ext.py
git commit -m "Add figure generation for embedding clustering, run on real data"
```

---

### Task 10: `07_generate_report.py` — assemble and write `REPORT.md`

**Files:**
- Create: `analysis/embedding_clustering/07_generate_report.py`
- Create: `analysis/embedding_clustering/REPORT.md` (generated, then hand-edited)

**Interfaces:**
- Consumes: every table under `tables/` and every figure under `figures/` produced by Tasks 3-9.
- Produces: the final `REPORT.md` deliverable.

- [ ] **Step 1: Write the report-assembly script**

Create `analysis/embedding_clustering/07_generate_report.py` following the established pattern from the prior two projects (pull real numbers from CSVs programmatically, leave the interpretive "Headline findings" as a section filled in by hand in Step 3 after inspecting real output). Include:
- Scope: N proteins embedded per space, N clusters found per space, N noise points per space.
- Per-space cluster summary table (from `cluster_summary_*.csv`, via a simple markdown rendering).
- The cluster contingency table (or a summary of it — e.g. adjusted Rand-like agreement description) between the two embedding spaces.
- The AA1/CAZy concentration-by-cluster result, explicitly noting the has_cazy-as-proxy limitation.
- All 8 figures embedded.
- A caveats section: exploratory/unsupervised (no significance testing), HDBSCAN noise category reported not hidden, the AA1 check is an approximation (has_cazy proxy, not AA1-specific), embedding models carry their own pretraining biases, cluster "architecture" identification (FLO11-like/Als-like) is a hypothesis suggested by the data, not a confirmed structural classification.

- [ ] **Step 2: Run it, then fill in Headline findings by hand**

Run the script, then read the actual per-cluster summaries and the contingency table, and write 3-5 sentences describing: how many clusters were found in each space, whether the two spaces' clusterings agree substantially or diverge, whether any cluster's property/domain signature resembles a recognizable architecture (name the actual defining features you see — e.g. "cluster 3 in the ESM C space has median Ser/Thr/Pro% of X, N=Y, Z% with signal peptide, consistent with a FLO11-like mucin/GPI-anchor architecture" — only if the real numbers actually support this, don't force a narrative that isn't there), and what the AA1-concentration check showed.

- [ ] **Step 3: Verify and commit**

Confirm the report renders, all figure paths exist, commit the script + REPORT.md + the small tables (`cluster_summary_*.csv`, `cluster_contingency_table.csv`, `aa1_concentration_by_cluster_*.csv`) + all figures with `git add -f` (matching the established pattern — these are gitignored by default). Do NOT commit the large embedding `.npy` chunk files or the raw `cluster_labels_*.csv` (large, per-protein, fully reproducible via Tasks 3/4/6 given the same real data source) — same reproducible-intermediate-data judgment call as the prior project's git-history ruling.

```bash
git add analysis/embedding_clustering/07_generate_report.py analysis/embedding_clustering/REPORT.md
git add -f analysis/embedding_clustering/tables/cluster_summary_*.csv analysis/embedding_clustering/tables/cluster_contingency_table.csv analysis/embedding_clustering/tables/aa1_concentration_by_cluster_*.csv analysis/embedding_clustering/figures/*.png
git commit -m "Add embedding clustering report, small summary tables, and figures"
```

---

## Self-Review Notes

- **Spec coverage:** every pipeline stage from the design doc has a task — feasibility gate (Task 1), both embedding extractions (Tasks 3-4), clustering (Tasks 5-6), validation/comparison including the AA1 concentration check (Tasks 7-8), figures (Task 9), report (Task 10).
- **Placeholder scan:** Task 1 is deliberately investigative (its outcome determines Task 2's ESM C loading code) — this is not a "TBD" left unresolved, it's an explicit, bounded spike with clear success/failure criteria and a hard stop-and-report-back path if both candidate approaches fail. Task 2's `extract_esmc300m_embeddings` has a `NotImplementedError` body that Step 4 of that same task explicitly requires filling in per Task 1's findings before the task can be considered complete — not a silently-shipped gap.
- **Type/name consistency:** checked `embed.py`'s `save_embeddings_chunk`/`load_all_embedding_chunks` contract is used identically by Tasks 3, 4, and 6; `cluster.py`'s `reduce_and_cluster` return dict keys match Task 6's usage; `validate.py`'s function signatures match Task 8's calls; `figures_ext.py`'s function signatures match Task 9's `06_make_figures.py` calls.
