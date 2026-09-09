#!/usr/bin/env python
"""Recovery script: embed the proteins that Task 3's real run silently
dropped from the ESM2-classifier embedding table, and save them as an
additional chunk (chunk_00150) alongside chunks 0-149.

Background (see task-3-report.md's dated addendum and progress.md's
"Task 3 CORRECTION" entries for the full investigation): the classifier's
own `get_esm_embeddings` (src/adhesion_predict/embeddings.py) wraps each
batch of sequences in a bare `except Exception: continue`, which drops
the ENTIRE BATCH -- not just the offending sequence -- if any sequence in
it fails to tokenize. The real 749,697-protein run was missing 176
proteins as a result. Diffing against Task 4's ESM-C run (which embedded
all 749,697 with zero missing, since ESM-C's tokenizer handles a broader
character set) identified exactly which 176. Re-running those 176 through
the identical model/function at `batch_size=1` (isolating each sequence
from whatever batch-mate originally caused the failure) recovered 159 of
them; the remaining 17 genuinely contain `J` or `*` characters that are
not in ESM-2's tokenizer alphabet at all, and are not recoverable this way
(confirmed: still fail even in isolation) -- these 17 proteins simply
cannot be embedded in ESM-2-classifier space, a real, permanent limitation
of that space, not a bug to route around further.

This script is idempotent: it recomputes the same "missing" set from the
current on-disk chunks each time (not off a hardcoded id list), so it is
safe to re-run if chunks 0-149 or the ESM-C reference set ever change.
`get_esm_embeddings` itself is deliberately left unmodified -- it is
pre-existing production code for the shipped classifier, out of scope for
this analysis project to alter, per this project's "replicate the
classifier's exact recipe, don't reimplement or alter it" principle.

Usage: run from the repo root (or this worktree's root) with the main
project's Python 3.9 environment -- this uses fair-esm, not ESM-C, so it
must NOT be run under .venv_esmc.
    python analysis/embedding_clustering/recover_missing_esm2_embeddings.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "adhesion_properties"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from adhesion_predict.embeddings import get_esm_embeddings  # noqa: E402
from db import connect, fetch_lengths  # noqa: E402
from embed import save_embeddings_chunk  # noqa: E402

ESM2_DIR = Path(__file__).resolve().parent / "tables" / "esm2_classifier"
ESMC_DIR = Path(__file__).resolve().parent / "tables" / "esmc300m"
RECOVERY_CHUNK_IDX = 150


def _load_ids(chunk_dir: Path) -> set:
    ids = set()
    for f in sorted(chunk_dir.glob("chunk_*_ids.txt")):
        ids.update(line.strip() for line in f.open() if line.strip())
    return ids


def main() -> None:
    esm2_ids = _load_ids(ESM2_DIR)
    esmc_ids = _load_ids(ESMC_DIR)
    missing = sorted(esmc_ids - esm2_ids)
    print(f"proteins present in ESM-C but missing from ESM2-classifier: {len(missing)}")

    if not missing:
        print("nothing to recover")
        return

    con = connect()
    df = fetch_lengths(con, missing)
    con.close()
    sequences = [
        {"id": row["protein_id"], "sequence": row["peptide"]}
        for _, row in df.iterrows()
        if row["peptide"]
    ]

    # batch_size=1 isolates each sequence from whatever batch-mate
    # originally triggered get_esm_embeddings' batch-level exception --
    # mathematically identical model output to any other batch size, only
    # the batching (and therefore which failures collide) differs.
    embeddings, ids = get_esm_embeddings(sequences, model_name="esm2_t12_35M_UR50D", batch_size=1)
    recovered = set(ids)
    still_missing = [pid for pid in missing if pid not in recovered]

    print(f"recovered: {len(ids)} / {len(missing)}")
    if still_missing:
        print(
            f"still unrecoverable (expected -- these contain 'J'/'*' characters "
            f"outside ESM-2's tokenizer alphabet): {len(still_missing)}"
        )
        for pid in still_missing:
            print(f"  {pid}")

    if not ids:
        # Nothing recovered this run -- e.g. a re-run where the only
        # remaining "missing" ids are the permanently-unrecoverable 17.
        # Do NOT call save_embeddings_chunk here: it would overwrite an
        # existing, valid chunk_00150 with an empty one. Only write when
        # there is something real to write.
        print("nothing recovered this run; leaving any existing chunk untouched")
        return

    save_embeddings_chunk(embeddings, ids, ESM2_DIR, RECOVERY_CHUNK_IDX)
    print(f"wrote chunk_{RECOVERY_CHUNK_IDX:05d} to {ESM2_DIR}")


if __name__ == "__main__":
    main()
