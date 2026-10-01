"""Tiny Phase B embedding plan for the J2 runner and assembly tests (CPU, ESM-2 8M)."""

import seqhash
import seqsets
import truth_table
from conftest import load_script

MODEL = "esm2_t6_8M_UR50D"
LONG = "MKLSTA" + "STPSSTSA" * 140  # 1,126 aa: the only sequence with a C-terminal window
SEQS = [
    "MKTLLVAGLLSSAAFA",
    "MSTTSSTTSTPSSTSA" * 4,
    "MKVLAAGIVALLLAAGCSSS" * 3,
    "MQRSLLLAVAALATPAFAAS" * 6,
    "MAEEKKAVEEVKSAGEW" * 2,
    LONG,
]


def make_plan(work, chunk_residues=400):
    """Write unique_sequences.tsv.gz and a dry chunk plan (assumed rate) under work/phaseb."""
    out = work / "phaseb"
    out.mkdir(parents=True)
    unique = seqsets.unique_rows({seqhash.seq_sha256(s): s for s in SEQS})
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    plan = load_script("06_plan_embedding")
    argv = ["--work-dir", str(work), "--rate", "100", "--models", MODEL]
    assert plan.main(argv + ["--chunk-residues", str(chunk_residues)]) == 0
    return out, unique


def run_cpu(work, scratch, **kw):
    import embed_chunks
    import torch

    cpu = torch.device("cpu")
    return embed_chunks.run_job(work, [MODEL], 0, 1, cpu, scratch, batch_size=2, **kw)


def change_row0_sequence(out):
    """Replace the sequence of row 0 in unique_sequences.tsv.gz and keep its stored hash.

    This is the reviewer's probe: a unique set that changed after 06 wrote the plan."""
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    unique[0]["sequence"] = "MKWWWWWWWWAAAA"
    unique[0]["length"] = str(len(unique[0]["sequence"]))
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
