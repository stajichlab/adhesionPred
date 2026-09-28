#!/usr/bin/env python
"""Pilot: embed reference proteomes with several PLMs and measure throughput.

For each (model, genome) this writes <out>/<model>/<genome>.npz with
  ids      protein ids (input order)
  length   residues after sanitizing
  full     mean over all residues, using overlapping windows for long proteins
  nterm    mean over the first --nterm residues (N-terminal adhesion-domain proxy)
and appends one JSON line of timing/memory stats to <out>/throughput.jsonl.

Pooling excludes BOS/EOS/padding. Chunks are length-sorted and batched by a token
budget; a batch that runs out of memory is split in half and retried.

Models: esmc_300m / esmc_600m (EvolutionaryScale `esm` SDK) and any
facebook/esm2_* checkpoint (HuggingFace transformers), final layer.
"""

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import torch
from Bio import SeqIO

RESIDUES = set("ACDEFGHIKLMNPQRSTVWYXBUZO")


def sanitize(seq):
    seq = seq.upper().replace("*", "").replace("-", "").replace(".", "")
    return "".join(c if c in RESIDUES else "X" for c in seq)


def read_fasta(path):
    handle = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    with handle:
        return [(r.id, sanitize(str(r.seq))) for r in SeqIO.parse(handle, "fasta")]


def windows(n, size, stride):
    if n <= size:
        return [(0, n)]
    starts = list(range(0, n - size, stride)) + [n - size]
    return [(s, s + size) for s in starts]


class HFEsm2:
    def __init__(self, name, device):
        from transformers import AutoTokenizer, EsmModel

        self.tok = AutoTokenizer.from_pretrained(f"facebook/{name}")
        self.model = EsmModel.from_pretrained(f"facebook/{name}", add_pooling_layer=False)
        self.model = self.model.to(device).eval()
        self.device = device
        self.special = torch.tensor(self.tok.all_special_ids, device=device)
        self.max_residues = 1022

    def residue_sums(self, seqs):
        enc = self.tok(seqs, return_tensors="pt", padding=True).to(self.device)
        hidden = self.model(**enc).last_hidden_state
        mask = ~torch.isin(enc["input_ids"], self.special)
        return (hidden * mask.unsqueeze(-1)).float().sum(1), mask.sum(1)


class EsmC:
    def __init__(self, name, device):
        from esm.models.esmc import ESMC

        self.model = ESMC.from_pretrained(name, device=device).eval()
        tok = self.model.tokenizer
        self.special = torch.tensor(
            [tok.cls_token_id, tok.eos_token_id, tok.pad_token_id], device=device
        )
        self.max_residues = 2046

    def residue_sums(self, seqs):
        tokens = self.model._tokenize(seqs)
        hidden = self.model(sequence_tokens=tokens).embeddings
        mask = ~torch.isin(tokens, self.special)
        return (hidden * mask.unsqueeze(-1)).float().sum(1), mask.sum(1)


def load_model(name, device):
    return EsmC(name, device) if name.startswith("esmc") else HFEsm2(name, device)


def embed_genome(model, proteins, nterm, token_budget, device):
    n = len(proteins)
    # chunk = (protein index, target 0=full/1=nterm, sequence text)
    chunks = []
    size = model.max_residues
    for i, (_, seq) in enumerate(proteins):
        for a, b in windows(len(seq), size, size // 2):
            chunks.append((i, 0, seq[a:b]))
        chunks.append((i, 1, seq[:nterm]))
    chunks.sort(key=lambda c: len(c[2]))

    dim = None
    sums = counts = None
    batches, cur = [], []
    for c in chunks:
        if cur and (len(cur) + 1) * (len(c[2]) + 2) > token_budget:
            batches.append(cur)
            cur = []
        cur.append(c)
    if cur:
        batches.append(cur)

    def run(batch):
        nonlocal dim, sums, counts
        try:
            with (
                torch.no_grad(),
                torch.autocast(device.type, dtype=torch.bfloat16, enabled=device.type == "cuda"),
            ):
                s, k = model.residue_sums([c[2] for c in batch])
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            if len(batch) == 1:
                raise
            half = len(batch) // 2
            run(batch[:half])
            run(batch[half:])
            return
        s, k = s.cpu().numpy(), k.cpu().numpy()
        if dim is None:
            dim = s.shape[1]
            sums = np.zeros((2, n, dim), dtype=np.float64)
            counts = np.zeros((2, n), dtype=np.int64)
        for (i, t, _), vec, cnt in zip(batch, s, k):
            sums[t, i] += vec
            counts[t, i] += cnt

    for batch in batches:
        run(batch)
    pooled = sums / np.maximum(counts, 1)[..., None]
    return pooled[0].astype(np.float16), pooled[1].astype(np.float16), len(batches)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--genomes", required=True, help="TSV: name<TAB>fasta path")
    ap.add_argument("--models", default="esmc_300m,esm2_t30_150M_UR50D,esm2_t33_650M_UR50D")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--nterm", type=int, default=300)
    ap.add_argument("--token-budget", type=int, default=24000)
    ap.add_argument("--limit", type=int, default=None, help="first N proteins per genome (testing)")
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    genomes = [line.rstrip("\n").split("\t") for line in open(args.genomes) if line.strip()]
    args.out.mkdir(parents=True, exist_ok=True)
    gpu = torch.cuda.get_device_name(0) if device.type == "cuda" else "cpu"

    for model_name in args.models.split(","):
        t0 = time.time()
        model = load_model(model_name, device)
        load_s = time.time() - t0
        (args.out / model_name).mkdir(exist_ok=True)
        for name, path in genomes:
            proteins = read_fasta(path)[: args.limit]
            if device.type == "cuda":
                torch.cuda.synchronize()
                torch.cuda.reset_peak_memory_stats()
            t1 = time.time()
            full, nterm, n_batches = embed_genome(
                model, proteins, args.nterm, args.token_budget, device
            )
            if device.type == "cuda":
                torch.cuda.synchronize()
            secs = time.time() - t1
            np.savez_compressed(
                args.out / model_name / f"{name}.npz",
                ids=np.array([p[0] for p in proteins]),
                length=np.array([len(p[1]) for p in proteins]),
                full=full,
                nterm=nterm,
            )
            stats = {
                "model": model_name,
                "genome": name,
                "gpu": gpu,
                "proteins": len(proteins),
                "residues": int(sum(len(p[1]) for p in proteins)),
                "seconds": round(secs, 1),
                "proteins_per_s": round(len(proteins) / secs, 1),
                "residues_per_s": round(sum(len(p[1]) for p in proteins) / secs),
                "batches": n_batches,
                "token_budget": args.token_budget,
                "model_load_s": round(load_s, 1),
                "peak_mem_gb": round(torch.cuda.max_memory_allocated() / 1e9, 2)
                if device.type == "cuda"
                else None,
            }
            print(json.dumps(stats), flush=True)
            with open(args.out / "throughput.jsonl", "a") as f:
                f.write(json.dumps(stats) + "\n")
        del model
        if device.type == "cuda":
            torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
