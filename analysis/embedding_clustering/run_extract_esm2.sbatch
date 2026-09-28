#!/bin/bash -l
# Partition: NOT short_gpu. A real 3-chunk (15,000-protein) smoke test on
# a short_gpu ada6000 node measured a steady-state ~41.8s/chunk DB fetch
# (analysis/adhesion_properties/db.py's fetch_lengths against the 278GB
# functional-annotation DuckDB file -- does not amortize across chunks,
# it's a real per-query cost) + ~64.5s/chunk ESM2 extraction, i.e.
# ~106.3s/chunk. At 150 chunks (749,697 proteins / 5,000 per chunk) that
# is ~266 min (~4.4h) steady-state -- genuinely over short_gpu's 2h limit,
# not just close to it. Using `gpu` (7-day limit, same node pool
# including a100/h100) instead, with a 6h wall-clock budget for margin.
# GPU type: NOT a bare `gpu:1`. A first real submission landed on gpu05
# (Tesla P100, CUDA capability sm_60) and every batch failed with "CUDA
# error: no kernel image is available for execution on the device" --
# the installed PyTorch build only supports sm_70+, so p100/k80 nodes in
# this partition's pool are unusable for this job. Pinning to a100 alone
# (`--gres=gpu:a100:1`) fixed correctness but left the job pending over
# an hour with no node (only 8 a100 GPUs total, heavily contended).
# `--constraint=gpu_latest` is a shared SLURM node feature covering every
# sm_70+ GPU generation this build supports (a100/ada6000/h100/
# blackwell6000 -- 44 GPUs total vs. 8 for a100 alone; ada6000 sm_86 was
# already confirmed compatible in the pre-submission smoke test), giving
# both correctness and a much larger schedulable pool.
#
# UPDATE: job 28148283 (above config, `-p gpu --gres=gpu:1 --constraint=
# gpu_latest`) sat PENDING(Priority) for ~10h with zero movement -- pure
# fairshare exhaustion on the shared `gpu` partition, not fixable by
# widening the GPU-type constraint further. Found a better path: node
# gpu12 (2x ada6000, already confirmed sm_86-compatible) lives in the
# `exfab` partition (`-A exfab -p exfab`), which this account has access
# to, is high-priority (PriorityTier=3), and at submission time had 0 of
# its 2 GPUs actually allocated despite other jobs running there (those
# are CPU-only). Switched to submitting there instead of `gpu`.
#SBATCH -p exfab
#SBATCH -A exfab
#SBATCH --gres=gpu:ada6000:1
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH --time=6:00:00
#SBATCH -J embed_esm2_classifier
#SBATCH -o logs/embed_esm2_classifier.%j.log
#SBATCH -e logs/embed_esm2_classifier.%j.log

# Directory resolution: SLURM copies the batch script into a per-job
# spool directory before running it, so resolving relative to "$0"
# (e.g. `dirname "$(readlink -f "$0")"`) resolves to
# /var/spool/slurmd/job<id>/.. -- NOT the repo -- and was verified broken
# on this cluster. SLURM_SUBMIT_DIR is set by sbatch to the directory
# sbatch was invoked from, and is the reliable way to get back to the
# repo root (this script must always be submitted via `sbatch` from the
# repo/worktree root).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR not set}"

python analysis/embedding_clustering/01_extract_embeddings_esm2_classifier.py
