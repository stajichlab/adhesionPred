#!/bin/bash -l
# Partition/account: `-p exfab -A exfab --gres=gpu:ada6000:1` (node gpu12,
# 2x ada6000) from the START -- this is the hard-won lesson carried
# forward from Task 3 (ESM2-classifier extraction), which burned through
# `short_gpu` (too slow, 2h limit), `gpu` bare `gpu:1` (landed on an
# incompatible P100 node, sm_60 -- "CUDA error: no kernel image is
# available"), `gpu --gres=gpu:a100:1` (>1h pending, only 8 a100 GPUs
# total), and `gpu --constraint=gpu_latest` (~10h pending, pure
# account-fairshare exhaustion on the shared `gpu` partition) before
# landing on `exfab`, which started RUNNING within seconds both times it
# was used. `scontrol show node gpu12` was checked before this job's
# submission: 1 of gpu12's 2 ada6000 GPUs was allocated to another
# exfab user at check time, 1 free -- proceeding with exfab per the
# task brief's guidance (try it first; only 2 GPUs total but
# PriorityTier=3 should queue-jump reasonably fast even if briefly
# contended).
#
# Model/memory: ESM Cambrian 300M is ~9x the classifier's 35M-parameter
# ESM2 model. Task 3's measured ESM2 rate was ~110-175s/chunk (150
# chunks, 5000 proteins/chunk). A smoke test (2-3 chunks) was run on an
# actual exfab/gpu12 node before sizing this job's wall-clock budget --
# see the SDD task-4-report.md for the measured per-chunk rate and the
# reasoning behind the --time value below. Requesting --mem=48G (vs.
# Task 3's 32G) since the larger model and its activations need more
# headroom; --time budgeted generously since even Task 3's smaller model
# needed one resume cycle after a 6h TIMEOUT -- the resume-skip logic in
# 02_extract_embeddings_esmc300m.py (hardened over Task 3's version,
# with an ids/embeddings row-count consistency check, not just file
# existence) makes a TIMEOUT-and-resubmit cycle a non-issue if it
# happens again here.
#SBATCH -p exfab
#SBATCH -A exfab
#SBATCH --gres=gpu:ada6000:1
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH -J embed_esmc300m
#SBATCH -o logs/embed_esmc300m.%j.log
#SBATCH -e logs/embed_esmc300m.%j.log

# Directory resolution: SLURM copies the batch script into a per-job
# spool directory before running it, so resolving relative to "$0"
# (e.g. `dirname "$(readlink -f "$0")"`) resolves to
# /var/spool/slurmd/job<id>/.. -- NOT the repo -- confirmed broken on
# this cluster during Task 3. SLURM_SUBMIT_DIR is set by sbatch to the
# directory sbatch was invoked from, and is the reliable way to get
# back to the repo root (this script must always be submitted via
# `sbatch` from the repo/worktree root).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR not set}"

# IMPORTANT: invoke the ESM-C-specific Python 3.11 venv explicitly, NOT
# bare `python` (which resolves to the project's main Python 3.9 env and
# lacks the `esm` PyPI package -- and even if it were installed there,
# fair-esm already occupies the `esm` module name in that env; see Task
# 1's report and embed.py's SCOPE NOTE docstring).
.venv_esmc/bin/python analysis/embedding_clustering/02_extract_embeddings_esmc300m.py
