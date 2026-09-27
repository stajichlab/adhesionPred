#!/usr/bin/env bash
# Build the pilot environment (run once on an HPCC login node, from this directory).
# EvolutionaryScale `esm` (ESM C) and HuggingFace `transformers` (ESM-2) share one env;
# fair-esm is not installed here because it also claims `import esm`.
set -euo pipefail
cd "$(dirname "$0")"
uv venv -q -p 3.11 .venv
# CUDA 12.6 wheels: the default PyPI torch may target a newer CUDA than the node driver supports
# httpx: imported by esm.sdk but not declared by the esm package
uv pip install -q -p .venv/bin/python esm httpx transformers biopython numpy
# esm pulls torch/torchvision from PyPI; replace both with a matched cu126 pair
uv pip install -q -p .venv/bin/python --reinstall-package torch --reinstall-package torchvision torch torchvision \
  --index-url https://download.pytorch.org/whl/cu126
mkdir -p data logs
[ -s data/S288C_orf_trans_all.fasta.gz ] || curl -sfL -o data/S288C_orf_trans_all.fasta.gz \
  http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz
# Download weights on the login node (compute nodes run with HF_HUB_OFFLINE=1).
# Download only: loading the models here exceeds the login-node memory cap.
.venv/bin/python - <<'PY'
from huggingface_hub import snapshot_download
for repo in ["facebook/esm2_t30_150M_UR50D", "facebook/esm2_t33_650M_UR50D",
             "EvolutionaryScale/esmc-300m-2024-12"]:
    print(snapshot_download(repo))
PY
