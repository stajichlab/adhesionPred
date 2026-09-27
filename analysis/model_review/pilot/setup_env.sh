#!/usr/bin/env bash
# Build the pilot environment (run once on an HPCC login node, from this directory).
# EvolutionaryScale `esm` (ESM C) and HuggingFace `transformers` (ESM-2) share one env;
# fair-esm is not installed here because it also claims `import esm`.
set -euo pipefail
cd "$(dirname "$0")"
uv venv -q -p 3.11 .venv
uv pip install -q -p .venv/bin/python "torch>=2.4" esm transformers biopython numpy
mkdir -p data logs
[ -s data/S288C_orf_trans_all.fasta.gz ] || curl -sfL -o data/S288C_orf_trans_all.fasta.gz \
  http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz
# pre-fetch weights on the login node so compute nodes don't need internet
.venv/bin/python - <<'PY'
from transformers import AutoTokenizer, EsmModel
for m in ["facebook/esm2_t30_150M_UR50D", "facebook/esm2_t33_650M_UR50D"]:
    AutoTokenizer.from_pretrained(m); EsmModel.from_pretrained(m, add_pooling_layer=False)
from esm.models.esmc import ESMC
ESMC.from_pretrained("esmc_300m", device="cpu")
print("weights cached")
PY
