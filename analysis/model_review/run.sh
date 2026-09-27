#!/usr/bin/env bash
# Reproduce the 2026-09-27 model review experiments (CPU is fine; ~1.5 h on a laptop).
# Outputs land in the current working directory (use a scratch dir).
#   Requires: python>=3.10 with torch fair-esm scikit-learn biopython pandas; mmseqs2 on PATH.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
[ -f scer.fasta.gz ] || curl -sL -o scer.fasta.gz \
  http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz
gzip -dc "$REPO"/data/positive/*.gz "$REPO"/data/negative/*.gz > train_all.fa
mmseqs easy-cluster train_all.fa clu tmpmm --min-seq-id 0.3 -c 0.5 --cov-mode 0 -v 1
python "$HERE/01_embed_esm2_8M.py" scer.fasta.gz
python "$HERE/02_cv_and_proteome_eval.py" | tee results.txt
