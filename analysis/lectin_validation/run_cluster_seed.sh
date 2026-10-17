#!/bin/bash -l
# Cluster the truth sequences (MMseqs2 easy-cluster, 30% identity, 50% coverage) and search them
# against the Pfam seed sequences of the 20 models (leakage). Login node, small input.
# Usage: run_cluster_seed.sh OUT_DIR SEQ_FAA SEEDS_FAA_GZ TMP_DIR
# No BASH_SOURCE is used. All paths come from arguments.
set -euo pipefail
OUT=${1:?out dir}
SEQ=${2:?sequences.faa}
SEEDS=${3:?seeds.faa.gz}
TMP=${4:?tmp dir}
module load MMseqs2/17-b804f
mkdir -p "$OUT" "$TMP"
mmseqs easy-cluster "$SEQ" "$OUT/clu" "$TMP/clu_tmp" --min-seq-id 0.3 -c 0.5 --threads 4 -v 1
# seeds: unzip to a plain file in TMP, make a search against it
zcat "$SEEDS" > "$TMP/seeds.faa"
# query = truth sequence, target = seed sequence. Identity at 0.3, target coverage 0.5.
mmseqs easy-search "$SEQ" "$TMP/seeds.faa" "$OUT/seed_hits.m8" "$TMP/search_tmp" \
  --min-seq-id 0.3 --cov-mode 1 -c 0.5 -s 7.5 -e 1e-3 --threads 4 -v 1 \
  --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"
mmseqs version > "$OUT/mmseqs_version.txt"
rm -rf "$TMP"
