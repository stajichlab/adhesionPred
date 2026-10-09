#!/bin/bash -l
# SignalP 6 (GPU, fast) for every proteome of kind e1 in RUN_LIST, in one job. Submit with -p short_gpu:
# sbatch -p short_gpu --export=ALL,PROJ_ROOT=...,RUN_LIST=... hydrophobin_signalp.sbatch
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH -J hyd_signalp
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${RUN_LIST:?export RUN_LIST}"
S="$PROJ_ROOT/scripts/sorting_hat"
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  [ "$kind" = e1 ] || continue
  echo "== $name $(date +%T)"
  FASTA="$fasta" WORKDIR="$workdir" bash -l "$S/signalp_gpu.sbatch"
done
echo "done $(date +%T)"
