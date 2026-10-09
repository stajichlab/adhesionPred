#!/bin/bash -l
# SignalP 6 (GPU, fast) on the reference sequences (T1/T2/T3, LP, hard negatives), then the R0 module table. One job.
# Submit: sbatch -p short_gpu --export=ALL,PROJ_ROOT=... r0_reference.sbatch
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH -J hyd_r0_ref
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
W="$PROJ_ROOT/analysis/hydrophobin_truth/r0_reference"
FASTA="$W/reference.faa" WORKDIR="$W" bash -l "$PROJ_ROOT/scripts/sorting_hat/signalp_gpu.sbatch"
