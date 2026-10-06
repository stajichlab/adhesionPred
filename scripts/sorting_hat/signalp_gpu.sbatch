#!/bin/bash -l
# SignalP 6 (GPU build, fast mode) on one proteome FASTA. Output: $WORKDIR/raw/signalp/prediction_results.txt
# Submit: sbatch --export=ALL,FASTA=/path/p.faa,WORKDIR=/path/work signalp_gpu.sbatch
# Then:   cellsurface_sorting_hat_module signalp --results $WORKDIR/raw/signalp/prediction_results.txt \
#           --signalp-version <contents of version.txt> --signalp-mode fast ...
# Measured rate on gpu12: about 192 proteins/s (analysis/step1_compare/jobs/j1_features.sh, job 29280458).
# The GPU build refuses to run without a CUDA device. Do not use the CPU build (about 2 s/protein).
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH -J csh_signalp
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_signalp.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$WORKDIR/raw/signalp"
module load signalp/6-gpu
command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH" >&2; exit 1; }
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
signalp6 --fastafile "$TMP/in.fasta" --organism eukarya --output_dir "$TMP/out" \
  --format none --mode fast --write_procs 4 --torch_num_threads 8
cp "$TMP/out/prediction_results.txt" "$WORKDIR/raw/signalp/.tmp.prediction_results.txt"
mv "$WORKDIR/raw/signalp/.tmp.prediction_results.txt" "$WORKDIR/raw/signalp/prediction_results.txt"
signalp6 --version > "$WORKDIR/raw/signalp/version.txt" 2>&1 || true
