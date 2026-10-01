#!/bin/bash -l
# J0 (spec 8): ESM-2 8M and 35M throughput on 2,000 truth sequences, then the spec 7
# GPU-CPU difference harness (which also repeats the GPU run to check repeatability).
# Submit (see the plan's run section):
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=... -o <log> -e <log> j0_pilot.sh
# PROJ_ROOT comes from the environment, never from the script location (that is wrong in
# SLURM spool directories). Work happens in node-local $SCRATCH; results are copied to
# $STEP1_WORKDIR/phaseb/j0/ before the job ends.
# exfab: one node (gpu12, 2x ada6000), usually less contended than short_gpu. The time limit
# is 2 h: the GPU rate is not measured yet (spec 8 assumes under 15 min; an estimate is 5 min
# at 29,500 residues/s and 48 min at 3,000 residues/s). J0 measures it.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH -J step1_j0
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_j0"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/jobs"
OUT="$STEP1_WORKDIR/phaseb/j0"
mkdir -p "$TMP" "$OUT"

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv > "$TMP/nvidia_smi.csv"
"$ENV_PY" "$S1/jobs/throughput_pilot.py" --work-dir "$STEP1_WORKDIR" \
  --out "$TMP/throughput.json" --device cuda --n 2000 --batch-sizes 8,16,32,64
"$ENV_PY" "$S1/jobs/gpu_cpu_diff.py" --work-dir "$STEP1_WORKDIR" \
  --out "$TMP/gpu_cpu_diff.json" --device-a cuda --device-b cpu --n 200 --n-long 20

cp "$TMP/throughput.json" "$TMP/gpu_cpu_diff.json" "$TMP/nvidia_smi.csv" "$OUT/"
echo "J0 done: $OUT"
