#!/bin/bash -l
# J2 (spec 8, D3): ESM-2 8M and 35M, layer 6, residue-mean pooling, N- and C-terminal windows.
# Submit (see the plan's run section) with the job count and time from job_plan.json:
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=...,J2_JOB_COUNT=<n_jobs> \
#          --array=0-<n_jobs-1> --time=<time_minutes> -o <log> -e <log> j2_embed.sh
# Array task i embeds chunks k with k mod J2_JOB_COUNT == i (06_plan_embedding.py). Each
# finished chunk is written to node-local $SCRATCH and copied at once to
# $STEP1_WORKDIR/phaseb/emb/<model>/, so a killed task loses at most the chunk in progress
# and a resubmission skips the finished chunks after a hash check.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH -J step1_j2
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
JOB_INDEX="${SLURM_ARRAY_TASK_ID:-0}"
JOB_COUNT="${J2_JOB_COUNT:-1}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_j2_${JOB_INDEX}"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/jobs"
LOGDIR="$STEP1_WORKDIR/phaseb/emb/logs"
mkdir -p "$TMP" "$LOGDIR"

nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used --format=csv,noheader -l 30 \
  > "$TMP/gpu_util.csv" &
SMI=$!
set +e
"$ENV_PY" "$S1/jobs/embed_chunks.py" --work-dir "$STEP1_WORKDIR" \
  --job-index "$JOB_INDEX" --job-count "$JOB_COUNT" --device cuda --scratch-dir "$TMP/emb"
STATUS=$?
set -e
kill "$SMI" || true
gzip -n -c "$TMP/gpu_util.csv" > "$LOGDIR/gpu_util.${SLURM_JOB_ID:-local}_${JOB_INDEX}.csv.gz"
echo "J2 task $JOB_INDEX of $JOB_COUNT: exit $STATUS"
exit "$STATUS"
