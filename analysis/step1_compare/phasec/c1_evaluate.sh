#!/bin/bash -l
# Phase C job C1 (spec 5): 09 (clusters and splits), 10 (fit and score), 11 (metrics) on CPU.
# Submit (see the plan's run section):
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=... -o <log> -e <log> c1_evaluate.sh
# One job, not one job per split: the three steps depend on each other, and the global rule
# says not to split work into jobs of minutes. No run time is measured on the real data yet;
# the first submission is the pilot. It prints `C1 step <n> wall_seconds=<s>` for every step,
# and the run section sizes later submissions from these lines (C1_STEPS selects steps).
# The time limit of 4 h is a bound for the pilot, not a measurement (estimate in the plan).
# Timings from the Task 8 and Task 9 reviews (single-process runs on this cluster; review
# estimates, NOT measurements of C1): step 09 25 min 20 s with 2 CPUs on c01; one fit_unit per
# candidate 61 to 81 s at 5,600 rows on 1 BLAS thread; step 11 about 41 s per 1,000 test rows
# at 2,000 resamples on 1 process (about 18 min for the real run, an estimate).
# Memory (estimate, not measured): step 10 runs WORKERS worker processes and each holds a full
# copy of the universe, so memory is (workers + 1) x universe size. Embeddings alone are
# 69,941 x (320 + 480) float32 = about 0.22 GB; with the feature table and the long-sequence
# arrays a reviewer measured about 1.2 GB peak RSS for step 11 at 7,150 rows. Estimate for
# step 10: (16 + 1) x 1.2 GB = about 20 GB. --mem=48G keeps a margin of more than 2 x.
# A job that runs an AVX2 tool must request a node feature that has AVX2 (--constraint=ryzen on
# partition epyc). Partition epyc and constraint ryzen (ruling C-15): the epyc nodes carry the
# features ryzen, amd, milan; their CPUs have AVX2, which the MMseqs2 module build needs (it
# stops with exit 132 on the abu_dhabi Opterons of partition batch). Step 11 runs one process
# with one BLAS thread, so its wall time is a one-core number, not a 16-core number.
# Step 09 gets --threads from the SLURM allocation. Every step stops with exit code 2 on a
# STOP; set -e passes that code on, so no step masks a failure of an earlier one.
# PROJ_ROOT comes from the environment, never from the script location. Temp files go to
# node-local $SCRATCH; the scripts write their outputs to $STEP1_WORKDIR/phasec/ (on /bigdata)
# with temp names and os.replace, so the results are on /bigdata when the job ends.
#SBATCH -p epyc
#SBATCH --constraint=ryzen
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH -J step1_c1
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_c1"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/phasec"
# one BLAS thread per process: 10 runs one process per CPU
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CPUS="${SLURM_CPUS_PER_TASK:-4}"
STEPS="${C1_STEPS:-09 10 11}"
LOGDIR="$STEP1_WORKDIR/phasec/logs"
mkdir -p "$TMP" "$LOGDIR"
: > "$TMP/wall.txt"
WALL_OUT="$LOGDIR/wall.${SLURM_JOB_ID:-local}.txt"
# copy the wall times on every exit path, also after a failed step (review M-1)
trap 'cp "$TMP/wall.txt" "$WALL_OUT"' EXIT

timed() {  # timed NAME CMD...: run CMD and record its wall time
  local name=$1 t0=$SECONDS
  shift
  "$@"
  echo "C1 step $name wall_seconds=$((SECONDS - t0))" | tee -a "$TMP/wall.txt"
}

for step in $STEPS; do
  case "$step" in
    09) timed 09 bash "$S1/phasec/09_cluster_and_split.sh" ;;
    10) timed 10 "$ENV_PY" "$S1/phasec/10_fit_and_score.py" --work-dir "$STEP1_WORKDIR" \
          --workers "$CPUS" --candidates "${PHASEC_CANDIDATES:-B0,B1,R0,R1,R2,M8,M35,M8-C,M35-C,H}" ;;
    11) timed 11 "$ENV_PY" "$S1/phasec/11_evaluate.py" --work-dir "$STEP1_WORKDIR" \
          --sets "${PHASEC_SETS:-$S1/sequence_sets.tsv}" \
          --species "${PHASEC_SPECIES:-$S1/species.tsv}" \
          --n-resamples "${PHASEC_N_RESAMPLES:-2000}" ;;
    *) echo "C1: unknown step $step" >&2; exit 2 ;;
  esac
done
echo "C1 done: $STEPS"
