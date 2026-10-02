#!/bin/bash -l
# Phase C step 2 (spec 5): MMseqs2 clustering, maximum-identity search and the split tables.
# c1_evaluate.sh runs this script inside its SLURM job; it is not submitted alone.
# PROJ_ROOT and STEP1_WORKDIR come from the environment, never from the script location.
# MMseqs2 temp files go to node-local $SCRATCH; 09_make_splits.py writes its outputs to
# $STEP1_WORKDIR/phasec/ (on /bigdata) with a temp name and os.replace.
# The module MMseqs2/17-b804f puts the AVX2 build on PATH. It stops with "Illegal instruction"
# (exit 132) on CPUs without AVX2 (for example the abu_dhabi nodes of partition batch), so
# c1_evaluate.sh requests partition epyc with --constraint=ryzen. STEP1_MMSEQS (a binary path)
# skips the module (tests pass the stub there). A job that runs an AVX2 tool must request a node
# feature that has AVX2; the two #SBATCH lines below keep this script on such nodes if it is
# ever submitted alone (c1_evaluate.sh requests the same).
#SBATCH -p epyc
#SBATCH --constraint=ryzen
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root)}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR}"
TMP="${SCRATCH:?SCRATCH is not set; run this inside a SLURM job}/phasec_09"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/phasec"
if [ -z "${STEP1_MMSEQS:-}" ]; then
  module load "${PHASEC_MMSEQS_MODULE:-MMseqs2/17-b804f}"
  STEP1_MMSEQS=$(command -v mmseqs)
fi
"$ENV_PY" "$S1/phasec/09_make_splits.py" --work-dir "$STEP1_WORKDIR" \
  --species "${PHASEC_SPECIES:-$S1/species.tsv}" --mmseqs "$STEP1_MMSEQS" \
  --tmp-dir "$TMP" --threads "${SLURM_CPUS_PER_TASK:-4}"
