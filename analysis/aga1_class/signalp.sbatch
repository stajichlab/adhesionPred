#!/bin/bash -l
# AGA1-class prediction step, job 1 of 3: SignalP 6 (eukarya, fast mode, CPU build) on the 850 staged proteins.
# DRAFT, NOT SUBMITTED. Nothing of this step has run on HPCC. Needs the owner's go-ahead for EACH submission.
# The CPU build is used because the repo's GPU job refuses to run without a CUDA device and the partition that holds GPUs for
# `short` is NOT CONFIRMED; 850 proteins at the measured CPU rate (about 2 s/protein, scripts/sorting_hat/signalp_gpu.sbatch)
# is about 30 min on one core, less with threads. Heavy work is a SLURM job, never on the login node.
#
# Inputs:  $IN/aga1_proteins.faa (+ .sha256), staged by analysis/aga1_class/stage_inputs.py and copied by the owner.
# Output:  $OUT/prediction_results.txt, $OUT/provenance.txt (tool versions, input sha256, commit). $OUT must NOT exist.
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/aga1_class/signalp.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH -J aga1_signalp
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with aga1_proteins.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/aga1_signalp.$$"
trap 'rm -rf "$TMP"' EXIT
PINNED="analysis/aga1_class"
N_EXP=850
SIGNALP_MODULE="signalp/6"

fail() { echo "FATAL: $*" >&2; exit 1; }

# The code and the manifest must be the reviewed ones.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$PINNED" \
    || fail "$PINNED changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$PINNED")" ] || fail "uncommitted changes in $PINNED"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
FA="$IN/aga1_proteins.faa"
[ -f "$FA" ] && [ -f "$FA.sha256" ] || fail "$FA or its .sha256 is missing"
[ "$(sha256sum < "$FA" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$FA.sha256")" ] || fail "$FA differs from its .sha256"
[ "$(grep -c '^>' "$FA")" -eq "$N_EXP" ] || fail "$FA does not hold $N_EXP sequences"
grep '^>' "$FA" | sed 's/^>//; s/[[:space:]].*//' | sort | cmp -s - "$PROJ_ROOT/$PINNED/manifest_accessions.txt" \
  || fail "the names of $FA differ from the committed manifest"
mkdir -p "$TMP"

module load "$SIGNALP_MODULE"
command -v signalp6 >/dev/null || fail "signalp6 not on PATH after module load $SIGNALP_MODULE"
signalp6 --fastafile "$FA" --organism eukarya --output_dir "$TMP/out" \
  --format none --mode fast --write_procs 4 --torch_num_threads "${SLURM_CPUS_PER_TASK:-8}"
[ -s "$TMP/out/prediction_results.txt" ] || fail "no prediction_results.txt"
N_ROWS=$(grep -vc '^#' "$TMP/out/prediction_results.txt" || true)
[ "$N_ROWS" -eq "$N_EXP" ] || fail "prediction_results.txt has $N_ROWS rows, not $N_EXP"

{
  echo "job: aga1_signalp (AGA1-class prediction step)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "module: $SIGNALP_MODULE"
  echo "signalp6_version: $(signalp6 --version 2>&1)"
  echo "settings: --organism eukarya --mode fast --format none"
  echo "input sha256 $(sha256sum < "$FA" | cut -d' ' -f1)"
  echo "rows: $N_ROWS"
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP/out/prediction_results.txt" "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
