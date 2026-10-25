#!/bin/bash -l
# AGA1-class proteome test: one SLURM array task = one chunk (at most 1000 proteins) of the staged UniProt proteomes, one tool.
# TOOL is signalp, predgpi or aiupred (same settings as the 850-protein run of 2026-10-24: SignalP 6 eukarya fast mode, PredGPI
# `-m gff3`, AIUPred 2.x defaults on a GPU). Heavy work is a SLURM job, never on the login node.
#
# Inputs:  $IN/chunks.tsv, $IN/chunks/<chunk>.faa, staged by analysis/aga1_class/stage_proteomes.py and copied by the owner.
# Output:  $OUT/<chunk>/{prediction_results.txt | predgpi_raw.gff3 | disorder_raw.txt} and provenance.txt. $OUT may already exist
#          (the array tasks share it); $OUT/<chunk> must NOT exist, nothing is overwritten.
#
# Defaults suit SignalP (-c 8). Override on the command line for the others (the GPU partition is not written here):
#   SignalP : sbatch --array=1-32%4 --export=ALL,PROJ_ROOT=..,CODE_COMMIT=..,IN=..,OUT=..,TOOL=signalp proteome_array.sbatch
#   PredGPI : sbatch -c 2 --mem=8G --time=2:00:00 --array=1-32%4 --export=...,TOOL=predgpi ...
#   AIUPred : sbatch -p <gpu-partition> --gres=gpu:1 -c 4 --mem=8G --time=1:30:00 --array=1-32%2 --export=...,TOOL=aiupred ...
#   (add -o $LOGS/%x.%A_%a.log -e $LOGS/%x.%A_%a.log; the log directory must exist)
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH -J aga1_prot
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the staged directory with chunks.tsv and chunks/) before sbatch}"
: "${OUT:?export OUT (the output base directory) before sbatch}"
: "${TOOL:?export TOOL=signalp, predgpi or aiupred}"
: "${SLURM_ARRAY_TASK_ID:?submit with --array}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/aga1_prot.${SLURM_ARRAY_JOB_ID:-x}.${SLURM_ARRAY_TASK_ID}.$$"
trap 'rm -rf "$TMP"' EXIT
PINNED="analysis/aga1_class"

fail() { echo "FATAL: $*" >&2; exit 1; }

check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$PINNED" \
    || fail "$PINNED changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$PINNED")" ] || fail "uncommitted changes in $PINNED"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

check_pin
case "$TOOL" in signalp|predgpi|aiupred) ;; *) fail "TOOL must be signalp, predgpi or aiupred, not $TOOL" ;; esac
[ -f "$IN/chunks.tsv" ] || fail "$IN/chunks.tsv is missing"
ROW=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR == n + 1 { print }' "$IN/chunks.tsv")
[ -n "$ROW" ] || fail "array task $SLURM_ARRAY_TASK_ID is beyond the $(($(wc -l < "$IN/chunks.tsv") - 1)) chunks"
CHUNK=$(echo "$ROW" | cut -f1); N_EXP=$(echo "$ROW" | cut -f3); SHA_EXP=$(echo "$ROW" | cut -f4)
FA="$IN/chunks/$CHUNK.faa"
[ -f "$FA" ] || fail "$FA is missing"
[ "$(sha256sum < "$FA" | cut -d' ' -f1)" = "$SHA_EXP" ] || fail "$FA differs from the sha256 in chunks.tsv"
[ "$(grep -c '^>' "$FA")" -eq "$N_EXP" ] || fail "$FA does not hold $N_EXP sequences"
D="$OUT/$CHUNK"
[ ! -e "$D" ] || fail "$D exists; nothing is overwritten"
mkdir -p "$TMP"
echo "task $SLURM_ARRAY_TASK_ID: $TOOL on $CHUNK ($N_EXP proteins)"

case "$TOOL" in
  signalp)
    TOOL_MODULE="signalp/6"
    module load "$TOOL_MODULE"
    command -v signalp6 >/dev/null || fail "signalp6 not on PATH after module load $TOOL_MODULE"
    signalp6 --fastafile "$FA" --organism eukarya --output_dir "$TMP/out" \
      --format none --mode fast --write_procs 4 --torch_num_threads "${SLURM_CPUS_PER_TASK:-8}"
    RAW="$TMP/out/prediction_results.txt"; RAWNAME=prediction_results.txt
    [ -s "$RAW" ] || fail "no prediction_results.txt"
    N_ROWS=$(grep -vc '^#' "$RAW" || true)
    [ "$N_ROWS" -eq "$N_EXP" ] || fail "prediction_results.txt has $N_ROWS rows, not $N_EXP"
    TOOL_INFO="signalp6_path: $(readlink -f "$(command -v signalp6)")
signalp6_version: $(signalp6 --version 2>&1)
settings: --organism eukarya --mode fast --format none"
    ;;
  predgpi)
    TOOL_MODULE="predgpi/202001"
    module load "$TOOL_MODULE"
    export PREDGPI_HOME=/opt/linux/rocky/8.x/x86_64/pkgs/predgpi/202001   # the module did not set it in a non-interactive shell
    [ -f "$PREDGPI_HOME/predgpi.py" ] || fail "$PREDGPI_HOME/predgpi.py is missing"
    predgpi.py -f "$FA" -o "$TMP/predgpi_raw.gff3" -m gff3
    RAW="$TMP/predgpi_raw.gff3"; RAWNAME=predgpi_raw.gff3
    [ -s "$RAW" ] || fail "PredGPI wrote no output"
    N_ROWS=$(wc -l < "$RAW")
    # one gff3 row per protein (verified on the 850-protein run); a different count is fatal here because every chunk is checked
    [ "$N_ROWS" -eq "$N_EXP" ] || fail "predgpi_raw.gff3 has $N_ROWS lines, not $N_EXP"
    TOOL_INFO="predgpi_py: $PREDGPI_HOME/predgpi.py
predgpi_py_sha256: $(sha256sum < "$PREDGPI_HOME/predgpi.py" | cut -d' ' -f1)
settings: -m gff3 (default omega site, no -c)"
    ;;
  aiupred)
    TOOL_MODULE="aiupred/2.x"
    module load "$TOOL_MODULE"
    AIUPRED_PY="$(command -v aiupred.py)" || fail "aiupred.py not on PATH after module load $TOOL_MODULE"
    python3 -c 'import torch, numpy, scipy' || fail "python3 cannot import torch, numpy and scipy"
    python3 -c 'import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)' \
      || fail "no CUDA device; submit to a GPU partition with --gres=gpu:1"
    python3 "$AIUPRED_PY" -i "$FA" -o "$TMP/disorder_raw.txt" -g 0
    RAW="$TMP/disorder_raw.txt"; RAWNAME=disorder_raw.txt
    [ -s "$RAW" ] || fail "AIUPred wrote no output"
    N_ROWS=$(grep -c '^#>' "$RAW" || true)   # AIUPred 2.x writes each record header as '#>ID'
    [ "$N_ROWS" -eq "$N_EXP" ] || fail "disorder_raw.txt has $N_ROWS '#>' records, not $N_EXP"
    TOOL_INFO="tool_path: $AIUPRED_PY
tool_script_sha256: $(sha256sum < "$AIUPRED_PY" | cut -d' ' -f1)
torch_version: $(python3 -c 'import torch; print(torch.__version__)' 2>&1 || true)
settings: defaults of aiupred.py, -g 0"
    ;;
esac

{
  echo "job: aga1_prot array task $SLURM_ARRAY_TASK_ID of ${SLURM_ARRAY_JOB_ID:-NA} (AGA1-class proteome test)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "tool: $TOOL"
  echo "module: $TOOL_MODULE"
  echo "module_list: $(module -t list 2>&1 | tr '\n' ' ')"
  echo "$TOOL_INFO"
  echo "chunk: $CHUNK ($N_EXP proteins)"
  echo "chunk sha256 $SHA_EXP"
  echo "output rows: $N_ROWS"
} > "$TMP/provenance.txt"

[ ! -e "$D" ] || fail "$D appeared while the job ran; nothing is overwritten"
mkdir -p "$D"
cp "$RAW" "$D/$RAWNAME"
cp "$TMP/provenance.txt" "$D/"
echo "wrote $D"
