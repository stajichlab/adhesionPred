#!/bin/bash -l
# AGA1-class proteome test: one SLURM array task = one chunk (at most 1000 proteins) of the staged UniProt proteomes, one tool.
# TOOL is signalp, predgpi, aiupred, pfam_family, pfam_full or tmhmm (the first three with the settings of the 850-protein run of
# 2026-10-24: SignalP 6 eukarya fast mode, PredGPI `-m gff3`, AIUPred 2.x defaults on a GPU). pfam_family searches only the models of
# data/sorting_hat/family_table.tsv (hmmsearch --cut_ga, what the sorting-hat `pfam` module reads), pfam_full searches all of Pfam-A
# (--cut_ga) for domain architecture with coordinates, tmhmm runs TMHMM 2.0c (-short). The two Pfam tools need FAMILY_TABLE (pfam_family only),
# PFAM_RELEASE and optionally PFAM_HMM exported. Heavy work is a SLURM job, never on the login node.
#
# Inputs:  $IN/chunks.tsv, $IN/chunks/<chunk>.faa, staged by analysis/aga1_class/stage_proteomes.py and copied by the owner.
# Output:  $OUT/<chunk>/{prediction_results.txt | predgpi_raw.gff3 | disorder_raw.txt} and provenance.txt. $OUT may already exist
#          (the array tasks share it); $OUT/<chunk> must NOT exist, nothing is overwritten.
#
# Defaults suit SignalP (-c 8). Override on the command line for the others (the GPU partition is not written here):
#   SignalP : sbatch --array=1-32%4 --export=ALL,PROJ_ROOT=..,CODE_COMMIT=..,IN=..,OUT=..,TOOL=signalp proteome_array.sbatch
#   PredGPI : sbatch -c 2 --mem=8G --time=2:00:00 --array=1-32%4 --export=...,TOOL=predgpi ...
#   AIUPred : sbatch -p <gpu-partition> --gres=gpu:1 -c 4 --mem=8G --time=1:30:00 --array=1-32%2 --export=...,TOOL=aiupred ...
#   Pfam    : sbatch -c 8 --mem=8G --time=2:00:00 --array=1-32%4 --export=...,TOOL=pfam_full,PFAM_RELEASE=38.2 ...   (pfam_family: add FAMILY_TABLE=...)
#   TMHMM   : sbatch -c 2 --mem=4G --time=2:00:00 --array=1-32%6 --export=...,TOOL=tmhmm ...
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
case "$TOOL" in signalp|predgpi|aiupred|pfam_family|pfam_full|tmhmm) ;; *) fail "TOOL must be signalp, predgpi, aiupred, pfam_family, pfam_full or tmhmm, not $TOOL" ;; esac
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
  pfam_family|pfam_full)
    : "${PFAM_RELEASE:?export PFAM_RELEASE (for example 38.2) for the Pfam tools}"
    PFAM_HMM="${PFAM_HMM:-/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm}"
    TOOL_MODULE="hmmer/3.4"
    module load "$TOOL_MODULE"
    RESOLVED="$(readlink -f "$PFAM_HMM")"
    case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
    if [ "$TOOL" = pfam_family ]; then
      : "${FAMILY_TABLE:?export FAMILY_TABLE (data/sorting_hat/family_table.tsv) for pfam_family}"
      tail -n +2 "$FAMILY_TABLE" | cut -f2 | sort -u > "$TMP/keys.txt"
      hmmfetch -f "$PFAM_HMM" "$TMP/keys.txt" > "$TMP/models.hmm"
      N_MODELS=$(grep -c '^NAME' "$TMP/models.hmm" || true)
      [ "$N_MODELS" -eq "$(wc -l < "$TMP/keys.txt")" ] || fail "hmmfetch found $N_MODELS of $(wc -l < "$TMP/keys.txt") models"
      awk '/^NAME/{n=$2} /^ACC/{split($2,a,"."); print a[1]"\t"n}' "$TMP/models.hmm" | sort > "$TMP/fetched.txt"
      tail -n +2 "$FAMILY_TABLE" | cut -f1,2 | sort > "$TMP/table.txt"
      diff -q "$TMP/fetched.txt" "$TMP/table.txt" > /dev/null || fail "model names in the database differ from the family table"
      MODELS="$TMP/models.hmm"
      SETTINGS="hmmsearch --cut_ga with the $N_MODELS models of the family table (sha256 of the table: $(sha256sum < "$FAMILY_TABLE" | cut -d' ' -f1))"
    else
      MODELS="$PFAM_HMM"
      SETTINGS="hmmsearch --cut_ga with all of Pfam-A"
    fi
    hmmsearch --cut_ga --cpu "${SLURM_CPUS_PER_TASK:-2}" --noali -o /dev/null --domtblout "$TMP/domtbl.txt" "$MODELS" "$FA"
    tail -n 3 "$TMP/domtbl.txt" | grep -qx '# \[ok\]' || fail "domain table has no '# [ok]' trailer"
    RAW="$TMP/domtbl.txt"; RAWNAME=domtbl.txt
    N_ROWS=$(grep -vc '^#' "$RAW" || true)
    TOOL_INFO="pfam_hmm: $PFAM_HMM
pfam_resolved: $RESOLVED
pfam_release: $PFAM_RELEASE
pfam_sha256: $(sha256sum < "$RESOLVED" | cut -d' ' -f1)
hmmer: $(hmmsearch -h | sed -n 2p)
settings: $SETTINGS"
    ;;
  tmhmm)
    TOOL_MODULE="tmhmm/2.0c"
    module load "$TOOL_MODULE"
    command -v tmhmm >/dev/null || fail "tmhmm not on PATH after module load $TOOL_MODULE"
    ( cd "$TMP" && tmhmm -short "$FA" > "$TMP/tmhmm.out" )   # TMHMM writes a TMHMM_<pid> directory into the current directory
    awk 'BEGIN{OFS="\t"; print "protein_id","len","exp_aa","first60","pred_hel","topology"}
      {for(i=2;i<=NF;i++){sub(/^[^=]*=/,"",$i)} print $1,$2,$3,$4,$5,$6}' "$TMP/tmhmm.out" > "$TMP/tmhmm.tsv"
    RAW="$TMP/tmhmm.tsv"; RAWNAME=tmhmm.tsv
    N_ROWS=$(($(wc -l < "$RAW") - 1))
    [ "$N_ROWS" -eq "$N_EXP" ] || fail "tmhmm.tsv has $N_ROWS rows, not $N_EXP"
    TOOL_INFO="settings: tmhmm -short (TMHMM 2.0c)"
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
