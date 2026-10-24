#!/bin/bash -l
# AGA1-class prediction step, job 3 of 3: per-residue disorder (AIUPred or IUPred3) on the 850 staged proteins.
# DRAFT, NOT SUBMITTED. Nothing of this step has run on HPCC. Needs the owner's go-ahead for EACH submission.
# Command lines confirmed read-only from the install directories (2026-10-24):
#   aiupred/2.x : `aiupred.py -i FASTA [-o OUT] [-g GPU] [--force-cpu]`; the module prepends .../pkgs/aiupred/2.x/AIUPred to PATH
#                 and puts NO python on PATH; it needs torch (~=2.0.1), numpy, scipy; CPU is documented as very slow.
#   iupred/3    : `iupred3.py SEQFILE {long|short|glob} [-a] [-s no|medium|strong]` (needs scipy). It reads ONE sequence per
#                 file, so this job calls it once per record.
# STILL UNVERIFIED: torch/scipy availability after the module load (the job checks and stops with a clear message), the
# CPU/GPU speed, and the raw output layouts (AIUPred: '>header' then 'pos<TAB>aa<TAB>score' per README; IUPred3: '#' comment
# lines then 'pos aa score' rows, from memory). A wrong record count is REPORTED (warning and provenance), not fatal, because
# the raw file is the product; the converter does the strict checks.
#
# GPU: AIUPred should run on a GPU. The partition that holds GPUs is NOT written here: pass it on the sbatch command line
# (`-p <gpu-partition> --gres=gpu:1`, see the README). On a node without CUDA the AIUPred job stops unless ALLOW_CPU=yes
# is exported (CPU AIUPred is very slow; the 1:30 limit may then not be enough).
#
# Inputs:  $IN/aga1_proteins.faa (+ .sha256), staged by analysis/aga1_class/stage_inputs.py and copied by the owner.
# Output:  $OUT/disorder_raw.txt (per-residue, '>ID' records), $OUT/provenance.txt. $OUT must NOT exist.
#
# Submit (after the go-ahead; the log directory must exist; TOOL is aiupred or iupred3):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,TOOL=aiupred \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/aga1_class/disorder.sbatch
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH -J aga1_disorder
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with aga1_proteins.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${TOOL:?export TOOL=aiupred or TOOL=iupred3}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/aga1_disorder.$$"
trap 'rm -rf "$TMP"' EXIT
PINNED="analysis/aga1_class"
N_EXP=850

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
case "$TOOL" in aiupred|iupred3) ;; *) fail "TOOL must be aiupred or iupred3, not $TOOL" ;; esac
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
FA="$IN/aga1_proteins.faa"
[ -f "$FA" ] && [ -f "$FA.sha256" ] || fail "$FA or its .sha256 is missing"
[ "$(sha256sum < "$FA" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$FA.sha256")" ] || fail "$FA differs from its .sha256"
[ "$(grep -c '^>' "$FA")" -eq "$N_EXP" ] || fail "$FA does not hold $N_EXP sequences"
grep '^>' "$FA" | sed 's/^>//; s/[[:space:]].*//' | sort | cmp -s - "$PROJ_ROOT/$PINNED/manifest_accessions.txt" \
  || fail "the names of $FA differ from the committed manifest"
mkdir -p "$TMP"

case "$TOOL" in
  aiupred)
    TOOL_MODULE="aiupred/2.x"
    module load "$TOOL_MODULE"
    AIUPRED_PY="$(command -v aiupred.py)" || fail "aiupred.py not on PATH after module load $TOOL_MODULE"
    python3 -c 'import torch, numpy, scipy' || fail "python3 cannot import torch, numpy and scipy (AIUPred needs them; no python module is loaded by $TOOL_MODULE)"
    if python3 -c 'import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)'; then
      DEVICE="gpu"
      python3 "$AIUPRED_PY" -i "$FA" -o "$TMP/disorder_raw.txt" -g 0
    else
      [ "${ALLOW_CPU:-no}" = "yes" ] || fail "no CUDA device; AIUPred on CPU is very slow. Submit to a GPU partition, or export ALLOW_CPU=yes"
      DEVICE="cpu"
      python3 "$AIUPRED_PY" -i "$FA" -o "$TMP/disorder_raw.txt" --force-cpu
    fi
    TOOL_PATH="$AIUPRED_PY"
    ;;
  iupred3)
    TOOL_MODULE="iupred/3"
    module load "$TOOL_MODULE"
    IUPRED_PY="$(command -v iupred3.py)" || IUPRED_PY=/opt/linux/rocky/8.x/x86_64/pkgs/iupred/3/iupred3.py
    [ -f "$IUPRED_PY" ] || fail "iupred3.py not found ($IUPRED_PY)"
    python3 -c 'import scipy' || fail "python3 cannot import scipy (IUPred3 needs it)"
    DEVICE="cpu"
    # IUPred3 reads one sequence per file: split the staged FASTA by record (the FASTA has one name per '>' line) and
    # write '>ID' followed by the tool's 'pos aa score' rows ('#' comment lines of the tool are dropped).
    mkdir "$TMP/one"
    : > "$TMP/disorder_raw.txt"
    awk -v d="$TMP/one" '/^>/ { name = substr($1, 2); f = d "/" name ".faa" } { print >> f; close(f) }' "$FA"
    while read -r acc; do
      echo ">$acc" >> "$TMP/disorder_raw.txt"
      python3 "$IUPRED_PY" "$TMP/one/$acc.faa" long | grep -v '^#' >> "$TMP/disorder_raw.txt"
    done < "$PROJ_ROOT/$PINNED/manifest_accessions.txt"
    TOOL_PATH="$IUPRED_PY"
    ;;
esac
[ -s "$TMP/disorder_raw.txt" ] || fail "$TOOL wrote no output"
N_REC=$(grep -c '^>' "$TMP/disorder_raw.txt" || true)
[ "$N_REC" -eq "$N_EXP" ] || echo "WARNING: disorder_raw.txt has $N_REC '>' records, not $N_EXP; read it before conversion" >&2

{
  echo "job: aga1_disorder (AGA1-class prediction step)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "tool: $TOOL"
  echo "module: $TOOL_MODULE"
  echo "tool_path: $TOOL_PATH"
  echo "tool_script_sha256: $(sha256sum < "$TOOL_PATH" | cut -d' ' -f1)"
  echo "device: $DEVICE"
  echo "torch_version: $(python3 -c 'import torch; print(torch.__version__)' 2>&1 || true)"
  echo "settings: $([ "$TOOL" = iupred3 ] && echo 'long mode' || echo 'defaults of aiupred.py')"
  echo "input sha256 $(sha256sum < "$FA" | cut -d' ' -f1)"
  echo "records: $N_REC of $N_EXP"
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP/disorder_raw.txt" "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
