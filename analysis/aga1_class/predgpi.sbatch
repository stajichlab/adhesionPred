#!/bin/bash -l
# AGA1-class prediction step, job 2 of 3: PredGPI on the 850 staged proteins.
# DRAFT, NOT SUBMITTED. Nothing of this step has run on HPCC. Needs the owner's go-ahead for EACH submission.
# PLACEHOLDER: the PredGPI command line on HPCC is NOT KNOWN to the author. Before submitting, check read-only on a login node:
#     module load predgpi/202001; echo "$PREDGPI_HOME"; ls "$PREDGPI_HOME"; <the main script> --help
# (module predgpi/202001 loads miniconda3, sets PREDGPI_HOME=/opt/linux/rocky/8.x/x86_64/pkgs/predgpi/202001 and puts that
# directory on PATH.) Then replace the line marked PLACEHOLDER-CMD with the confirmed command and submit with CMD_CONFIRMED=yes.
# Without CMD_CONFIRMED=yes the job stops after the pin check and the input checks, before any tool is run.
#
# Inputs:  $IN/aga1_proteins.faa (+ .sha256), staged by analysis/aga1_class/stage_inputs.py and copied by the owner.
# Output:  $OUT/predgpi_raw.txt (the tool's own output, untouched), $OUT/provenance.txt. $OUT must NOT exist.
#          The converter wants the NORMALISED table (see convert_predictions.py); adapting the raw output to it is done
#          after the first real output has been read (README, "to confirm").
#
# Submit (after the go-ahead; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,CMD_CONFIRMED=yes \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/aga1_class/predgpi.sbatch
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J aga1_predgpi
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with aga1_proteins.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/aga1_predgpi.$$"
trap 'rm -rf "$TMP"' EXIT
PINNED="analysis/aga1_class"
N_EXP=850
TOOL_MODULE="predgpi/202001"

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
[ "${CMD_CONFIRMED:-no}" = "yes" ] || fail "PLACEHOLDER command not confirmed: check the PredGPI CLI read-only, edit the PLACEHOLDER-CMD line, then export CMD_CONFIRMED=yes"
mkdir -p "$TMP"

module load "$TOOL_MODULE"
: "${PREDGPI_HOME:?PREDGPI_HOME not set by module $TOOL_MODULE}"
# PLACEHOLDER-CMD (TO CONFIRM: module load predgpi/202001; <script> --help): the next line is a GUESS of the form, not a verified call.
predgpi.py -f "$FA" -o "$TMP/predgpi_raw.txt"
[ -s "$TMP/predgpi_raw.txt" ] || fail "PredGPI wrote no output"

{
  echo "job: aga1_predgpi (AGA1-class prediction step)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "module: $TOOL_MODULE"
  echo "predgpi_home: $PREDGPI_HOME"
  echo "predgpi_dir_listing: $(ls "$PREDGPI_HOME" | tr '\n' ' ')"
  echo "tool_version: PLACEHOLDER (record what the confirmed --help or a version flag prints)"
  echo "input sha256 $(sha256sum < "$FA" | cut -d' ' -f1)"
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP/predgpi_raw.txt" "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
