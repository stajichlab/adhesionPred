#!/bin/bash -l
# AGA1-class prediction step, job 3 of 3: per-residue disorder (AIUPred or IUPred3) on the 850 staged proteins.
# DRAFT, NOT SUBMITTED. Nothing of this step has run on HPCC. Needs the owner's go-ahead for EACH submission.
# PLACEHOLDER: the command lines of aiupred/2.x and iupred/3 on HPCC are NOT KNOWN to the author. Before submitting, check
# read-only on a login node:   module load aiupred; <aiupred script> --help     or     module load iupred/3; <iupred3 script> --help
# (also: `module avail aiupred` for the exact 2.x version, and what the tool prints per residue). Then edit the PLACEHOLDER-CMD
# line of the chosen tool and submit with TOOL=aiupred or TOOL=iupred3 and CMD_CONFIRMED=yes. Without CMD_CONFIRMED=yes the job
# stops after the pin check and the input checks, before any tool is run.
#
# Inputs:  $IN/aga1_proteins.faa (+ .sha256), staged by analysis/aga1_class/stage_inputs.py and copied by the owner.
# Output:  $OUT/disorder_raw.txt (the tool's own per-residue output, untouched), $OUT/provenance.txt. $OUT must NOT exist.
#          The converter summarises it (mean, fraction >= 0.5, longest run); it expects `>ID` records with `pos aa score` rows.
#
# Submit (after the go-ahead; the log directory must exist; TOOL is aiupred or iupred3):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,TOOL=aiupred,CMD_CONFIRMED=yes \
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
[ "${CMD_CONFIRMED:-no}" = "yes" ] || fail "PLACEHOLDER command not confirmed: check the CLI of $TOOL read-only, edit the PLACEHOLDER-CMD line, then export CMD_CONFIRMED=yes"
mkdir -p "$TMP"

case "$TOOL" in
  aiupred)
    TOOL_MODULE="aiupred"   # PLACEHOLDER: the exact module name/version (aiupred/2.x), see `module avail aiupred`
    module load "$TOOL_MODULE"
    # PLACEHOLDER-CMD (TO CONFIRM: module load aiupred; <script> --help): a GUESS of the form, not a verified call.
    aiupred.py -i "$FA" -o "$TMP/disorder_raw.txt"
    ;;
  iupred3)
    TOOL_MODULE="iupred/3"
    module load "$TOOL_MODULE"
    # PLACEHOLDER-CMD (TO CONFIRM: module load iupred/3; <script> --help): a GUESS of the form, not a verified call.
    # IUPred3 usually scores ONE sequence per call: if so, loop over the records here and write '>ID' + 'pos aa score' rows.
    iupred3.py "$FA" long > "$TMP/disorder_raw.txt"
    ;;
esac
[ -s "$TMP/disorder_raw.txt" ] || fail "$TOOL wrote no output"
N_REC=$(grep -c '^>' "$TMP/disorder_raw.txt" || true)
[ "$N_REC" -eq "$N_EXP" ] || fail "disorder_raw.txt has $N_REC '>' records, not $N_EXP (adapt the output layout to the converter's)"

{
  echo "job: aga1_disorder (AGA1-class prediction step)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "tool: $TOOL"
  echo "module: $TOOL_MODULE"
  echo "tool_version: PLACEHOLDER (record what the confirmed --help or a version flag prints)"
  echo "input sha256 $(sha256sum < "$FA" | cut -d' ' -f1)"
  echo "records: $N_REC"
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP/disorder_raw.txt" "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
