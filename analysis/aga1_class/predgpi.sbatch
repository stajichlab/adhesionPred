#!/bin/bash -l
# AGA1-class prediction step, job 2 of 3: PredGPI on the 850 staged proteins.
# DRAFT, NOT SUBMITTED. Nothing of this step has run on HPCC. Needs the owner's go-ahead for EACH submission.
# Command line confirmed read-only from the install directory (2026-10-24): `predgpi.py -f FASTA -o OUTF -m json|gff3`
# (-m default gff3; -c selects the conservative omega site). `module load predgpi/202001` did NOT set PREDGPI_HOME in a
# non-interactive shell, so this job exports it and stops if the directory or predgpi.py is missing.
# STILL UNVERIFIED: the raw gff3 layout and its key names (the source prints 'Estimated False Positive Rate =',
# 'HMM log(probability) =', 'GPI-Anchor length=' and one of Highly Probable / Probable / Weakly GPI-Anchored protein /
# Improbable GPI-Anchored protein), and the run time. The raw file is therefore kept untouched; the converter wants the
# NORMALISED table (convert_predictions.py), so an adapter is written after the first real output is read.
#
# Inputs:  $IN/aga1_proteins.faa (+ .sha256), staged by analysis/aga1_class/stage_inputs.py and copied by the owner.
# Output:  $OUT/predgpi_raw.gff3 (the tool's own output, untouched), $OUT/provenance.txt. $OUT must NOT exist.
#
# Submit (after the go-ahead; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/aga1_class/predgpi.sbatch
#SBATCH -p short
#SBATCH -c 2
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
PREDGPI_DIR=/opt/linux/rocky/8.x/x86_64/pkgs/predgpi/202001

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

module load "$TOOL_MODULE"
export PREDGPI_HOME="$PREDGPI_DIR"   # the module did not set it in a non-interactive shell
[ -d "$PREDGPI_HOME" ] || fail "$PREDGPI_HOME is missing"
[ -f "$PREDGPI_HOME/predgpi.py" ] || fail "$PREDGPI_HOME/predgpi.py is missing"
PREDGPI_PY="$PREDGPI_HOME/predgpi.py"
predgpi.py -f "$FA" -o "$TMP/predgpi_raw.gff3" -m gff3
[ -s "$TMP/predgpi_raw.gff3" ] || fail "PredGPI wrote no output"
N_LINES=$(wc -l < "$TMP/predgpi_raw.gff3")
N_IDS=0
while read -r a; do grep -qF -- "$a" "$TMP/predgpi_raw.gff3" && N_IDS=$((N_IDS + 1)); done < "$PROJ_ROOT/$PINNED/manifest_accessions.txt" || true
# The raw layout is UNVERIFIED: a short count is reported, not fatal (the raw file is the product; the adapter checks it).
[ "$N_IDS" -eq "$N_EXP" ] || echo "WARNING: only $N_IDS of $N_EXP accessions appear in the raw output; read it before conversion" >&2

{
  echo "job: aga1_predgpi (AGA1-class prediction step)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "module: $TOOL_MODULE"
  echo "predgpi_home: $PREDGPI_HOME"
  echo "predgpi_py: $PREDGPI_PY"
  echo "predgpi_py_sha256: $(sha256sum < "$PREDGPI_PY" | cut -d' ' -f1)"
  echo "predgpi_dir_listing: $(ls "$PREDGPI_HOME" | tr '\n' ' ')"
  echo "settings: -m gff3 (default omega site, no -c)"
  echo "input sha256 $(sha256sum < "$FA" | cut -d' ' -f1)"
  echo "raw_lines: $N_LINES; manifest accessions found in raw output: $N_IDS of $N_EXP"
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP/predgpi_raw.gff3" "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
