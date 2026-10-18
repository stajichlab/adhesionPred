#!/bin/bash -l
# Task 3.10, job 1 (run on HPCC on 2026-10-17; results: docs/reports/2026-10-18-task-3-10-refresh-results.md): `calibrate phasec` for step1_rule@R0 in the four work directories
# that had an R0 status file (r0=yes in inputs.tsv). The new file carries leakage_exception, the reworded
# note, negative_basis per stratum (N-int and N-sec curated, PM-TM assumed) and the per-stratum rates.
# `calibrate phasec` has no --negative-basis option: the curated basis is set inside the command (D16).
# Inputs: Phase C metrics.json (10 MB), eval_table.tsv.gz (1.5 MB), clusters.tsv.gz (1.5 MB), names.dmp (316 MB),
#   nodes.dmp (232 MB), data/sorting_hat/phasec_set_species.tsv, and modules/step1_rule@R0.json of each new work directory.
# The command runs on a scratch copy of the module records. The result is copied to
#   $NEW/<proteome>/status/step1_rule@R0.json (+ .sha256) under a temporary name and renamed. Needs job 0.
# Submit: sbatch --export=ALL,PROJ_ROOT=/path/to/clone -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log 01_phasec_r0.sbatch
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --time=0:45:00
#SBATCH -J csh_r310_phasec
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
. "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/lib.sh"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_r310_phasec.$$"
trap 'rm -rf "$TMP" "$NEW"/*/status/.*.tmp.$$' EXIT
mkdir -p "$TMP"
check_pin
make_venv "$TMP"
for f in "$PHASEC_DIR/metrics.json" "$PHASEC_DIR/eval_table.tsv.gz" "$PHASEC_DIR/clusters.tsv.gz" "$TAXDUMP/names.dmp" "$TAXDUMP/nodes.dmp"; do
  [ -s "$f" ] || fail "$f is missing or empty"
done
# build every result in scratch first; copy only when all four succeeded
for p in $(prots r0); do
  echo "== $p"
  [ -f "$NEW/$p/modules/step1_rule@R0.json" ] || fail "$NEW/$p/modules/step1_rule@R0.json not found; run job 0"
  [ ! -e "$NEW/$p/status/step1_rule@R0.json" ] || fail "$NEW/$p/status/step1_rule@R0.json exists; nothing is overwritten"
  mkdir -p "$TMP/$p/modules" "$TMP/$p/status"
  cp "$NEW/$p"/modules/*.json "$TMP/$p/modules/"
  "$PY" -m cellsurface_sorting_hat.calibration.cli phasec --workdir "$TMP/$p" \
    --metrics "$PHASEC_DIR/metrics.json" --clusters "$PHASEC_DIR/clusters.tsv.gz" --eval-table "$PHASEC_DIR/eval_table.tsv.gz" \
    --phasec-signalp-module signalp/6-gpu --phasec-signalp-mode fast \
    --set-species "$PROJ_ROOT/data/sorting_hat/phasec_set_species.tsv" \
    --names-dmp "$TAXDUMP/names.dmp" --nodes-dmp "$TAXDUMP/nodes.dmp"
  "$PY" "$R/refresh_checks.py" r0-check "$TMP/$p/status/step1_rule@R0.json"
  # for the log: how the new file differs from the old one (the old file stays where it is)
  diff <("$PY" -m json.tool "$OLD/$p/status/step1_rule@R0.json") <("$PY" -m json.tool "$TMP/$p/status/step1_rule@R0.json") | head -120 || true
done
for p in $(prots r0); do
  mkdir -p "$NEW/$p/status"
  cp "$TMP/$p/status/step1_rule@R0.json.sha256" "$NEW/$p/status/.step1_rule@R0.json.sha256.tmp.$$"
  mv "$NEW/$p/status/.step1_rule@R0.json.sha256.tmp.$$" "$NEW/$p/status/step1_rule@R0.json.sha256"
  cp "$TMP/$p/status/step1_rule@R0.json" "$NEW/$p/status/.step1_rule@R0.json.tmp.$$"
  mv "$NEW/$p/status/.step1_rule@R0.json.tmp.$$" "$NEW/$p/status/step1_rule@R0.json"
  echo "wrote $NEW/$p/status/step1_rule@R0.json"
done
