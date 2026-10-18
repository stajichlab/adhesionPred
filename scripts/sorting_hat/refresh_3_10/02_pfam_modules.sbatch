#!/bin/bash -l
# Task 3.10, job 2 (run on HPCC on 2026-10-17; results: docs/reports/2026-10-18-task-3-10-refresh-results.md): rebuild the four Pfam modules (pfam_adhesion, pfam_hydrophobin,
# pfam_hsba, pfam_allergen) for the 14 proteomes of inputs.tsv, so that their identity (artefact hash with
# the family digest, params hash with narrow_families) matches the new family table.
# It does NOT run hmmsearch. It converts the domain tables that already exist (raw/pfam/domtbl.txt, made on
# 2026-10-08 with all 24 models of the family table; the model set of the table has not changed since commit
# 3819f00 of 2026-10-08 09:00). To search again, use scripts/sorting_hat/pfam_hmmsearch.sbatch (not part of this job).
# Checks before anything is copied:
#   - every accession in each domain table is a family of the table, and the table ends with '# [ok]' and shows --cut_ga;
#   - for the 8 proteomes with a committed pfam_hydrophobin status file, the rebuilt pfam_hydrophobin identity
#     (module, version, params_hash, artefact_hash) equals the one in that file, so the file stays valid.
#   - revision after the run (issue #95): the model list in raw/pfam/provenance.json (written by pfam_hmmsearch.sbatch) must contain
#     every row of the family table (accession and name). Provenance without a list (the 14 directories of 2026-10-08) fails
#     unless ALLOW_UNRECORDED_MODELS=1 is exported; then a WARNING line is logged for each proteome and the job goes on.
#   If one check fails nothing is copied to $NEW. The old work directories are never touched.
# Output: $NEW/<proteome>/modules/pfam_{adhesion,allergen,hydrophobin,hsba}.{json,tsv.gz} (+ .sha256). Needs job 0.
# Submit: sbatch --export=ALL,PROJ_ROOT=/path/to/clone -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log 02_pfam_modules.sbatch
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH -J csh_r310_pfam
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
. "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/lib.sh"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_r310_pfam.$$"
trap 'rm -rf "$TMP" "$NEW"/*/modules/.*.tmp.$$' EXIT
mkdir -p "$TMP"
check_pin
make_venv "$TMP"
module load hmmer/3.4   # only to read the version for the module record; hmmsearch is not run
HMMER_VERSION="$(hmmsearch -h | sed -n 2p | awk '{print $3}')"
[ "$HMMER_VERSION" = 3.4 ] || fail "hmmer module is $HMMER_VERSION; the domain tables were made with 3.4"
tail -n +2 "$FAMILY_TABLE" | cut -f1 | sort -u > "$TMP/table_acc.txt"
MODULES="pfam_adhesion pfam_allergen pfam_hydrophobin pfam_hsba"
# 1. build in scratch
for p in $(tail -n +2 "$INPUTS" | cut -f1); do
  echo "== $p"
  W="$NEW/$p"
  [ -f "$W/raw/pfam/domtbl.txt" ] || fail "$W/raw/pfam/domtbl.txt not found; run job 0"
  for m in $MODULES; do [ ! -e "$W/modules/$m.json" ] || fail "$W/modules/$m.json exists; nothing is overwritten"; done
  D="$W/raw/pfam/domtbl.txt"
  tail -n 3 "$D" | grep -qx '# \[ok\]' || fail "$D has no '# [ok]' trailer"
  grep -q -e '--cut_ga' "$D" || fail "$D does not show --cut_ga"
  grep -v '^#' "$D" | awk '{split($5,a,"."); print a[1]}' | sort -u > "$TMP/seen_acc.txt"
  [ -s "$TMP/seen_acc.txt" ] || fail "$D has no hit"
  stray=$(comm -13 "$TMP/table_acc.txt" "$TMP/seen_acc.txt" | tr '\n' ' ')
  [ -z "$stray" ] || fail "$D has hits of accessions that are not in the family table: $stray"
  # every model of the family table must have been searched (a model never searched looks like a clean zero)
  rc=0
  "$PY" "$R/refresh_checks.py" models-searched "$W/raw/pfam/provenance.json" "$FAMILY_TABLE" || rc=$?
  case "$rc" in
    0) ;;
    4) if [ "${ALLOW_UNRECORDED_MODELS:-}" = 1 ]; then
         echo "WARNING: $p: $W/raw/pfam/provenance.json has no list of the searched models; ALLOW_UNRECORDED_MODELS=1, so it is NOT known that every model of the family table was searched"
       else
         fail "$p: $W/raw/pfam/provenance.json has no list of the searched models. Export ALLOW_UNRECORDED_MODELS=1 (owner decision) to go on, or search again with pfam_hmmsearch.sbatch"
       fi ;;
    *) fail "$p: a model of the family table is not in the model list of $W/raw/pfam/provenance.json" ;;
  esac
  check_fasta "$p"
  mkdir -p "$TMP/$p/modules"
  cp "$W"/modules/* "$TMP/$p/modules/"
  "$PY" -m cellsurface_sorting_hat.modules.cli pfam --fasta "$OLD/$p.faa" --workdir "$TMP/$p" --domtbl "$D" \
    --family-table "$FAMILY_TABLE" --pfam-release "$PFAM_RELEASE" --pfam-sha256 "$PFAM_HMM_SHA256" \
    --hmmer-version "$HMMER_VERSION" --sp-module 'step1_rule@R0' --tm-module tm
  for m in $MODULES; do
    [ -f "$TMP/$p/modules/$m.json" ] || fail "$p: $m.json was not written"
    # old identity (from the old work directory) beside the new one, for the log
    echo "$m old: $(grep -h '"artefact_hash"\|"params_hash"' "$OLD/$p/modules/$m.json" | tr -d ' \n')"
    echo "$m new: $(grep -h '"artefact_hash"\|"params_hash"' "$TMP/$p/modules/$m.json" | tr -d ' \n')"
  done
done
# 2. the hydrophobin status files must stay valid
for p in $(prots hyd); do
  "$PY" "$R/refresh_checks.py" identity-equal "$TMP/$p/modules/pfam_hydrophobin.json" "$CALL_STATUS_DIR/hydrophobin/pfam_hydrophobin.$p.json" \
    || fail "$p: the rebuilt pfam_hydrophobin has another identity than the committed status file; the file would be stale. Nothing was copied."
done
# 3. copy (temporary name, then rename)
for p in $(tail -n +2 "$INPUTS" | cut -f1); do
  for m in $MODULES; do
    for ext in json.sha256 tsv.gz.sha256 json tsv.gz; do
      t="$NEW/$p/modules/.$m.$ext.tmp.$$"
      cp "$TMP/$p/modules/$m.$ext" "$t" && mv "$t" "$NEW/$p/modules/$m.$ext"
    done
  done
  echo "wrote four Pfam modules in $NEW/$p/modules"
done
