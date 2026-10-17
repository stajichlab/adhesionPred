#!/bin/bash -l
# Task 3.10, job 4 (DRAFT, not submitted): `calibrate truth --call-status` for tandem_repeat_protein, as in M2d and M5b2
# (owner definition of 2026-10-08, truth_v3). An ARRAY job: task 0 = Scer_S288C (taxon 4932), task 1 = Calb_SC5314 (taxon 5476).
# Same inputs and parameters as the committed files docs/reports/data/sorting_hat/call_status/<proteome>/
# tandem_repeat_protein.owner_definition.json: the truth file truth_v3.<proteome>.tsv, the calibration set name
# repeat_motif_definition_<S288C|Calb>, --leakage tuned_on_truth, --negative-basis assumed (the negatives are
# assumed from an absent annotation, so 'curated' would be wrong), --n-boot 2000 --seed 1 (the defaults), and the
# same --notes text, which is read from the committed file.
# Needs: the run $NEW/<proteome>/out_pre (job 3), the module records of repeat02 and repeat14 in the new work directory.
# The command runs on a scratch work directory; the result is copied to $NEW/<proteome>/status/calls/ (temporary name, rename).
# The job compares the result with the numbers of the old file (S288C 0.435 and 0.987; C. albicans 0.462 and 0.976) and exits 3
# after the copy if they differ (stop and report, see the plan of 2026-10-08).
# Submit: sbatch --array=0-1 --export=ALL,PROJ_ROOT=/path/to/clone -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log 04_truth_repeat_call.sbatch
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH -J csh_r310_truth
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
: "${SLURM_ARRAY_TASK_ID:?submit with --array=0-1}"
. "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/lib.sh"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_r310_truth.$$"
trap 'rm -rf "$TMP" "$NEW"/*/status/calls/.*.tmp.$$' EXIT
mkdir -p "$TMP"
check_pin
case "$SLURM_ARRAY_TASK_ID" in
  0) p=Scer_S288C; taxon=4932; set_name=repeat_motif_definition_S288C; want_sens=0.435; want_spec=0.987 ;;
  1) p=Calb_SC5314; taxon=5476; set_name=repeat_motif_definition_Calb; want_sens=0.462; want_spec=0.976 ;;
  *) fail "array index must be 0 or 1" ;;
esac
echo "== $p"
W="$NEW/$p"
OLDFILE="$CALL_STATUS_DIR/$p/tandem_repeat_protein.owner_definition.json"
TRUTH="$TRUTH_DIR/truth_v3.$p.tsv"
[ -f "$OLDFILE" ] || fail "$OLDFILE not found"
[ -f "$TRUTH" ] || fail "$TRUTH not found"
for f in calls.long.tsv.gz run.json proteins.tsv.gz; do [ -f "$W/out_pre/$f" ] || fail "$W/out_pre/$f missing; run job 3 with OUTNAME=out_pre"; done
[ ! -e "$W/status/calls/tandem_repeat_protein.json" ] || fail "$W/status/calls/tandem_repeat_protein.json exists; nothing is overwritten"
make_venv "$TMP"
NOTES="$("$PY" "$R/refresh_checks.py" notes "$OLDFILE")"
echo "notes: $NOTES"
mkdir -p "$TMP/work/modules" "$TMP/work/status" "$TMP/out"
cp "$W"/modules/*.json "$TMP/work/modules/"
cp "$W/out_pre/calls.long.tsv.gz" "$W/out_pre/run.json" "$W/out_pre/proteins.tsv.gz" "$TMP/out/"
rc=0
"$PY" -m cellsurface_sorting_hat.calibration.cli truth --workdir "$TMP/work" --call-status --config "$CONFIG" \
  --calls-long "$TMP/out/calls.long.tsv.gz" --call tandem_repeat_protein --truth "$TRUTH" \
  --calibration-set "$set_name" --taxa "$taxon" --nodes-dmp "$TAXDUMP/nodes.dmp" \
  --notes "$NOTES" --leakage tuned_on_truth --negative-basis assumed
"$PY" "$R/refresh_checks.py" measure "$TMP/work/status/calls/tandem_repeat_protein.json" "$want_sens" "$want_spec" || rc=$?
mkdir -p "$W/status/calls"
for ext in json.sha256 json; do
  t="$W/status/calls/.tandem_repeat_protein.$ext.tmp.$$"
  cp "$TMP/work/status/calls/tandem_repeat_protein.$ext" "$t" && mv "$t" "$W/status/calls/tandem_repeat_protein.$ext"
done
echo "wrote $W/status/calls/tandem_repeat_protein.json (compare exit code $rc)"
exit "$rc"
