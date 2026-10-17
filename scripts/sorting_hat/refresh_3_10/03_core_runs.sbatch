#!/bin/bash -l
# Task 3.10, job 3 (DRAFT, not submitted): the classifier run (`cellsurface_sorting_hat`) on one proteome, with the
# new config and the rebuilt Pfam modules. An ARRAY job: task i runs the i-th proteome with run=yes in inputs.tsv
# (0 Afum_A1163, 1 Afum_Af293_Fungi5k, 2 Afum_Af293_UniProt, 3 Afum_W72310, 4 Bder_ER3, 5 Calb_SC5314,
# 6 Cimm_RS, 7 Scer_S288C, 8 Scer_S288C_nodubious): the nine runs of docs/reports/2026-10-07-classifier-runs-and-scale.md.
# It is submitted TWICE (each time needs the owner's go-ahead):
#   OUTNAME=out_pre    after jobs 0 to 2; the repeat call status is measured on this run (job 4);
#   OUTNAME=out_final  after job 4; the same run with every status file in place. Use it for the counts and statuses.
# Reads the new work directory; writes only $NEW/<proteome>/$OUTNAME (an existing directory stops the job).
# Config: the packaged categories.yaml of the pinned commit. Taxon and protein count per proteome: inputs.tsv.
# Submit: sbatch --array=0-8%3 --export=ALL,PROJ_ROOT=/path/to/clone,OUTNAME=out_pre \
#           -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log 03_core_runs.sbatch
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J csh_r310_core
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
: "${OUTNAME:?export OUTNAME=out_pre or out_final}"
: "${SLURM_ARRAY_TASK_ID:?submit with --array=0-8}"
case "$OUTNAME" in out_pre|out_final) ;; *) echo "FATAL: OUTNAME must be out_pre or out_final" >&2; exit 1 ;; esac
. "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/lib.sh"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_r310_core.$$"
trap 'rm -rf "$TMP" "$NEW"/*/.*.tmp.$$' EXIT
mkdir -p "$TMP"
check_pin
mapfile -t LIST < <(prots run)
p="${LIST[$SLURM_ARRAY_TASK_ID]:?no proteome for array index $SLURM_ARRAY_TASK_ID}"
echo "== $p $OUTNAME"
W="$NEW/$p"
[ ! -e "$W/$OUTNAME" ] || fail "$W/$OUTNAME exists; nothing is overwritten"
for m in pfam_adhesion pfam_allergen pfam_hydrophobin pfam_hsba; do [ -f "$W/modules/$m.json" ] || fail "$W/modules/$m.json missing; run job 2"; done
if [ "$OUTNAME" = out_final ] && { [ "$p" = Scer_S288C ] || [ "$p" = Calb_SC5314 ]; }; then
  [ -f "$W/status/calls/tandem_repeat_protein.json" ] || fail "$W/status/calls/tandem_repeat_protein.json missing; run job 4 first"
fi
check_fasta "$p"
make_venv "$TMP"
"$PY" -m cellsurface_sorting_hat.cli --fasta "$OLD/$p.faa" --taxon "$(field "$p" taxon)" --taxdump "$TAXDUMP/nodes.dmp" \
  --workdir "$W" --out "$TMP/out" --config "$CONFIG"
[ "$(grep -o '"n_proteins": [0-9]*' "$TMP/out/run.json" | grep -o '[0-9]*$')" = "$(field "$p" n_proteins)" ] || fail "n_proteins in run.json differs from inputs.tsv"
grep -o '"config_sha256": "[0-9a-f]*"\|"version": "[^"]*"\|"taxonomy_sha256": "[0-9a-f]*"' "$TMP/out/run.json"
cp -a "$TMP/out" "$W/.$OUTNAME.tmp.$$"
mv "$W/.$OUTNAME.tmp.$$" "$W/$OUTNAME"
echo "wrote $W/$OUTNAME"
