#!/bin/bash -l
# One step of the module pipeline for rows START..END of RUN_LIST, in one job (calls the existing job scripts inline).
# STEP is pfam, signalp, tmhmm, repeats or allergen. RUN_LIST columns: name fasta taxon kind workdir (header first).
# Submit: STEP=tmhmm START=1 END=10 PROJ_ROOT=... RUN_LIST=... sbatch --export=ALL scan_step.sbatch
# The signalp step needs: sbatch -p short_gpu --gres=gpu:1 --export=ALL scan_step.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:58:00
#SBATCH -J fscan_step
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${RUN_LIST:?export RUN_LIST}" "${STEP:?export STEP}" "${START:?export START}" "${END:?export END}"
S="$PROJ_ROOT/scripts/sorting_hat"
tail -n +2 "$RUN_LIST" | sed -n "${START},${END}p" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  echo "== $STEP $name $(date +%T)"
  case "$STEP" in
    pfam) FASTA="$fasta" WORKDIR="$workdir" FAMILY_TABLE="$PROJ_ROOT/data/sorting_hat/family_table.tsv" PFAM_RELEASE=38.2 bash -l "$S/pfam_hmmsearch.sbatch" ;;
    signalp) FASTA="$fasta" WORKDIR="$workdir" bash -l "$S/signalp_gpu.sbatch" ;;
    tmhmm) FASTA="$fasta" WORKDIR="$workdir" bash -l "$S/tmhmm.sbatch" ;;
    repeats) PROJ_ROOT="$PROJ_ROOT" FASTA="$fasta" WORKDIR="$workdir" bash -l "$S/repeats.sbatch" ;;
    allergen) FASTA="$fasta" WORKDIR="$workdir" ALLERGEN_FASTA="$PROJ_ROOT/_workdir/sorting_hat/iuis_fungal_allergens.faa" bash -l "$S/blast_allergen.sbatch" ;;
    *) echo "unknown STEP $STEP" >&2; exit 2 ;;
  esac
done
echo "done $(date +%T)"
