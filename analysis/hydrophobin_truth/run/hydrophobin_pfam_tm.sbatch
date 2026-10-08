#!/bin/bash -l
# Pfam hmmsearch (full family table) for every proteome in RUN_LIST, and TMHMM for the proteomes of kind e1.
# One job for all proteomes (jobs of 1 hour or more are better than one job per proteome).
# Submit: sbatch --export=ALL,PROJ_ROOT=...,RUN_LIST=...,FAMILY_TABLE=...,PFAM_RELEASE=38.2 hydrophobin_pfam_tm.sbatch
# RUN_LIST columns: name fasta taxon kind workdir (tab separated, with header). Calls the existing job scripts inline.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J hyd_pfam_tm
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${RUN_LIST:?export RUN_LIST}" "${FAMILY_TABLE:?export FAMILY_TABLE}" "${PFAM_RELEASE:?export PFAM_RELEASE}"
S="$PROJ_ROOT/scripts/sorting_hat"
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  echo "== $name $(date +%T)"
  FASTA="$fasta" WORKDIR="$workdir" FAMILY_TABLE="$FAMILY_TABLE" PFAM_RELEASE="$PFAM_RELEASE" bash -l "$S/pfam_hmmsearch.sbatch"
  if [ "$kind" = e1 ]; then
    FASTA="$fasta" WORKDIR="$workdir" bash -l "$S/tmhmm.sbatch"
  fi
done
echo "done $(date +%T)"
