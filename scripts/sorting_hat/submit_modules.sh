#!/bin/bash
# Submit the module jobs for one proteome. Usage:
#   PROJ_ROOT=... FASTA=... WORKDIR=... FAMILY_TABLE=... PFAM_RELEASE=38.2 ALLERGEN_FASTA=... bash submit_modules.sh
# Set DRY_RUN=1 to print the sbatch commands without submitting.
# Prints one job ID per line. Wait for the jobs (squeue / sacct), then run the converters
# (module commands of Task 13) and the core command. All paths come from PROJ_ROOT.
# Relative paths are made absolute here: the sbatch jobs get absolute paths only (tmhmm.sbatch changes directory).
# The allergen converter must pass --max-target-seqs 200 --evalue 1 --seg no (the blastp options in
# blast_allergen.sbatch); the cap must exceed the database size so no hit is lost.
# The repeat converters must pass --script <detector path> (see repeats.sbatch).
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
: "${FASTA:?export FASTA}"
: "${WORKDIR:?export WORKDIR}"
: "${FAMILY_TABLE:?export FAMILY_TABLE}"
: "${PFAM_RELEASE:?export PFAM_RELEASE}"
: "${ALLERGEN_FASTA:?export ALLERGEN_FASTA}"
PROJ_ROOT="$(realpath "$PROJ_ROOT")"
FASTA="$(realpath "$FASTA")"
WORKDIR="$(realpath -m "$WORKDIR")"
FAMILY_TABLE="$(realpath "$FAMILY_TABLE")"
ALLERGEN_FASTA="$(realpath "$ALLERGEN_FASTA")"
export PROJ_ROOT FASTA WORKDIR FAMILY_TABLE PFAM_RELEASE ALLERGEN_FASTA
S="$PROJ_ROOT/scripts/sorting_hat"
[ "${DRY_RUN:-0}" = 1 ] || mkdir -p "$WORKDIR/logs"
for job in signalp_gpu pfam_hmmsearch repeats blast_allergen tmhmm; do
  if [ "${DRY_RUN:-0}" = 1 ]; then
    echo "sbatch --parsable --export=ALL -o $WORKDIR/logs/$job.%j.log -e $WORKDIR/logs/$job.%j.log $S/$job.sbatch"
  else
    sbatch --parsable --export=ALL -o "$WORKDIR/logs/$job.%j.log" -e "$WORKDIR/logs/$job.%j.log" "$S/$job.sbatch"
  fi
done
