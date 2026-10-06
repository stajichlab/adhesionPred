#!/bin/bash -l
# Run both repeat detectors on one proteome FASTA. Output: $WORKDIR/raw/repeats/{repeat02,repeat14}.tsv
# Submit: sbatch --export=ALL,PROJ_ROOT=...,FASTA=...,WORKDIR=... repeats.sbatch
# Then, to make the module tables (--script is part of the module identity, so a changed detector gives a new identity):
#   cellsurface_sorting_hat_module repeat02 --table $WORKDIR/raw/repeats/repeat02.tsv --script $PROJ_ROOT/analysis/cocci_repeats/02_repeat_profile.py ...
#   cellsurface_sorting_hat_module repeat14 --table $WORKDIR/raw/repeats/repeat14.tsv --script $PROJ_ROOT/analysis/cocci_repeats/14_repeat_detect_general.py ...
# Run time is not measured for a whole proteome; record it from the job (sacct) in the run report.
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J csh_repeats
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_repeats.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/repeats"
mkdir -p "$TMP" "$OUT"
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
/usr/bin/python3.12 "$PROJ_ROOT/analysis/cocci_repeats/02_repeat_profile.py" "$TMP/in.fasta" --out "$TMP/repeat02.tsv"
/usr/bin/python3.12 "$PROJ_ROOT/analysis/cocci_repeats/14_repeat_detect_general.py" "$TMP/in.fasta" --out "$TMP/repeat14.tsv"
for m in repeat02 repeat14; do cp "$TMP/$m.tsv" "$OUT/.tmp.$m.tsv" && mv "$OUT/.tmp.$m.tsv" "$OUT/$m.tsv"; done
