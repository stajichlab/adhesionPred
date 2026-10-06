#!/bin/bash -l
# TMHMM 2.0c on one proteome FASTA. Output: $WORKDIR/raw/tmhmm/tmhmm.tsv (columns as analysis/cocci_spherule/01b_tmhmm.sh)
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=... tmhmm.sbatch
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH -J csh_tmhmm
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_tmhmm.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/tmhmm"
mkdir -p "$TMP" "$OUT"
module load tmhmm/2.0c
cd "$TMP"  # TMHMM writes a TMHMM_<pid> directory into the current directory
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
tmhmm -short "$TMP/in.fasta" > "$TMP/tmhmm.out"
awk 'BEGIN{OFS="\t"; print "protein_id","len","exp_aa","first60","pred_hel","topology"}
  {for(i=2;i<=NF;i++){sub(/^[^=]*=/,"",$i)} print $1,$2,$3,$4,$5,$6}' "$TMP/tmhmm.out" > "$TMP/tmhmm.tsv"
cp "$TMP/tmhmm.tsv" "$OUT/.tmp.tmhmm.tsv" && mv "$OUT/.tmp.tmhmm.tsv" "$OUT/tmhmm.tsv"
