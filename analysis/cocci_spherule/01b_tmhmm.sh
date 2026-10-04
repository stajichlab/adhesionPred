#!/bin/bash -l
# TMHMM 2.0c on the C. immitis RS RefSeq proteins (short format, one line per protein).
# Usage: sbatch 01b_tmhmm.sh <step1_workdir> <out_tsv>   (CPU job; node-local $SCRATCH)
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J cocci_tmhmm
set -euo pipefail
WORK="${1:?step1 workdir}"
OUT="${2:?output tsv}"
TMP="${SCRATCH:?SCRATCH is not set; run as a SLURM job}/cocci_tmhmm"
mkdir -p "$TMP"
module load tmhmm/2.0c
zcat "$WORK/phaseb/unique_sequences.fasta.gz" > "$TMP/all.fasta"
# TMHMM reads a FASTA; the short format prints: id len= ExpAA= First60= PredHel= Topology=
tmhmm -short "$TMP/all.fasta" > "$TMP/tmhmm.out"
awk 'BEGIN{OFS="\t"; print "protein_id","len","exp_aa","first60","pred_hel","topology"}
  {for(i=2;i<=NF;i++){sub(/^[^=]*=/,"",$i)} print $1,$2,$3,$4,$5,$6}' "$TMP/tmhmm.out" > "$OUT"
wc -l "$OUT"
