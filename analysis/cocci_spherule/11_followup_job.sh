#!/bin/bash -l
# Follow-up checks for the specific extreme and Cys-rich spherule-up genes (CPU job).
#   1. phmmer of the follow-up proteins against the search proteomes (Onygenales + outgroups)
#   2. diamond (ultra-sensitive) against the pangenome RS proteome (gene-model check)
#   3. hmmscan of all spherule-up proteins against Pfam-A (gathering thresholds)
#   4. Phobius on all spherule-up proteins (second signal-peptide and TM call)
# Usage: sbatch --export=ALL,FOLLOWUP=<dir from 10_select_followup.py> 11_followup_job.sh
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH -J cocci_followup
set -euo pipefail
: "${FOLLOWUP:?export FOLLOWUP (output dir of 10_select_followup.py)}"
TMP="${SCRATCH:?SCRATCH is not set; run as a SLURM job}/cocci_followup"
PFAM=/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm
RSPAN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input/Coccidioides_immitis_RS.proteins.fa
mkdir -p "$TMP" "$FOLLOWUP/results"

echo "== 1. search database"
: > "$TMP/search.fa"
tail -n +2 "$FOLLOWUP/search_proteomes.tsv" | cut -f1 | while read -r f; do
  stem=$(basename "$f" .proteins.fa)
  awk -v s="$stem" '/^>/{print ">" s "|" substr($1,2); next} {print}' "$f" >> "$TMP/search.fa"
done
grep -c '^>' "$TMP/search.fa"

echo "== 1. phmmer"
module load hmmer/3.4
phmmer --cpu 16 -E 1e-3 --tblout "$FOLLOWUP/results/phmmer.tblout" \
  "$FOLLOWUP/followup_proteins.fa" "$TMP/search.fa" > "$TMP/phmmer.out"

echo "== 2. diamond against the pangenome RS proteome"
module load diamond/2.1.24
diamond makedb --in "$RSPAN" -d "$TMP/rs_pan" --quiet
diamond blastp --ultra-sensitive -q "$FOLLOWUP/followup_proteins.fa" -d "$TMP/rs_pan" \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \
  --max-target-seqs 5 --evalue 1e-3 --threads 16 --quiet -o "$FOLLOWUP/results/diamond_rs_pangenome.tsv"

echo "== 3. hmmscan Pfam-A"
hmmscan --cpu 16 --cut_ga --noali --domtblout "$FOLLOWUP/results/pfam.domtblout" \
  "$PFAM" "$FOLLOWUP/up_proteins.fa" > "$TMP/hmmscan.out"

echo "== 4. Phobius"
module load phobius/1.01
phobius.pl -short "$FOLLOWUP/up_proteins.fa" > "$FOLLOWUP/results/phobius_short.txt"
wc -l "$FOLLOWUP/results/"*
echo done
