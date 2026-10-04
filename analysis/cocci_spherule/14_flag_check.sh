#!/bin/bash -l
# Re-test the antigen ranking's specificity flag on the follow-up genes with diamond.
# The flag (cocci_antigens/01_build_inputs.sh) used mmseqs easy-search against ONE proteome per
# confounder genus (first file by name) and human, keeping hits with identity >= 0.3 and
# coverage >= 0.5 of query and target. Here: the same proteomes, diamond very-sensitive,
# no identity filter in the search, so the thresholds can be applied afterwards.
# Usage: sbatch --export=ALL,FOLLOWUP=<dir> 14_flag_check.sh
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -J cocci_flagchk
set -euo pipefail
: "${FOLLOWUP:?}"
TMP="${SCRATCH:?SCRATCH is not set}/flagchk"
F5K=/bigdata/stajichlab/shared/projects/Fungi_5k/input
HUMAN=/bigdata/stajichlab/jstajich/projects/adhesionPred_review/cocci_antigens/human.faa
mkdir -p "$TMP" "$FOLLOWUP/results"
module load diamond/2.1.24
OUT="$FOLLOWUP/results/flag_check.tsv"
printf 'group\tproteome\tqseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\n' > "$OUT"
for pat in Histoplasma Blastomyces Paracoccidioides Aspergillus_fumigatus; do
  fa=$(ls "$F5K"/${pat}*.proteins.fa | head -1)
  diamond makedb --in "$fa" -d "$TMP/$pat" --quiet
  diamond blastp --very-sensitive -q "$FOLLOWUP/followup_proteins.fa" -d "$TMP/$pat" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --max-target-seqs 5 --evalue 1e-3 --threads 8 --quiet \
    | awk -v g="$pat" -v p="$(basename "$fa")" 'BEGIN{OFS="\t"}{print g,p,$0}' >> "$OUT"
done
diamond makedb --in "$HUMAN" -d "$TMP/human" --quiet
diamond blastp --very-sensitive -q "$FOLLOWUP/followup_proteins.fa" -d "$TMP/human" \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
  --max-target-seqs 5 --evalue 1e-3 --threads 8 --quiet \
  | awk 'BEGIN{OFS="\t"}{print "human","human.faa",$0}' >> "$OUT"
wc -l "$OUT"
