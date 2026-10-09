#!/bin/bash -l
# Pre-freeze scoring of the TUNING data only: tuning proteomes, tuning hard-negative parts, and the T2 sequences (used only to
# choose between the two search options). hmmsearch of the seven-model template, default filters and --nobias, reporting every
# sequence with a score of at least 0 bits (--tblout). One job.
# Submit: sbatch --export=ALL,PROJ_ROOT=... prefreeze_scores.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH -J hyd_prefreeze
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
module load hmmer/3.4
H="$PROJ_ROOT/analysis/hydrophobin_truth"
OUT="$H/scores_pre_freeze"
mkdir -p "$OUT"
TUNING=$(awk -F'\t' 'NR>1 && $4=="tuning"{print $1}' "$H/proteome_split.tsv")
run() { # option flags, name, fasta
  local opt="$1" flags="$2" name="$3" faa="$4"
  hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$OUT/${name}.${opt}.tblout" "$H/relaxed/template.hmm" "$faa"
}
for opt in default nobias; do
  flags=""; [ "$opt" = nobias ] && flags="--nobias"
  for name in $TUNING; do
    faa=$(awk -F'\t' -v n="$name" 'NR>1 && $1==n{print $2}' "$H/run_list.tsv")
    run "$opt" "$flags" "$name" "$faa"
  done
  run "$opt" "$flags" hn_tuning "$OUT/hn_tuning.faa"
  run "$opt" "$flags" t2_positives "$OUT/t2_positives.faa"
done
echo done
