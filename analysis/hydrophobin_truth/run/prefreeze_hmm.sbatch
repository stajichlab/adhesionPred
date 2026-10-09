#!/bin/bash -l
# Build leave-one-cluster-out HMMs and the all-cluster HMM (mafft L-INS-i, hmmbuild) for the folds in FOLDS (comma list of fold numbers, or 'all'), check the conserved cysteine columns,
# and score the TUNING data only (tuning proteomes and tuning hard-negative parts), default filters and --nobias. One job.
# Submit (commas in --export are separators, so pass FOLDS in the environment): FOLDS=1,2,3 PROJ_ROOT=... sbatch --export=ALL prefreeze_hmm.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:50:00
#SBATCH -J hyd_hmm_pre
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${FOLDS:?export FOLDS}"
module load hmmer/3.4 mafft/7.505
H="$PROJ_ROOT/analysis/hydrophobin_truth"
M="$H/hmm_models"
SC="$H/scores_pre_freeze_hmm"
mkdir -p "$M" "$SC"
/usr/bin/python3.12 "$H/lco.py" build --dir "$H" --out "$M" --threads 8 --only "$FOLDS"
TUNING=$(awk -F'\t' 'NR>1 && $4=="tuning"{print $1}' "$H/proteome_split.tsv")
for f in ${FOLDS//,/ }; do
  if [ "$f" = all ]; then base=fold_all; else base=$(printf "fold_%02d" "$f"); fi
  hmm="$M/$base.hmm"
  for opt in default nobias; do
    flags=""; [ "$opt" = nobias ] && flags="--nobias"
    for name in $TUNING; do
      faa=$(awk -F'\t' -v n="$name" 'NR>1 && $1==n{print $2}' "$H/run_list.tsv")
      hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$SC/${base}.${name}.${opt}.tblout" "$hmm" "$faa"
    done
    hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$SC/${base}.hn_tuning.${opt}.tblout" "$hmm" "$H/scores_pre_freeze/hn_tuning.faa"
  done
done
echo done
