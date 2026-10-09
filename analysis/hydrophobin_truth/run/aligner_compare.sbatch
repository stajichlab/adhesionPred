#!/bin/bash -l
# Exploratory aligner comparison (after L5; does not change the frozen v1). For one aligner (ALIGNER=mafft|famsa|muscle5):
# build the 28 leave-one-cluster-out HMMs and the all-cluster HMM, score the tuning data (for cutoffs), the held-out members,
# the test hard negatives and the scored LP proteins. One job per aligner.
# Submit: ALIGNER=famsa PROJ_ROOT=... sbatch --export=ALL aligner_compare.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:55:00
#SBATCH -J hyd_aligner
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${ALIGNER:?export ALIGNER}"
module load hmmer/3.4 mafft/7.505 famsa/2.4.1 muscle/5.1
H="$PROJ_ROOT/analysis/hydrophobin_truth"
M="$H/aligner_compare/hmm_$ALIGNER"
SC="$H/aligner_compare/scores_$ALIGNER"
F="$H/scores_post_freeze/fasta"
mkdir -p "$M" "$SC"
/usr/bin/python3.12 "$H/lco.py" build --dir "$H" --out "$M" --threads 8 --aligner "$ALIGNER"
TUNING=$(awk -F'\t' 'NR>1 && $4=="tuning"{print $1}' "$H/proteome_split.tsv")
for i in $(seq -w 1 28) all; do
  base="fold_$i"
  for opt in default nobias; do
    flags=""; [ "$opt" = nobias ] && flags="--nobias"
    for name in $TUNING; do
      faa=$(awk -F'\t' -v n="$name" 'NR>1 && $1==n{print $2}' "$H/run_list.tsv")
      hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$SC/${base}.${name}.${opt}.tblout" "$M/$base.hmm" "$faa"
    done
    hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$SC/${base}.hn_tuning.${opt}.tblout" "$M/$base.hmm" "$H/scores_pre_freeze/hn_tuning.faa"
    if [ "$i" = all ]; then
      for set in hn_test lp_scored; do
        hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$SC/hmm.fold_all.$set.$opt.tblout" "$M/fold_all.hmm" "$F/$set.faa"
      done
    else
      hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$SC/hmm.fold_$i.test.$opt.tblout" "$M/$base.hmm" "$F/fold_${i}_test.faa"
    fi
  done
done
echo done
