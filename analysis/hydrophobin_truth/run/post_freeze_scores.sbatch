#!/bin/bash -l
# POST-freeze scoring for L5: refuses to run unless the freeze commit is an ancestor of HEAD and freeze.json is unchanged.
# relaxed template (frozen option --nobias, scores at or above 0; the cutoff is applied offline) on the test proteomes, the test
# hard-negative parts, the scored LP proteins and T3; strict Pfam GA (--cut_ga) on T2, T3, LP and the test hard negatives;
# per-fold HMMs on their held-out members; the all-cluster HMM on the test hard negatives, LP and the test proteomes.
# Submit: PROJ_ROOT=... sbatch --export=ALL post_freeze_scores.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH -J hyd_post_freeze
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
module load hmmer/3.4
H="$PROJ_ROOT/analysis/hydrophobin_truth"
cd "$PROJ_ROOT"
FREEZE=$(git log -1 --format=%H -- analysis/hydrophobin_truth/freeze.json)
git merge-base --is-ancestor "$FREEZE" HEAD || { echo "FATAL: the freeze commit is not an ancestor of HEAD" >&2; exit 1; }
git diff --quiet "$FREEZE" -- analysis/hydrophobin_truth/freeze.json || { echo "FATAL: freeze.json changed after the freeze" >&2; exit 1; }
OUT="$H/scores_post_freeze"
F="$OUT/fasta"
mkdir -p "$OUT"
TESTP=$(awk -F'\t' 'NR>1 && $4=="test"{print $1}' "$H/proteome_split.tsv")
T="$H/relaxed/template.hmm"
rel() { hmmsearch --nobias -T 0 --cpu 8 --noali -o /dev/null --tblout "$OUT/relaxed.$1.tblout" "$T" "$2"; }
ga() { hmmsearch --cut_ga --cpu 8 --noali -o /dev/null --tblout "$OUT/gagate.$1.tblout" "$T" "$2"; }
for name in $TESTP; do
  faa=$(awk -F'\t' -v n="$name" 'NR>1 && $1==n{print $2}' "$H/run_list.tsv")
  rel "$name" "$faa"
done
for set in hn_test lp_scored t3; do rel "$set" "$F/$set.faa"; done
for set in t2 t3 lp_scored hn_test; do ga "$set" "$F/$set.faa"; done
M="$H/hmm_models"
for i in $(seq -w 1 28); do
  for opt in default nobias; do
    flags=""; [ "$opt" = nobias ] && flags="--nobias"
    hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$OUT/hmm.fold_$i.test.$opt.tblout" "$M/fold_$i.hmm" "$F/fold_${i}_test.faa"
  done
done
for opt in default nobias; do
  flags=""; [ "$opt" = nobias ] && flags="--nobias"
  for set in hn_test lp_scored; do
    hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$OUT/hmm.fold_all.$set.$opt.tblout" "$M/fold_all.hmm" "$F/$set.faa"
  done
  for name in $TESTP; do
    faa=$(awk -F'\t' -v n="$name" 'NR>1 && $1==n{print $2}' "$H/run_list.tsv")
    hmmsearch $flags -T 0 --cpu 8 --noali -o /dev/null --tblout "$OUT/hmm.fold_all.$name.$opt.tblout" "$M/fold_all.hmm" "$faa"
  done
done
echo done
