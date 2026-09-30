#!/bin/bash -l
# Issue #9 plan, step 1: does scaling matter, and how much does the pooling change move the
# shipped model's scores? Runs on node-local $SCRATCH, copies results to RESULTS_DIR.
#   sbatch analysis/model_review/step1_pooling_scaler.sbatch
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH -J adh_step1
#SBATCH -o /bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1/slurm.%j.log
set -euo pipefail
module load MMseqs2/17-b804f
PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
REPO=${ADH_REPO:?set ADH_REPO to the repo checkout to run}
RESULTS_DIR=${RESULTS_DIR:-/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1}
HERE=$REPO/analysis/model_review
W=${SCRATCH:?SCRATCH is not set}/step1
mkdir -p "$W" "$RESULTS_DIR"; cd "$W"
echo "repo=$REPO git=$(git -C "$REPO" rev-parse --short HEAD) host=$(hostname) start=$(date -Is)"

curl -sfL -o scer.fasta.gz \
  http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz \
  || cp "$RESULTS_DIR/scer.fasta.gz" scer.fasta.gz
gzip -dc "$REPO"/data/positive/*.gz "$REPO"/data/negative/*.gz > train_all.fa
mmseqs easy-cluster train_all.fa clu tmpmm --min-seq-id 0.3 -c 0.5 --cov-mode 0 -v 1
$PY "$HERE/01_embed_esm2_8M.py" scer.fasta.gz

for scaled in 1 0; do
  d=cv_scaler$scaled; mkdir -p $d
  ln -s "$W"/train_*.npy "$W"/train_meta.csv "$W"/scer_*.npy "$W"/scer_meta.csv "$W"/clu_cluster.tsv $d/
  (cd $d && ADH_SCALER=$scaled $PY "$HERE/02_cv_and_proteome_eval.py" | tee results.txt)
done
$PY "$HERE/03_pooling_mismatch.py"

cp scer.fasta.gz "$RESULTS_DIR"/
cp -rL cv_scaler1 cv_scaler0 pooling_mismatch_*.tsv train_meta.csv "$RESULTS_DIR"/
cp train_legacy.npy train_masked.npy scer_legacy.npy scer_masked.npy "$RESULTS_DIR"/
echo "done $(date -Is)"
