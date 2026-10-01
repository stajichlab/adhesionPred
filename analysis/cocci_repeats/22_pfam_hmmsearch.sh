#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 16 --mem 32gb --time 2:00:00
#SBATCH -J pfamconf -o logs/22_pfam_hmmsearch.%A.log
# Pfam-A over the four protein sets built by 21_pfam_sets.py.
#
# The Pfam gathering threshold (--cut_ga) is used. It is the curated per-family cutoff that
# Pfam itself uses to define family membership, so no E-value is invented here.
#
# Also dumps NAME/ACC/LENG/DESC for every Pfam model, because the model length is what the
# detector's period estimate is checked against in 23_pfam_confirm.py.
#
# Submit from analysis/cocci_repeats:   sbatch 22_pfam_hmmsearch.sh
#
# NOTE: do not use $(dirname "${BASH_SOURCE[0]}") here. Under sbatch the script is a spooled
# copy and that path is wrong. $SLURM_SUBMIT_DIR is the submit directory.

set -euo pipefail

DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$DIR"
mkdir -p logs

module load hmmer/3.4

PFAM=/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm
WORK="${SCRATCH:?SCRATCH is not set}/pfamconf.$$"
mkdir -p "$WORK"

# node-local NVMe for the 1.9 GB profile database and the tabular output
cp "$PFAM" "$PFAM".h3f "$PFAM".h3i "$PFAM".h3m "$PFAM".h3p "$WORK"/
cp pfam_sets.fa "$WORK"/

hmmsearch --cpu "${SLURM_CPUS_PER_TASK:-8}" --cut_ga \
  --domtblout "$WORK/pfam_sets.domtbl" \
  --tblout "$WORK/pfam_sets.tbl" \
  -o /dev/null \
  "$WORK/Pfam-A.hmm" "$WORK/pfam_sets.fa"

# model metadata: accession, name, model length, description
awk '/^NAME /{n=$2} /^ACC /{a=$2} /^LENG /{l=$2}
     /^DESC /{d=""; for(i=2;i<=NF;i++) d=d (i>2?" ":"") $i}
     /^\/\/$/{print a "\t" n "\t" l "\t" d; a=n=l=d=""}' \
  "$WORK/Pfam-A.hmm" > "$WORK/pfam_models.tsv"

gzip -c "$WORK/pfam_sets.domtbl" > pfam_sets.domtbl.gz
gzip -c "$WORK/pfam_models.tsv" > pfam_models.tsv.gz
cp "$WORK/pfam_sets.domtbl" pfam_sets.domtbl
cp "$WORK/pfam_models.tsv" pfam_models.tsv

wc -l pfam_sets.domtbl pfam_models.tsv
rm -rf "$WORK"
