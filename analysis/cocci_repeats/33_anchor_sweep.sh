#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 2 --mem 32gb --time 1:30:00
#SBATCH -J anchswp -o logs/33_anchor_sweep.%A.log
# Anchored search applied to every class-2a candidate family (32_class2a_anchor_sweep.py).
# One job, single core for the background motif count, about 25 min.
#
# Submit from analysis/cocci_repeats:   sbatch 33_anchor_sweep.sh
#
# NOTE: no $(dirname "${BASH_SOURCE[0]}") -- wrong under sbatch. Use $SLURM_SUBMIT_DIR.

set -euo pipefail

DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$DIR"
mkdir -p logs

for F in class2a_candidates_general.tsv class2a_candidates.tsv \
         repeat_profile_longread.tsv repeat_profile_reference.tsv \
         repeat_general_longread.tsv repeat_general_reference.tsv; do
  [ -f "$F" ] || gunzip -kc "$F.gz" > "$F"
done

/usr/bin/python3.12 32_class2a_anchor_sweep.py \
    --per-family 12 \
    --out-families class2a_anchor_families.tsv \
    --out-hits class2a_anchor_hits.tsv

gzip -kf class2a_anchor_families.tsv class2a_anchor_hits.tsv
