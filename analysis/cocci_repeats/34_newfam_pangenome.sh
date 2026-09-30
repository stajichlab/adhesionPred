#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 4 --mem 32gb --time 1:30:00
#SBATCH -J newfampan -o logs/34_newfam_pangenome.%A.log
# Anchored search for the PTGIPTEWP repeat family (type protein CIMG_04070) across the
# 493-strain Coccidioides pangenome. Same sharding pattern as 31_anchor_pangenome.sh.
#
# This family was found by 32_class2a_anchor_sweep.py: secreted, 9 aa Pro/Thr repeat,
# one orthologue in each of the 7 long-read/reference proteomes, called by NEITHER
# periodicity detector. Question here: how prevalent is it, and does unit count vary?
#
# NOTE: no $(dirname "${BASH_SOURCE[0]}") -- wrong under sbatch. Use $SLURM_SUBMIT_DIR.
set -euo pipefail
DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$DIR"; mkdir -p logs
PAN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input
WORK="${SCRATCH:?SCRATCH is not set}/newfampan.$$"; mkdir -p "$WORK"
ls "$PAN"/*.proteins.fa > "$WORK/all.txt"
split -n l/8 -d "$WORK/all.txt" "$WORK/shard."
for S in "$WORK"/shard.*; do
  /usr/bin/python3.12 30_anchor_family_search.py search \
      --motif IPTEWP --reference newfam_CIMG_04070.faa --reference-id CIMG_04070 \
      --period 9 --proteomes $(cat "$S") --out "$S.tsv" > "$S.log" 2>&1 &
done
wait
head -1 "$WORK/shard.00.tsv" > newfam_anchored_pangenome.tsv
for S in "$WORK"/shard.*.tsv; do tail -n +2 "$S" >> newfam_anchored_pangenome.tsv; done
gzip -kf newfam_anchored_pangenome.tsv
wc -l newfam_anchored_pangenome.tsv
rm -rf "$WORK"
