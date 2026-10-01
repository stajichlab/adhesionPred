#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 4 --mem 32gb --time 1:30:00
#SBATCH -J anchpan -o logs/31_anchor_pangenome.%A.log
# Anchored SOWgp search over the 493-strain pangenome proteome set, in 8 shards run in
# parallel in one job. About 10 min of CPU in total.
#
# This answers the pangenome question in the 2026-09-30 anchored-search report: the 489
# SOWgp models in sowgp_pangenome.tsv came from the OrthoFinder ORTHOGROUPS. An anchored
# search over the proteomes themselves does not depend on orthogroup assignment, so it can
# recover members the orthogroup approach dropped, and it can tell split pairs from
# genuine fragments among the 284 models below 250 aa.
#
# Submit from analysis/cocci_repeats:   sbatch 31_anchor_pangenome.sh
#
# NOTE: no $(dirname "${BASH_SOURCE[0]}") -- wrong under sbatch. Use $SLURM_SUBMIT_DIR.

set -euo pipefail

DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$DIR"
mkdir -p logs

PAN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input
WORK="${SCRATCH:?SCRATCH is not set}/anchpan.$$"
mkdir -p "$WORK"

MOTIF="${ANCHOR_MOTIF:-KKYGDC}"
REF=sowgp_seed.fa
REFID=SOWgp58_Cocci_immitis_published

[ -f "$REF" ] || gunzip -kc "$REF.gz" > "$REF"

ls "$PAN"/*.proteins.fa > "$WORK/all.txt"
split -n l/8 -d "$WORK/all.txt" "$WORK/shard."

for S in "$WORK"/shard.*; do
  /usr/bin/python3.12 30_anchor_family_search.py search \
      --motif "$MOTIF" --reference "$REF" --reference-id "$REFID" --period 47 \
      --proteomes $(cat "$S") --out "$S.tsv" > "$S.log" 2>&1 &
done
wait
tail -n 3 "$WORK"/shard.*.log

head -1 "$WORK/shard.00.tsv" > sowgp_anchored_pangenome.tsv
for S in "$WORK"/shard.*.tsv; do tail -n +2 "$S" >> sowgp_anchored_pangenome.tsv; done

gzip -kf sowgp_anchored_pangenome.tsv
wc -l sowgp_anchored_pangenome.tsv
rm -rf "$WORK"
