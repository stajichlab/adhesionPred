#!/bin/bash -l
# Map curated adhesins/hard negatives (data/curated.fa, UniProt) onto the pilot proteomes to build
# per-genome reference inventories for tier T5. Same-species matches only: >=90% identity and
# >=80% coverage of both query and target. Run from this directory:
#   srun -p epyc -c 8 --mem 16G -t 30 ./map_inventory.sh
set -euo pipefail
module load MMseqs2
mkdir -p inventory
: > inventory/hits.m8
while IFS=$'\t' read -r name fasta; do
  mmseqs easy-search data/curated.fa "$fasta" inventory/$name.m8 inventory/tmp_$name \
    --min-seq-id 0.9 -c 0.8 --cov-mode 0 --threads ${SLURM_CPUS_PER_TASK:-4} -v 1 \
    --format-output query,target,pident,qcov,tcov,evalue
  awk -v g="$name" 'BEGIN{OFS="\t"}{print g,$0}' inventory/$name.m8 >> inventory/hits.m8
  rm -rf inventory/tmp_$name
done < genomes.tsv
wc -l inventory/hits.m8
