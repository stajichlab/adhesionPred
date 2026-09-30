#!/usr/bin/bash -l
# Align SOWgp protein sets with MAFFT (module mafft).
# Inputs:
#   sowgp_seed.fa         8 sequences: 5 genome copies + published SOWgp58/66/82 (05_sowgp_pangenome.py --stage seed)
#   sowgp_pangenome.fa    489 annotated SOWgp copies from the 496-genome pangenome (05_sowgp_pangenome.py --stage extract)
# Outputs:
#   sowgp_seed.msa.fa
#   sowgp_pangenome.msa.fa
# Run: sbatch -N 1 -c 8 --mem 8G --time 2:0:0 -p short 06_sowgp_msa.sh
# or interactively: bash -lc 'module load mafft && bash 06_sowgp_msa.sh'

module load mafft
set -euo pipefail
cd "$(dirname "$0")"

mafft --localpair --thread 8 sowgp_seed.fa > sowgp_seed.msa.fa
mafft --localpair --thread 8 sowgp_pangenome.fa > sowgp_pangenome.msa.fa

echo "aligned sowgp_seed.msa.fa"
echo "aligned sowgp_pangenome.msa.fa"
