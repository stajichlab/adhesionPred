#!/usr/bin/bash
#SBATCH -p epyc -N 1 -n 1 -c 32 --mem 64gb --time 3:00:00
#SBATCH -J asmref -o logs/50_asm_vs_ref.%A.log
# Align every short-read Coccidioides assembly to a chromosome-level long-read reference and
# record, per 10 kb bin of the reference, how much of it the assembly contains.
#
# Question: how much sequence does a short-read assembly lack compared with a long-read one?
# References: C. posadasii -> Silveira 2022 (Nanopore, 9 contigs); C. immitis -> CiB10637
# (canu PacBio, 14 contigs). Each species is aligned to a reference of its own species.
#
# Needs asm_vs_ref_strains.tsv (built from Assembly/asm_stats.tsv). Output: asm_vs_ref_bins.tsv.gz
# NOTE: no $(dirname "${BASH_SOURCE[0]}") -- wrong under sbatch. Use $SLURM_SUBMIT_DIR.
set -euo pipefail
DIR="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$DIR"; mkdir -p logs
module load minimap2/2.30

G=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Genotyping/genome
U=/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/2026/annotation/assemblies/canu
POS=$G/FungiDB-68_CposadasiiSilveira2022_Genome.fasta
IMM=$U/CiB10637.sorted.fasta
BIN=10000

WORK="${SCRATCH:?SCRATCH is not set}/asmref.$$"; mkdir -p "$WORK/out"
lens () { /usr/bin/python3.12 - "$1" <<'PY'
import sys
n=None;c=0
for l in open(sys.argv[1]):
    if l.startswith(">"):
        if n: print(n,c)
        n=l[1:].split()[0]; c=0
    else: c+=len(l.strip())
print(n,c)
PY
}
lens "$POS" > "$WORK/pos.lens"; lens "$IMM" > "$WORK/imm.lens"
cp "$POS" "$WORK/pos.fa"; cp "$IMM" "$WORK/imm.fa"

run_one () {
  strain=$1; species=$2; fasta=$3
  if [ "$species" = posadasii ]; then ref=$WORK/pos.fa; lens=$WORK/pos.lens; rn=Silveira2022
  else ref=$WORK/imm.fa; lens=$WORK/imm.lens; rn=CiB10637; fi
  minimap2 -x asm10 --secondary=no -t 1 "$ref" "$fasta" 2>/dev/null \
    | /usr/bin/python3.12 51_paf_bins.py "$strain" "$species" "$rn" "$BIN" "$lens" > "$WORK/out/$strain.tsv"
}
export -f run_one; export WORK BIN

tail -n +2 asm_vs_ref_strains.tsv | awk -F'\t' '{print $1"\t"$2"\t"$3}' \
  | xargs -P 30 -L 1 bash -c 'run_one "$0" "$1" "$2"'

( echo -e "strain\tspecies\tref\tcontig\tbin_start\tbin_end\tcovered\tdepth"
  cat "$WORK"/out/*.tsv ) | gzip > asm_vs_ref_bins.tsv.gz
ls "$WORK/out" | wc -l
rm -rf "$WORK"
