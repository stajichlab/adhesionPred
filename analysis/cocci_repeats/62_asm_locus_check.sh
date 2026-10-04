#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 32 --mem 32gb --time 2:00:00
#SBATCH -J asmlocus -o logs/62_asm_locus_check.%A.log
# Is the RS CIMG_04613 (SOWgp) locus present in each strain's short-read ASSEMBLY?
#
# Query : RS GG704914:968,094-973,690 (CIMG_04613 gene +/- 2 kb, 5,597 bp).
# Target: three sequence sets per strain
#   asm  = Assembly/asm/AAFTF/<strain>.sorted.fasta                     (AAFTF contigs)
#   mito = Assembly/asm/AAFTF/<strain>.vecscreen.mitochondria.fasta (contigs the AAFTF
#          mitochondrial screen removed from the nuclear assembly)
#   ann  = Assembly/annotation/<strain>/annotate_results/*.scaffolds.fa  (the genome the
#          pangenome gene models were predicted on)
# Tools : minimap2 -x asm20 (assembly contigs mapped onto the query), and blastn
#         (query vs assembly, -subject) as a second method.
# Output: <strain>.<asm|ann|mito>.paf, .blastn.tsv, .mito.len, .ann.nrun.tsv packed into asm_locus_raw.tar.gz;
#         parse with 63_asm_locus_parse.py.
# Run from analysis/cocci_repeats:  sbatch 62_asm_locus_check.sh
# No $(dirname "${BASH_SOURCE[0]}"): wrong under sbatch. Use $SLURM_SUBMIT_DIR.
set -euo pipefail
DIR="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$DIR"; mkdir -p logs
module load minimap2/2.30
module load ncbi-blast/2.14.0+
module load samtools
module load parallel

P=/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci
REF=$P/Genotyping/all_C_immitis_ref_RS/genome/FungiDB-68_CimmitisRS_Genome.fasta
AAFTF=$P/Assembly/asm/AAFTF
ANN=$P/Assembly/annotation
CPU=${SLURM_CPUS_ON_NODE:-2}

WORK="${SCRATCH:?SCRATCH is not set}/asmlocus.$$"; mkdir -p "$WORK/out"
samtools faidx "$REF" GG704914:968094-973690 | sed '1s/.*/>RS_CIMG_04613_pm2kb/' > "$WORK/query.fa"
cp "$WORK/query.fa" asm_locus_query.fa

run_one () {
  s=$1
  a="$AAFTF/$s.sorted.fasta"
  if [ -s "$a" ]; then
    minimap2 -x asm20 -c --cs -t 1 "$WORK/query.fa" "$a" > "$WORK/out/$s.asm.paf" 2>/dev/null
    blastn -query "$WORK/query.fa" -subject "$a" -evalue 1e-20 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
      > "$WORK/out/$s.asm.blastn.tsv" 2>/dev/null
  fi
  # contigs that the AAFTF mitochondrial screen removed from the nuclear assembly
  mt="$AAFTF/$s.vecscreen.mitochondria.fasta"
  if [ -s "$mt" ]; then
    minimap2 -x asm20 -c --cs -t 1 "$WORK/query.fa" "$mt" > "$WORK/out/$s.mito.paf" 2>/dev/null
    grep -v '>' "$mt" | tr -d '\n' | wc -c > "$WORK/out/$s.mito.len"
  fi
  n=$(ls "$ANN/$s"/annotate_results/*.scaffolds.fa 2>/dev/null | head -1 || true)
  if [ -n "$n" ] && [ -s "$n" ]; then
    minimap2 -x asm20 -c --cs -t 1 "$WORK/query.fa" "$n" > "$WORK/out/$s.ann.paf" 2>/dev/null
    blastn -query "$WORK/query.fa" -subject "$n" -evalue 1e-20 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
      > "$WORK/out/$s.ann.blastn.tsv" 2>/dev/null
    # N content of the annotated scaffold region that aligns to the query, for scaffold-gap checks
    # PAF columns 1,3,4 = scaffold name, start, end (the scaffold is the minimap2 query here)
    awk '{print $1"\t"$3"\t"$4}' "$WORK/out/$s.ann.paf" | sort -u > "$WORK/out/$s.ann.regions" || true
    # (copy to $SCRATCH first: faidx writes a .fai, and the shared tree must not be touched)
    if [ -s "$WORK/out/$s.ann.regions" ]; then
      cp "$n" "$WORK/$s.ann.fa"
      while read -r c st en; do
        printf "%s\t%s\t%s\t" "$c" "$st" "$en"
        samtools faidx "$WORK/$s.ann.fa" "$c:$((st+1))-$en" 2>/dev/null | tail -n +2 | tr -d '\n' | tr -cd 'Nn' | wc -c
      done < "$WORK/out/$s.ann.regions" > "$WORK/out/$s.ann.nrun.tsv"
      rm -f "$WORK/$s.ann.fa" "$WORK/$s.ann.fa.fai"
    fi
  fi
}
export -f run_one
export WORK AAFTF ANN

{ ls "$AAFTF"/*.sorted.fasta | sed 's#.*/##; s/\.sorted\.fasta$//'; ls "$ANN"; } | sort -u > "$WORK/strains.txt"
echo "strains: $(wc -l < "$WORK/strains.txt")"
parallel -j "$CPU" run_one :::: "$WORK/strains.txt"

tar -C "$WORK" -czf asm_locus_raw.tar.gz out
echo "done: asm_locus_raw.tar.gz"
rm -rf "$WORK"
