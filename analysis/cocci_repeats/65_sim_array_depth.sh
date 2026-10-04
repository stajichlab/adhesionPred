#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 32 --mem 48gb --time 2:00:00
#SBATCH -J simarray -o logs/65_sim_array_depth.%A.log
# Read simulation with KNOWN SOWgp unit counts, pushed through the same steps as the real data:
#   bwa-mem2 2.2.1 mem (defaults) -> samtools 1.19.2 fixmate/sort -> Picard 2.26.11
#   MarkDuplicates -> mosdepth 0.3.12 -n --by <array BED> (all reads, and -Q 20).
# These are the tool versions in the CRAM @PG lines and in pipeline/11c_mosdepth_regions.sh.
#
# Grid
#   RS edits (u = 2..6, 2-3 edit designs per u, cycled over replicates): 6 reps
#   Real long-read alleles (Silveira 2022 u4, CiB10637 u5, Cpos1038 u5, VFC140 u2): 3 reps
#   x read length / median insert: 100/420, 150/360, 250/510, 300/470 (from cram_readlen.tsv)
#   x coverage 30x and 75x; substitution error 0.2%; insert sd 120.
# Output: sim_array_depth.regions.tsv.gz (run, allele, L, cov, rep, mapq, BED region, depth)
#         sim_alleles.tsv; summarise with 66_sim_array_summary.py
# Run from analysis/cocci_repeats:  sbatch 65_sim_array_depth.sh
set -euo pipefail
DIR="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$DIR"; mkdir -p logs
module unload java 2>/dev/null || true
module load picard/2.26.11
module load samtools/1.19.2
module load bwa-mem2/2.2.1
module load mosdepth/0.3.12
module load parallel
module load workspace/scratch   # 'module unload java' also unloads this, which unsets $SCRATCH
PY=/usr/bin/python3.12
R=/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS
BED=$R/coverage/mosdepth_gene/CIMG_04613.array_regions.bed
CPU=${SLURM_CPUS_ON_NODE:-2}

WORK="${SCRATCH:?SCRATCH is not set}/simarray.$$"; mkdir -p "$WORK/ref" "$WORK/run"
cp "$R"/genome/FungiDB-68_CimmitisRS_Genome.fasta* "$WORK/ref/"
cp "$BED" "$WORK/array.bed"
$PY 65_sim_array_reads.py build --out "$WORK/sim_alleles"
cp "$WORK/sim_alleles.tsv" sim_alleles.tsv
export WORK PY DIR

# task list: allele L ins cov rep seed
$PY - "$WORK/sim_alleles.tsv" > "$WORK/tasks.txt" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
LI = [(100, 420), (150, 360), (250, 510), (300, 470)]
seed = 0
by_u = {}
for r in rows:
    if r["source"] == "RS edit":
        by_u.setdefault(int(r["true_units"]), []).append(r["allele"])
for u, al in sorted(by_u.items()):
    for L, ins in LI:
        for cov in (30, 75):
            for rep in range(6):
                seed += 1
                print(al[rep % len(al)], L, ins, cov, rep, seed)
for r in rows:
    if r["source"] != "RS edit":
        for L, ins in LI:
            for cov in (30, 75):
                for rep in range(3):
                    seed += 1
                    print(r["allele"], L, ins, cov, rep, seed)
PY
echo "tasks: $(wc -l < "$WORK/tasks.txt")"

one () {
  allele=$1; L=$2; ins=$3; cov=$4; rep=$5; seed=$6
  id="${allele}.L${L}.c${cov}.r${rep}"; d="$WORK/run/$id"; mkdir -p "$d"
  $PY "$DIR/65_sim_array_reads.py" sim --fasta "$WORK/sim_alleles.fa" --allele "$allele" \
      --len "$L" --cov "$cov" --ins "$ins" --seed "$seed" --out "$d/r" 2>/dev/null
  bwa-mem2 mem -t 2 -R "@RG\tID:$id\tSM:$id\tLB:$id\tPL:illumina" "$WORK/ref/FungiDB-68_CimmitisRS_Genome.fasta" \
      "$d/r_R1.fq.gz" "$d/r_R2.fq.gz" 2>/dev/null > "$d/a.sam"
  samtools fixmate -u -O BAM "$d/a.sam" "$d/fm.bam"
  samtools sort -O BAM -o "$d/srt.bam" -T "$d/tmp" "$d/fm.bam" 2>/dev/null
  java -Xmx2g -jar "$PICARD" MarkDuplicates -I "$d/srt.bam" -O "$d/dd.bam" -METRICS_FILE "$d/dd.metrics" \
      -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT 2>/dev/null
  mosdepth -n --by "$WORK/array.bed" -t 1 "$d/all" "$d/dd.bam"
  mosdepth -n --by "$WORK/array.bed" -Q 20 -t 1 "$d/q20" "$d/dd.bam"
  for m in all q20; do
    zcat "$d/$m.regions.bed.gz" | awk -v id="$id" -v a="$allele" -v L="$L" -v c="$cov" -v r="$rep" -v m="$m" \
      'BEGIN{OFS="\t"}{print id,a,L,c,r,m,$1,$2,$3,$4,$5}'
  done > "$WORK/run/$id.tsv"
  rm -rf "$d"
}
export -f one
parallel -j $((CPU/2)) --colsep ' ' one {1} {2} {3} {4} {5} {6} :::: "$WORK/tasks.txt"

{ printf "run\tallele\tread_len\tcov\trep\tmapq\tchrom\tstart\tend\tregion\tdepth\n"
  cat "$WORK"/run/*.tsv; } | gzip -c > sim_array_depth.regions.tsv.gz
echo "runs written: $(ls "$WORK"/run/*.tsv | wc -l)"
rm -rf "$WORK"
