#!/bin/bash -l
# Align the six C. immitis RS RNA-seq samples (UCSD mycelia / spherule) to the RefSeq genome with STAR.
# No annotation is given to STAR, so junctions are found from the reads alone and the annotated gene
# models can be compared with them (20_ ... 21_exon_support.py). Two-pass mode finds new junctions.
# Intron length limits are set for a fungus. Outputs: <OUT>/<sample>/{Aligned.sortedByCoord.out.bam,SJ.out.tab}.
# A sample whose BAM exists is skipped, so a rerun resumes.
# Usage: sbatch --export=ALL,OUT=<dir> 20_star_align.sh
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH -J cocci_star
set -euo pipefail
: "${OUT:?export OUT (output dir)}"
TMP="${SCRATCH:?SCRATCH is not set}/cocci_star"
REF=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/cocci_spherule/ref/GCF_000149335.2_ASM14933v2_genomic.fna.gz
FQ=/bigdata/stajichlab/shared/projects/Coccidioides/csSeq/UCSD_SpheruleMycelium/fastq
SAMPLES=/bigdata/stajichlab/shared/projects/Coccidioides/csSeq/UCSD_SpheruleMycelium/samples.csv
mkdir -p "$TMP" "$OUT"
module load star/2.7.11b samtools/1.19.2

zcat "$REF" > "$TMP/genome.fa"
STAR --runMode genomeGenerate --genomeDir "$TMP/index" --genomeFastaFiles "$TMP/genome.fa" \
  --genomeSAindexNbases 11 --runThreadN 16 --outFileNamePrefix "$TMP/idx_" > "$TMP/idx.log"

tail -n +2 "$SAMPLES" | while IFS=, read -r SAMPLE COND REP BASE; do
  dest="$OUT/${COND}.r${REP}"
  if [ -s "$dest/Aligned.sortedByCoord.out.bam" ]; then echo "skip $dest"; continue; fi
  mkdir -p "$TMP/$SAMPLE"
  STAR --genomeDir "$TMP/index" --runThreadN 16 \
    --readFilesIn "$FQ/${BASE}_R1_001.fastq.gz" "$FQ/${BASE}_R2_001.fastq.gz" --readFilesCommand zcat \
    --twopassMode Basic --alignIntronMin 20 --alignIntronMax 5000 --alignMatesGapMax 5000 \
    --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM \
    --outFileNamePrefix "$TMP/$SAMPLE/" > "$TMP/$SAMPLE/star.log"
  samtools index -@ 8 "$TMP/$SAMPLE/Aligned.sortedByCoord.out.bam"
  mkdir -p "$dest"
  cp "$TMP/$SAMPLE/SJ.out.tab" "$TMP/$SAMPLE/Log.final.out" "$dest/"
  cp "$TMP/$SAMPLE/Aligned.sortedByCoord.out.bam.bai" "$dest/"
  cp "$TMP/$SAMPLE/Aligned.sortedByCoord.out.bam" "$dest/"
  rm -rf "$TMP/$SAMPLE"
  echo "done $dest"
done
echo finished
