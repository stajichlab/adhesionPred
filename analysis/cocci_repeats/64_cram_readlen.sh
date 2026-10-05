#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 32 --mem 16gb --time 1:00:00
#SBATCH -J readlen -o logs/64_cram_readlen.%A.log
# Read length and insert size per strain, from the RS-reference CRAMs.
# Sample: primary, mapped reads in GG704914:900,000-1,040,000 (the 140 kb around CIMG_04613).
# Output: cram_readlen.tsv  strain, n_reads, modal_len, median_len, max_len,
#         median_insert (proper pairs, TLEN > 0), iqr_insert
# Used as covariates (read length) in 61, and to set read lengths for the simulation (65).
# Run from analysis/cocci_repeats:  sbatch 64_cram_readlen.sh
set -euo pipefail
DIR="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$DIR"; mkdir -p logs
module load samtools/1.19.2
module load parallel
R=/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS
export REF=$R/genome/FungiDB-68_CimmitisRS_Genome.fasta
CPU=${SLURM_CPUS_ON_NODE:-2}

one () {
  c=$1; s=$(basename "$c" .cram)
  samtools view -T "$REF" -F 0x904 "$c" GG704914:900000-1040000 | \
  awk -v s="$s" 'BEGIN{OFS="\t"}
    { L=length($10); n++; len[n]=L; cnt[L]++; if (L>mx) mx=L
      if (and($2,2) && $9>0 && $9<5000) { m++; ins[m]=$9 } }
    END{
      best=0; for (k in cnt) if (cnt[k]>best) {best=cnt[k]; mode=k}
      asort(len); asort(ins)
      med = n ? len[int((n+1)/2)] : "NA"
      mi  = m ? ins[int((m+1)/2)] : "NA"
      iqr = m ? ins[int(m*0.75)+1]-ins[int(m*0.25)+1] : "NA"
      print s, n, mode, med, mx, mi, iqr }'
}
export -f one
{ printf "strain\tn_reads\tmodal_len\tmedian_len\tmax_len\tmedian_insert\tiqr_insert\n"
  parallel -j "$CPU" one ::: "$R"/aln/*.cram | sort; } > cram_readlen.tsv
echo "wrote cram_readlen.tsv: $(($(wc -l < cram_readlen.tsv)-1)) strains"
