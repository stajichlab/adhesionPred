#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 8 --mem 24gb --time 1:30:00
#SBATCH -J repgen -o logs/17_repeat_detect_real.%A.log
# Run the family-agnostic detector (14) over the same 71,044 Coccidioides proteins that
# 02_repeat_profile.py was run on: the 5 UArizona long-read proteomes plus the C. immitis RS
# and C. posadasii Silveira references.
#
# Submit from analysis/cocci_repeats:   sbatch 17_repeat_detect_real.sh
#
# NOTE: do not use $(dirname "${BASH_SOURCE[0]}") here. Under sbatch the script is a spooled
# copy and that path is wrong. $SLURM_SUBMIT_DIR is the submit directory.

set -euo pipefail

DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$DIR"
mkdir -p logs

PAN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input_run2
LR=/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc
WORK="${SCRATCH:?SCRATCH is not set}/repgen.$$"
mkdir -p "$WORK"

# one process per proteome; the long-read set and the references are reported separately,
# as in 02, because tandem arrays collapse in short-read assemblies.
for FA in "$LR"/*/*.proteins.fa; do
  B=$(basename "$FA" .proteins.fa)
  /usr/bin/python3.12 14_repeat_detect_general.py "$FA" --out "$WORK/$B.tsv" \
    >"$WORK/$B.err" 2>&1 &
done
for FA in "$PAN"/CimmitisRS_FungiDB.fasta "$PAN"/CposadasiiSilveira2022_FungiDB.fasta; do
  B=$(basename "$FA" .fasta)
  /usr/bin/python3.12 14_repeat_detect_general.py "$FA" --out "$WORK/$B.tsv" \
    >"$WORK/$B.err" 2>&1 &
done
wait
cat "$WORK"/*.err

# merge: long-read set and reference set, matching the two files 02 produced
head -1 "$WORK/Coccidioides_immitis_CiB10637.tsv" > repeat_general_longread.tsv
for B in Coccidioides_immitis_CiB10637 Coccidioides_immitis_CiB10992 Coccidioides_immitis_VFC140 \
         Coccidioides_posadasii_Cpos1038 Coccidioides_posadasii_Cpos3700; do
  tail -n +2 "$WORK/$B.tsv" >> repeat_general_longread.tsv
done
head -1 "$WORK/CimmitisRS_FungiDB.tsv" > repeat_general_reference.tsv
tail -n +2 "$WORK/CimmitisRS_FungiDB.tsv" >> repeat_general_reference.tsv
tail -n +2 "$WORK/CposadasiiSilveira2022_FungiDB.tsv" >> repeat_general_reference.tsv

gzip -kf repeat_general_longread.tsv repeat_general_reference.tsv
wc -l repeat_general_longread.tsv repeat_general_reference.tsv
rm -rf "$WORK"
