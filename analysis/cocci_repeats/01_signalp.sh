#!/bin/bash -l
# SignalP 6 over the long-read Coccidioides proteomes + the pangenome references.
# Long-read assemblies matter here: tandem repeat arrays are collapsed or broken in
# short-read assemblies, and repeat copy number is the feature we are trying to measure.
#
# Runs on short_gpu with the GPU build of SignalP 6 (short has a 2 h wall limit).
#   sbatch 01_signalp.sh
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu_latest
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:55:00
#SBATCH -J cocci_signalp
#SBATCH -o logs/signalp.%j.log
#SBATCH -e logs/signalp.%j.log
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs signalp
module load signalp/6-gpu 2>/dev/null || module load signalp/6

LR=/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc
PAN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input_run2

for FA in "$LR"/*/*.proteins.fa "$PAN"/CimmitisRS_FungiDB.fasta "$PAN"/CposadasiiSilveira2022_FungiDB.fasta; do
  [ -s "$FA" ] || continue
  NAME=$(basename "$FA" | sed -E 's/\.(proteins\.fa|fasta)$//')
  OUT=signalp/$NAME
  if [ -s "$OUT/prediction_results.txt" ]; then echo "  skip $NAME (done)"; continue; fi
  mkdir -p "$OUT"
  echo "  SignalP: $NAME ($(grep -c '>' "$FA") proteins)"
  signalp6 --fastafile "$FA" --organism eukarya --output_dir "$OUT" \
           --format none --mode fast --write_procs 8 2>&1 | tail -2
done
echo; for d in signalp/*/; do
  n=$(awk '$2!="OTHER" && $1!~/^#/' "$d/prediction_results.txt" 2>/dev/null | wc -l)
  echo "  $(basename "$d"): $n proteins with a signal peptide"
done
