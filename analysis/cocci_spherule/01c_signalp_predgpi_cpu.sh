#!/bin/bash -l
# SignalP 6 (CPU build, fast mode) and PredGPI on the C. immitis RS RefSeq proteins.
# CPU variant of analysis/step1_compare/jobs/j1_features.sh, for a small input (9,910 proteins)
# when no GPU is free. Same output layout, so parse_signalp and parse_predgpi_scores read it.
# Measured earlier (step1 j1 notes, 2026-10-01): CPU SignalP --mode fast about 2.1 s per protein
# on 2 threads. 8 parts of about 1,240 proteins, 8 parts in parallel, 2 threads each: about 45 min.
# Usage: sbatch --export=ALL,PROJ_ROOT=<repo>,STEP1_WORKDIR=<workdir> 01c_signalp_predgpi_cpu.sh
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH -J cocci_sp_cpu
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR}"
TMP="${SCRATCH:?SCRATCH is not set; run as a SLURM job}/cocci_sp_cpu"
S1="$PROJ_ROOT/analysis/step1_compare"
OUT="$STEP1_WORKDIR/phaseb"
PARTS=8
mkdir -p "$TMP/in" "$TMP/sp" "$TMP/gpi" "$OUT/signalp" "$OUT/predgpi"
zcat "$OUT/unique_sequences.fasta.gz" > "$TMP/all.fasta"
N=$(grep -c '^>' "$TMP/all.fasta")
awk -v n="$N" -v p="$PARTS" -v dir="$TMP/in" '
  /^>/ { k++; part = int((k - 1) * p / n) }
  { printf "%s\n", $0 > sprintf("%s/part_%03d.fasta", dir, part) }' "$TMP/all.fasta"
echo "$N sequences in $PARTS parts"

run_sp() {  # $1 = part name
  local part="$1" dest="$OUT/signalp/$1"
  if [ -s "$dest/prediction_results.txt.gz" ]; then echo "  SignalP skip $part"; return 0; fi
  signalp6 --fastafile "$TMP/in/$part.fasta" --organism eukarya --output_dir "$TMP/sp/$part" \
    --format none --mode fast --write_procs 2 --torch_num_threads 2 > "$TMP/sp/$part.log" 2>&1
  mkdir -p "$dest"
  for name in prediction_results.txt output.gff3 region_output.gff3; do
    gzip -n -c "$TMP/sp/$part/$name" > "$dest/$name.gz"
  done
  echo "  SignalP done $part"
}
export -f run_sp
export TMP OUT

(
  module load signalp/6
  command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH" >&2; exit 1; }
  ls "$TMP/in" | sed 's/\.fasta$//' | xargs -P "$PARTS" -I{} bash -c 'set -euo pipefail; run_sp "$1"' _ {}
)

(
  module load predgpi/202001
  for f in "$TMP"/in/part_*.fasta; do
    part=$(basename "$f" .fasta)
    python "$S1/jobs/predgpi_scores.py" --fasta "$f" --out "$TMP/gpi/$part.tsv"
    gzip -n -c "$TMP/gpi/$part.tsv" > "$OUT/predgpi/$part.tsv.gz"
    echo "  PredGPI done $part"
  done
)
echo "done"
