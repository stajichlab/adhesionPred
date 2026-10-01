#!/bin/bash -l
# J1 (spec 8, D2): SignalP 6 (GPU build, fast mode) and PredGPI on every unique sequence.
# Submit (see the plan's run section):
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=... -o <log> -e <log> j1_features.sh
# Input: $STEP1_WORKDIR/phaseb/unique_sequences.fasta.gz (05). The FASTA is split into
# J1_PARTS parts (default 8; record k of N goes to part floor(k*parts/N)). SignalP runs the
# parts one after another on the GPU; PredGPI (CPU only, single-threaded) runs the parts in
# parallel at the same time. Outputs go to $STEP1_WORKDIR/phaseb/signalp/part_NNN/ and
# $STEP1_WORKDIR/phaseb/predgpi/part_NNN.tsv.gz. A part whose output exists and whose
# recorded input SHA-256 equals the current part is skipped, so a rerun resumes.
# One job, not one job per part: measured rates (SignalP 6 fast on gpu12 about 192 proteins/s,
# job 29280458: 72,326 proteins in 6 min 16 s; PredGPI about 31 proteins/s per core on c01)
# give minutes per part, and the global job-size rule says not to split work into jobs of
# minutes. GPU build, not CPU build: the CPU build (signalp/6, weights
# distilled_model_signalp6.pt installed 2026-09-30) ran --mode fast at about 2.1 s per protein
# on 2 threads of c01 (2026-10-01), so 69,941 sequences would take about 41 h. The GPU build
# refuses to run without a CUDA device (verified 2026-10-01).
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH -J step1_j1
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_j1"
S1="$PROJ_ROOT/analysis/step1_compare"
OUT="$STEP1_WORKDIR/phaseb"
PARTS="${J1_PARTS:-8}"
SIGNALP_MODULE="${J1_SIGNALP_MODULE:-signalp/6-gpu}"  # tests pass a stub modulefile path
PREDGPI_MODULE="${J1_PREDGPI_MODULE:-predgpi/202001}"
mkdir -p "$TMP/in" "$TMP/sp" "$TMP/gpi" "$OUT/signalp" "$OUT/predgpi"

zcat "$OUT/unique_sequences.fasta.gz" > "$TMP/all.fasta"
N=$(grep -c '^>' "$TMP/all.fasta")
awk -v n="$N" -v p="$PARTS" -v dir="$TMP/in" '
  /^>/ { k++; part = int((k - 1) * p / n) }
  { printf "%s\n", $0 > sprintf("%s/part_%03d.fasta", dir, part) }' "$TMP/all.fasta"
echo "J1: $N sequences in $PARTS parts"

part_sha() { sha256sum "$TMP/in/$1.fasta" | cut -d' ' -f1; }

# copy_atomic SRC DEST: copy to a temp name next to DEST, then rename.
copy_atomic() { cp "$1" "$(dirname "$2")/.tmp.$(basename "$2")" && \
  mv "$(dirname "$2")/.tmp.$(basename "$2")" "$2"; }

run_predgpi() {  # $1 = part name; runs in a subshell with only the PredGPI module
  local part="$1" sha dest="$OUT/predgpi/$1.tsv.gz"
  sha=$(part_sha "$part")
  if [ -s "$dest" ] && [ "$(cat "$OUT/predgpi/$part.input.sha256" 2>/dev/null)" = "$sha" ]; then
    echo "  PredGPI skip $part (done)"; return 0
  fi
  python "$S1/jobs/predgpi_scores.py" --fasta "$TMP/in/$part.fasta" --out "$TMP/gpi/$part.tsv"
  gzip -n -c "$TMP/gpi/$part.tsv" > "$TMP/gpi/$part.tsv.gz"
  echo "$sha" > "$TMP/gpi/$part.input.sha256"
  copy_atomic "$TMP/gpi/$part.input.sha256" "$OUT/predgpi/$part.input.sha256"
  copy_atomic "$TMP/gpi/$part.tsv.gz" "$dest"
  echo "  PredGPI done $part"
}
export -f run_predgpi part_sha copy_atomic
export TMP OUT S1

set +e
(
  set -euo pipefail
  module load "$PREDGPI_MODULE"
  ls "$TMP/in" | sed 's/\.fasta$//' | xargs -P "$PARTS" -I{} bash -c 'set -euo pipefail; run_predgpi "$1"' _ {}
) > "$TMP/predgpi.log" 2>&1 &
GPI_PID=$!

(
  set -euo pipefail
  module load "$SIGNALP_MODULE"
  command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH" >&2; exit 1; }
  for f in "$TMP"/in/part_*.fasta; do
    part=$(basename "$f" .fasta)
    sha=$(part_sha "$part")
    dest="$OUT/signalp/$part"
    if [ -s "$dest/prediction_results.txt.gz" ] && \
       [ "$(cat "$dest/input.sha256" 2>/dev/null)" = "$sha" ]; then
      echo "  SignalP skip $part (done)"; continue
    fi
    rm -rf "$TMP/sp/$part"
    signalp6 --fastafile "$f" --organism eukarya --output_dir "$TMP/sp/$part" \
      --format none --mode fast --write_procs 8 --torch_num_threads 8 2>&1 | tail -2
    mkdir -p "$dest"
    for name in prediction_results.txt output.gff3 region_output.gff3; do
      gzip -n -c "$TMP/sp/$part/$name" > "$TMP/sp/$part/$name.gz"
    done
    echo "$sha" > "$TMP/sp/$part/input.sha256"
    copy_atomic "$TMP/sp/$part/output.gff3.gz" "$dest/output.gff3.gz"
    copy_atomic "$TMP/sp/$part/region_output.gff3.gz" "$dest/region_output.gff3.gz"
    copy_atomic "$TMP/sp/$part/input.sha256" "$dest/input.sha256"
    # prediction_results.txt.gz last: it is the done marker of the part
    copy_atomic "$TMP/sp/$part/prediction_results.txt.gz" "$dest/prediction_results.txt.gz"
    echo "  SignalP done $part"
  done
)
SP_STATUS=$?

wait "$GPI_PID"
GPI_STATUS=$?
set -e
cat "$TMP/predgpi.log"
echo "J1 done: SignalP status $SP_STATUS, PredGPI status $GPI_STATUS"
exit $(( SP_STATUS != 0 || GPI_STATUS != 0 ))
