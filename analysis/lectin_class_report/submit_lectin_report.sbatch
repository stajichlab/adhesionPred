#!/bin/bash -l
# Optional SLURM wrapper for lectin_report.py. Do not submit without checking the queue.
# Submit:
#   sbatch --export=ALL,REPORT_SCRIPT=/abs/path/lectin_report.py,OUT_DIR=/abs/out,FASTA_LIST=/abs/list.txt \
#     submit_lectin_report.sbatch
# FASTA_LIST: one FASTA path per line (plain, .gz or .zst). Optional: SIGNALP_FILE, PFAM_HMM.
# No BASH_SOURCE is used. All paths come from the environment.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J lectin_report
set -euo pipefail
: "${REPORT_SCRIPT:?export REPORT_SCRIPT (absolute path of lectin_report.py)}"
: "${OUT_DIR:?export OUT_DIR}"
: "${FASTA_LIST:?export FASTA_LIST (file with one FASTA path per line)}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/lectin_report.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP"
module load hmmer/3.4
ARGS=(--out-dir "$TMP/out" --tmpdir "$TMP" --cpu "${SLURM_CPUS_PER_TASK:-8}")
[ -n "${PFAM_HMM:-}" ] && ARGS+=(--pfam-hmm "$PFAM_HMM")
[ -n "${SIGNALP_FILE:-}" ] && ARGS+=(--signalp "$SIGNALP_FILE")
/usr/bin/python3.12 "$REPORT_SCRIPT" "${ARGS[@]}" --fasta-list "$FASTA_LIST"
# copy back with a temporary name, then rename
mkdir -p "$OUT_DIR"
for f in "$TMP"/out/*; do
  b=$(basename "$f")
  cp "$f" "$OUT_DIR/.$b.tmp"
  mv "$OUT_DIR/.$b.tmp" "$OUT_DIR/$b"
done
