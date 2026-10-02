#!/bin/bash -l
# hmmsearch of the PF28404 (ARB_05178) model against every proteome of proteomes.tsv.
# Inputs by environment (no script-relative paths; this runs under SLURM):
#   WORK  directory that holds proteomes.tsv (from 00_select_proteomes.py); outputs go to $WORK
#   HMM   hmm file with the model (default: $WORK/pf28404.hmm, made here from KNOWN_HMM)
#   KNOWN_HMM  cys_candidates known_families.hmm (used only when $HMM does not exist)
# Submit: sbatch -p highclock -c 16 --mem=16G -t 2:00:00 -J pf28404 --export=ALL,WORK=... 01_hmmsearch.sh
# One job: ~830 proteomes, each a few seconds. Temp files go to node-local $SCRATCH.
#SBATCH -p highclock
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -J pf28404
set -euo pipefail
: "${WORK:?export WORK}"
TMP="${SCRATCH:?SCRATCH is not set; run as a SLURM job}/pf28404"
mkdir -p "$TMP" "$WORK/domtbl"
module load hmmer/3.4
HMM="${HMM:-$WORK/pf28404.hmm}"
if [ ! -s "$HMM" ]; then
  : "${KNOWN_HMM:?export KNOWN_HMM or HMM}"
  hmmfetch "$KNOWN_HMM" ARB_05178 > "$HMM"
fi
grep -E '^(NAME|ACC|GA) ' "$HMM" > "$WORK/pf28404_model_info.txt"
sha256sum "$HMM" >> "$WORK/pf28404_model_info.txt"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export HMM WORK TMP
run_one() {
  label=$1; fa=$2
  out="$WORK/domtbl/$label.domtbl"
  [ -s "$out" ] && return 0
  hmmsearch --cut_ga --noali --cpu 1 -o /dev/null --domtblout "$TMP/$label.domtbl" "$HMM" "$fa"
  mv "$TMP/$label.domtbl" "$out"
}
export -f run_one
awk -F'\t' 'NR>1 {printf "%s%c%s%c", $1, 0, $2, 0}' "$WORK/proteomes.tsv" | xargs -0 -n2 -P "$CPUS" bash -c 'run_one "$0" "$1"'
echo "done: $(ls "$WORK/domtbl" | wc -l) domtbl files"
