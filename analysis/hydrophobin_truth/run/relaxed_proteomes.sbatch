#!/bin/bash -l
# L7a: hmmsearch --nobias --cut_ga of the frozen hydrophobin_relaxed.hmm on every proteome of RUN_LIST (one job), writing
# <workdir>/raw/hydrophobin_relaxed/domtbl.txt. The module wrapper reads it afterwards.
# Submit: PROJ_ROOT=... sbatch --export=ALL relaxed_proteomes.sbatch
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J hyd_relaxed_prot
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
module load hmmer/3.4
HMM="$PROJ_ROOT/data/sorting_hat/hydrophobin_relaxed.hmm"
tail -n +2 "$PROJ_ROOT/analysis/hydrophobin_truth/run_list.tsv" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  out="$workdir/raw/hydrophobin_relaxed"
  mkdir -p "$out"
  hmmsearch --nobias --cut_ga --cpu 8 --noali -o /dev/null --domtblout "$out/.tmp.domtbl.txt" "$HMM" "$fasta"
  tail -n 3 "$out/.tmp.domtbl.txt" | grep -qx '# \[ok\]' || { echo "FATAL: no ok trailer for $name" >&2; exit 1; }
  mv "$out/.tmp.domtbl.txt" "$out/domtbl.txt"
  echo "done $name"
done
