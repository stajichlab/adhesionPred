#!/bin/bash -l
# SignalP 6 over the long-read Coccidioides proteomes + the pangenome references.
# Long-read assemblies matter here: tandem repeat arrays are collapsed or broken in
# short-read assemblies, and repeat copy number is the feature we are trying to measure.
#
# Runs on short_gpu with the GPU build of SignalP 6 (short has a 2 h wall limit).
#   sbatch analysis/cocci_repeats/01_signalp.sh
#
# Output location: $PROJ_ROOT/analysis/cocci_repeats/signalp/<proteome>/ -- an absolute
# path, NOT a path relative to the submit directory. Before 2026-09-30 this script wrote
# to ./signalp relative to $SLURM_SUBMIT_DIR, so the output landed wherever the job was
# submitted from and never reached the repository. 03_repeat_surface_candidates.py then
# had no SignalP directory to read, which is why the secreted subset could not be
# recomputed for the 14_repeat_detect_general.py calls.
# No --constraint=gpu_latest. SignalP 6 in fast mode is a small distilled model and runs on
# any CUDA GPU in this partition; requiring the newest cards pushed the scheduled start out by
# ~5 h on 2026-09-30 for no gain. Re-add the constraint only if a specific card is needed.
# Partition: exfab, not short_gpu. exfab serves GPUs when asked with --gres=gpu:1 and is
# usually far less contended -- on 2026-09-30 short_gpu scheduled ~4 h out with 4 nodes down,
# while exfab had 0 jobs pending. A comma-separated partition list is NOT usable here: this
# account has a partition set in its SLURM association, so multi-partition requests are
# rejected with "Multiple partition job request not supported when a partition is set in the
# association". Switch this single value to short_gpu if exfab is busy instead.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:55:00
#SBATCH -J cocci_signalp
#SBATCH -o logs/signalp.%j.log
#SBATCH -e logs/signalp.%j.log
set -uo pipefail

PROJ_ROOT="${PROJ_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
source "$PROJ_ROOT/analysis/_common/paths.sh"

HERE=$PROJ_ROOT/analysis/cocci_repeats
SPDIR=$HERE/signalp
mkdir -p "$PROJ_ROOT/logs" "$SPDIR"
cd "$HERE"

# Load the GPU build and check it. The CPU build (signalp/6, signalp/6.0i) ships with an
# EMPTY model_weights/ directory on this cluster -- verified 2026-09-30 -- so it cannot run
# any mode, and --mode fast fails with:
#   FileNotFoundError: Fast mode requires model to be installed at .../distilled_model_signalp6.pt
# The previous `module load signalp/6-gpu || module load signalp/6` fallback therefore led to
# a build that cannot work. Fail loudly instead of falling back.
module load signalp/6-gpu 2>/dev/null || module load signalp/6.0i-gpu
command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH after module load" >&2; exit 1; }
SP_WEIGHTS=$(dirname "$(dirname "$(command -v signalp6)")")/lib/python*/site-packages/signalp/model_weights/distilled_model_signalp6.pt
# shellcheck disable=SC2086
if ! ls $SP_WEIGHTS >/dev/null 2>&1; then
  echo "FATAL: signalp6 at $(command -v signalp6) has no distilled model weights." >&2
  echo "       --mode fast cannot run. This partition needs the GPU build (signalp/6-gpu)." >&2
  exit 1
fi

LR=$COCCI_LONGREAD
PAN=$COCCI_PANGENOME/input_run2

for FA in "$LR"/*/*.proteins.fa "$PAN"/CimmitisRS_FungiDB.fasta "$PAN"/CposadasiiSilveira2022_FungiDB.fasta; do
  [ -s "$FA" ] || continue
  NAME=$(basename "$FA" | sed -E 's/\.(proteins\.fa|fasta)$//')
  OUT=$SPDIR/$NAME
  if [ -s "$OUT/prediction_results.txt" ]; then echo "  skip $NAME (done)"; continue; fi
  mkdir -p "$OUT"
  echo "  SignalP: $NAME ($(grep -c '>' "$FA") proteins)"
  signalp6 --fastafile "$FA" --organism eukarya --output_dir "$OUT" \
           --format none --mode fast --write_procs 8 2>&1 | tail -2
done

# Per-proteome counts, and a tracked summary so the run is recorded even though the
# per-protein tables are gitignored.
# prediction_results.txt is TAB separated with a "# ID<TAB>Prediction<TAB>..." header. Parse
# with FS='\t' explicitly -- a default-whitespace split mis-columns the "CS pos: 20-21. Pr:"
# field and can make every row look like a signal peptide. The per-class distribution is
# printed as well, so a mis-parse is visible instead of silently producing a wrong count.
SUMMARY=$HERE/signalp_summary.tsv
printf 'proteome\tn_proteins\tn_signal_peptide\tfrac\n' > "$SUMMARY"
for d in "$SPDIR"/*/; do
  [ -s "$d/prediction_results.txt" ] || continue
  name=$(basename "$d")
  read -r tot n < <(awk -F'\t' '
      $0 !~ /^#/ && NF >= 2 { t++; if ($2 != "OTHER") s++ }
      END { printf "%d %d\n", t, s }' "$d/prediction_results.txt")
  frac=$(awk -v a="$n" -v b="$tot" 'BEGIN{ printf (b>0 ? "%.3f" : "NA"), (b>0 ? a/b : 0) }')
  printf '%s\t%s\t%s\t%s\n' "$name" "$tot" "$n" "$frac" >> "$SUMMARY"
  echo "  $name: $n / $tot with a signal peptide (frac $frac)"
  echo "    class distribution:"
  awk -F'\t' '$0 !~ /^#/ && NF >= 2 { c[$2]++ } END { for (k in c) printf "      %-16s %d\n", k, c[k] }' \
      "$d/prediction_results.txt"
done
echo
echo "Sanity check: for a fungal proteome expect roughly 5-12% with a signal peptide."
echo "A fraction near 0.00 or near 1.00 means the parse or the run is wrong -- do not"
echo "feed such a directory to 03_repeat_surface_candidates.py."

# prediction_results.txt is one line per protein and reaches tens of MB across 7
# proteomes. Keep a compressed copy next to each (see CLAUDE.md storage guidance).
for d in "$SPDIR"/*/; do
  [ -s "$d/prediction_results.txt" ] || continue
  gzip -kf "$d/prediction_results.txt"
done
echo "summary: $SUMMARY"
