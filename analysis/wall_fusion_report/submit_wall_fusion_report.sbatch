#!/bin/bash -l
# Optional SLURM wrapper for wall_fusion_report.py. Do not submit without checking the queue.
# Submit:
#   sbatch --export=ALL,REPORT_SCRIPT=/abs/path/wall_fusion_report.py,MODELS_TSV=/abs/wall_fusion_models.tsv,OUT_DIR=/abs/out,FASTA_LIST=/abs/list.txt \
#     submit_wall_fusion_report.sbatch
# FASTA_LIST: one FASTA path per line (plain, .gz or .zst). Optional: GENE_SET (wall_remodeling,
# set names, comma separated, or all; default all), PFAM_HMM.
# DRAFT extension (custom Ecm33 model, report only; nothing has been submitted or run): set
# CUSTOM_MODELS_TSV=/abs/wall_fusion_custom_models.tsv to also search the custom model, with
# CUSTOM_MODEL_HMM=/abs/data/custom_models/ecm33_all_2026-10-23/model.hmm. The job refuses to run
# if the sha256 of CUSTOM_MODEL_HMM is not the pinned one below (the script checks the table's
# sha256 again).
# No BASH_SOURCE is used. All paths come from the environment.
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -J wall_fusion_report
set -euo pipefail
: "${REPORT_SCRIPT:?export REPORT_SCRIPT (absolute path of wall_fusion_report.py)}"
: "${MODELS_TSV:?export MODELS_TSV (absolute path of wall_fusion_models.tsv)}"
: "${OUT_DIR:?export OUT_DIR}"
: "${FASTA_LIST:?export FASTA_LIST (file with one FASTA path per line)}"
ECM33_PINNED_SHA256=9a6e2c812540b0b1892f081611bf4044ec674bc524d68f5c6fc3d16001c50937
if [ -n "${CUSTOM_MODELS_TSV:-}" ]; then
  : "${CUSTOM_MODEL_HMM:?export CUSTOM_MODEL_HMM (absolute path of the custom model.hmm)}"
  got=$(sha256sum "$CUSTOM_MODEL_HMM" | cut -d' ' -f1)
  if [ "$got" != "$ECM33_PINNED_SHA256" ]; then
    echo "ERROR: sha256 of $CUSTOM_MODEL_HMM is $got, pinned $ECM33_PINNED_SHA256. Refusing." >&2
    exit 3
  fi
fi
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/wall_fusion_report.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP"
module load hmmer/3.4
ARGS=(--out-dir "$TMP/out" --tmpdir "$TMP" --cpu "${SLURM_CPUS_PER_TASK:-4}" --models-tsv "$MODELS_TSV" --set "${GENE_SET:-all}")
[ -n "${PFAM_HMM:-}" ] && ARGS+=(--pfam-hmm "$PFAM_HMM")
[ -n "${CUSTOM_MODELS_TSV:-}" ] && ARGS+=(--custom-models-tsv "$CUSTOM_MODELS_TSV")
/usr/bin/python3.12 "$REPORT_SCRIPT" "${ARGS[@]}" --fasta-list "$FASTA_LIST"
# copy back with a temporary name, then rename
mkdir -p "$OUT_DIR"
for f in "$TMP"/out/*; do
  b=$(basename "$f")
  cp "$f" "$OUT_DIR/.$b.tmp"
  mv "$OUT_DIR/.$b.tmp" "$OUT_DIR/$b"
done
