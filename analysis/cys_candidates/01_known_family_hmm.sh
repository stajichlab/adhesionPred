#!/bin/bash -l
# Cut five Pfam models and search them against the Coccidioides proteomes (--cut_ga).
# Models: PF05730 (CFEM), PF04681 (Bys1), PF01185 and PF06766 (hydrophobins), PF28404 (PRA3 family).
# Small job (a few minutes on one CPU). Not submitted by this script.
# Environment (no path is derived from the script location):
#   PROJ_ROOT  primary checkout, default /bigdata/stajichlab/jstajich/projects/adhesionPred
#   OUT_DIR    output directory (domtbl.gz per proteome and known_families.hmm)
#   MANIFEST   TSV with columns name, fasta (extra columns are ignored); default
#              $OUT_DIR/proteomes.tsv
#   PFAM_HMM   default /bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm
#   CPU        default 4
# Optional SLURM header, only if you run this as a job:
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --time=1:00:00
set -euo pipefail

PROJ_ROOT="${PROJ_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
OUT_DIR="${OUT_DIR:?STOP: OUT_DIR is not set}"
PFAM_HMM="${PFAM_HMM:-/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm}"
MANIFEST="${MANIFEST:-$OUT_DIR/proteomes.tsv}"
CPU="${CPU:-4}"

[ -s "$PFAM_HMM" ] || { echo "STOP: Pfam-A.hmm not found: $PFAM_HMM" >&2; exit 2; }
[ -s "$MANIFEST" ] || { echo "STOP: manifest not found: $MANIFEST" >&2; exit 2; }
module load hmmer/3.4
mkdir -p "$OUT_DIR"

MODELS="$OUT_DIR/known_families.hmm"
MISSING="$OUT_DIR/missing_models.txt"
: > "$MODELS.tmp"
: > "$MISSING.tmp"
for acc in PF05730 PF04681 PF01185 PF06766 PF28404; do
  # hmmfetch needs the exact accession with version; find it in the file header lines.
  full=$(grep -m1 -E "^ACC +${acc}\." "$PFAM_HMM" | awk '{print $2}' || true)
  if [ -z "$full" ]; then
    # PF28404 is not in Pfam 38.0 (2025-07), the release installed here. Do not stop: the
    # other four models are still useful. Tell cys_candidates.py with --missing-models.
    echo "WARN: $acc not in $PFAM_HMM; skipped (use --missing-models for it)" >&2
    echo "$acc" >> "$MISSING.tmp"
    continue
  fi
  hmmfetch "$PFAM_HMM" "$full" >> "$MODELS.tmp" || { rm -f "$MODELS.tmp"; exit 2; }
done
mv "$MODELS.tmp" "$MODELS"
mv "$MISSING.tmp" "$MISSING"
echo "models: $(grep -c '^NAME' "$MODELS")"

while IFS=$'\t' read -r name fasta _rest; do
  case "$name" in ''|'#'*) continue ;; esac
  [ -s "$fasta" ] || { echo "STOP: FASTA not found for $name: $fasta" >&2; exit 2; }
  tmp="$OUT_DIR/$name.domtbl.tmp"
  hmmsearch --cut_ga --cpu "$CPU" --noali -o /dev/null --domtblout "$tmp" "$MODELS" "$fasta"
  gzip -n -c "$tmp" > "$OUT_DIR/$name.domtbl.gz.tmp"
  mv "$OUT_DIR/$name.domtbl.gz.tmp" "$OUT_DIR/$name.domtbl.gz"
  rm -f "$tmp"
  echo "$name: $(zcat "$OUT_DIR/$name.domtbl.gz" | grep -vc '^#') domain hits"
done < "$MANIFEST"
