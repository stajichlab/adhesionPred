#!/bin/bash
# Write the two manifests for the 7 Coccidioides proteomes.
#   proteomes.tsv  3 columns: name, fasta, signalp   (input of 01_known_family_hmm.sh)
#   manifest.tsv   4 columns: name, fasta, signalp, domtbl   (input of cys_candidates.py)
# Environment (no path is derived from the script location):
#   OUT_DIR  required
#   PROJ_ROOT  default /bigdata/stajichlab/jstajich/projects/adhesionPred (SignalP results)
#   LR   long-read proteome directory (default For_Marc, see config/site.yaml cocci_longread)
#   PAN  pangenome input_run2 directory (see config/site.yaml cocci_pangenome)
set -euo pipefail
stop() { echo "STOP: $*" >&2; exit 2; }
[ -n "${OUT_DIR:-}" ] || stop "OUT_DIR is not set"
PROJ_ROOT="${PROJ_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
LR="${LR:-/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc}"
PAN="${PAN:-/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input_run2}"
SP="$PROJ_ROOT/analysis/cocci_repeats/signalp"
mkdir -p "$OUT_DIR" || stop "cannot create $OUT_DIR"
{
  printf 'CimmitisRS_FungiDB\t%s\n' "$PAN/CimmitisRS_FungiDB.fasta"
  printf 'Coccidioides_immitis_CiB10637\t%s\n' "$LR/CiB10637/Coccidioides_immitis_CiB10637.proteins.fa"
  printf 'Coccidioides_immitis_CiB10992\t%s\n' "$LR/CiB10992/Coccidioides_immitis_CiB10992.proteins.fa"
  printf 'Coccidioides_immitis_VFC140\t%s\n' "$LR/VFC140/Coccidioides_immitis_VFC140.proteins.fa"
  printf 'Coccidioides_posadasii_Cpos1038\t%s\n' "$LR/Cpos1038/Coccidioides_posadasii_Cpos1038.proteins.fa"
  printf 'Coccidioides_posadasii_Cpos3700\t%s\n' "$LR/Cpos3700/Coccidioides_posadasii_Cpos3700.proteins.fa"
  printf 'CposadasiiSilveira2022_FungiDB\t%s\n' "$PAN/CposadasiiSilveira2022_FungiDB.fasta"
} > "$OUT_DIR/names_fasta.tmp"
: > "$OUT_DIR/proteomes.tsv.tmp"
: > "$OUT_DIR/manifest.tsv.tmp"
while IFS=$'\t' read -r name fasta; do
  [ -s "$fasta" ] || { rm -f "$OUT_DIR"/*.tmp; stop "FASTA not found: $fasta"; }
  sp="$SP/$name/prediction_results.txt"
  [ -s "$sp" ] || { rm -f "$OUT_DIR"/*.tmp; stop "SignalP file not found: $sp"; }
  printf '%s\t%s\t%s\n' "$name" "$fasta" "$sp" >> "$OUT_DIR/proteomes.tsv.tmp"
  printf '%s\t%s\t%s\t%s\n' "$name" "$fasta" "$sp" "$OUT_DIR/$name.domtbl.gz" >> "$OUT_DIR/manifest.tsv.tmp"
done < "$OUT_DIR/names_fasta.tmp"
rm -f "$OUT_DIR/names_fasta.tmp"
mv "$OUT_DIR/proteomes.tsv.tmp" "$OUT_DIR/proteomes.tsv"
mv "$OUT_DIR/manifest.tsv.tmp" "$OUT_DIR/manifest.tsv"
echo "wrote $OUT_DIR/proteomes.tsv and $OUT_DIR/manifest.tsv"
