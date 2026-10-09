#!/bin/bash
# Turn the raw tool outputs of every proteome in RUN_LIST into module tables and run the core command (needs the csh venv on PATH).
# Usage: PROJ_ROOT=... RUN_LIST=... bash scan_converters.sh [START END]
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${RUN_LIST:?export RUN_LIST}"
START="${1:-1}"; END="${2:-100000}"
S="$PROJ_ROOT/_workdir/sorting_hat"
T="$S/taxdump"
M=cellsurface_sorting_hat_module
tail -n +2 "$RUN_LIST" | sed -n "${START},${END}p" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  echo "== $name"
  PFAMJ="$workdir/raw/pfam/provenance.json"
  SHA=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['sha256'])")
  HMMER=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['hmmer'].split()[2])")
  $M signalp --fasta "$fasta" --workdir "$workdir" --results "$workdir/raw/signalp/prediction_results.txt" --signalp-version "$(paste -sd';' "$workdir/raw/signalp/version.txt")" > /dev/null
  $M tm --fasta "$fasta" --workdir "$workdir" --table "$workdir/raw/tmhmm/tmhmm.tsv" > /dev/null
  $M pfam --fasta "$fasta" --workdir "$workdir" --domtbl "$workdir/raw/pfam/domtbl.txt" --family-table "$PROJ_ROOT/data/sorting_hat/family_table.tsv" \
    --pfam-release 38.2 --pfam-sha256 "$SHA" --hmmer-version "$HMMER" --sp-module step1_rule@R0 --tm-module tm > /dev/null
  $M repeat02 --fasta "$fasta" --workdir "$workdir" --table "$workdir/raw/repeats/repeat02.tsv" --script "$PROJ_ROOT/analysis/cocci_repeats/02_repeat_profile.py" > /dev/null
  $M repeat14 --fasta "$fasta" --workdir "$workdir" --table "$workdir/raw/repeats/repeat14.tsv" --script "$PROJ_ROOT/analysis/cocci_repeats/14_repeat_detect_general.py" > /dev/null
  $M allergen --fasta "$fasta" --workdir "$workdir" --blast "$workdir/raw/allergen/blast.tsv" --allergen-fasta "$S/iuis_fungal_allergens.faa" \
    --blast-version "$(head -1 "$workdir/raw/allergen/version.txt")" --evalue 1 --seg no --max-target-seqs 200 > /dev/null
  rm -rf "$workdir/out"
  cellsurface_sorting_hat --fasta "$fasta" --taxon "$taxon" --taxdump "$T/nodes.dmp" --workdir "$workdir" --out "$workdir/out" > /dev/null
  echo "done $name"
done
echo "all done"
