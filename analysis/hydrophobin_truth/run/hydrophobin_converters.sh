#!/bin/bash
# Convert raw outputs to module tables and run the core command for every proteome in RUN_LIST (task H5 of the
# hydrophobin plan). In-scope proteomes keep their old signalp, tm, repeat and allergen modules; the Pfam modules,
# cys8_pattern and the core run are redone (the family table and the categories changed). Proteomes of kind e1 get
# signalp, tm, pfam and cys8_pattern; their repeat, allergen and antigen modules are not run, so those calls are
# not assessable there.
# Usage: PROJ_ROOT=... RUN_LIST=... bash hydrophobin_converters.sh   (needs the csh venv on PATH)
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${RUN_LIST:?export RUN_LIST}"
S="$PROJ_ROOT/_workdir/sorting_hat"
T="$S/taxdump"
M=cellsurface_sorting_hat_module
tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  echo "== $name"
  PFAMJ="$workdir/raw/pfam/provenance.json"
  SHA=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['sha256'])")
  HMMER=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['hmmer'].split()[2])")
  if [ "$kind" = e1 ]; then
    $M signalp --fasta "$fasta" --workdir "$workdir" --results "$workdir/raw/signalp/prediction_results.txt" --signalp-version "$(paste -sd';' "$workdir/raw/signalp/version.txt")"
    $M tm --fasta "$fasta" --workdir "$workdir" --table "$workdir/raw/tmhmm/tmhmm.tsv"
  fi
  $M pfam --fasta "$fasta" --workdir "$workdir" --domtbl "$workdir/raw/pfam/domtbl.txt" --family-table "$PROJ_ROOT/data/sorting_hat/family_table.tsv" \
    --pfam-release 38.2 --pfam-sha256 "$SHA" --hmmer-version "$HMMER" --sp-module step1_rule@R0 --tm-module tm
  $M cys8 --fasta "$fasta" --workdir "$workdir" --spacing "$PROJ_ROOT/data/sorting_hat/cys8_spacing.yaml" --sp-module step1_rule@R0
  rm -rf "$workdir/out_hyd"
  cellsurface_sorting_hat --fasta "$fasta" --taxon "$taxon" --taxdump "$T/nodes.dmp" --workdir "$workdir" --out "$workdir/out_hyd"
done
echo "all done"
