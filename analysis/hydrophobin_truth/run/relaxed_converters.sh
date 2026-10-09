#!/bin/bash
# L7a: turn the domain tables into the hydrophobin_relaxed module for every proteome of RUN_LIST.
# Usage: PROJ_ROOT=... bash relaxed_converters.sh   (needs the csh venv on PATH)
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
HMMER=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PROJ_ROOT/_workdir/sorting_hat/Afum_A1163/raw/pfam/provenance.json'))['hmmer'].split()[2])")
tail -n +2 "$PROJ_ROOT/analysis/hydrophobin_truth/run_list.tsv" | while IFS=$'\t' read -r name fasta taxon kind workdir; do
  cellsurface_sorting_hat_module hydrophobin_relaxed --fasta "$fasta" --workdir "$workdir" --domtbl "$workdir/raw/hydrophobin_relaxed/domtbl.txt" \
    --hmm "$PROJ_ROOT/data/sorting_hat/hydrophobin_relaxed.hmm" --search-option nobias --hmmer-version "$HMMER" | cut -c1-120
done
