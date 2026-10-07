#!/bin/bash
# Convert the raw tool outputs of one proteome to module tables and run the core command
# (Plan 2, Task 12 step 3). Run on a login or compute node after the five module jobs finished.
# Usage: PROJ_ROOT=... PROT=name FASTA=path TAXON=id [STATUS=1] bash run_converters.sh
# Needs the csh venv on PATH (cellsurface_sorting_hat, cellsurface_sorting_hat_module).
# STATUS=1 also writes the Phase C status source for step1_rule@R0 (only for the Af293 UniProt run).
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}" "${PROT:?export PROT}" "${FASTA:?export FASTA}" "${TAXON:?export TAXON}"
S="$PROJ_ROOT/_workdir/sorting_hat"
WORKDIR="$S/$PROT"
T="$S/taxdump"
M=cellsurface_sorting_hat_module
FQ="$S/$PROT.faa"
case "$FASTA" in *.gz) zcat "$FASTA" > "$FQ" ;; *) cp "$FASTA" "$FQ" ;; esac
PFAMJ="$WORKDIR/raw/pfam/provenance.json"
SHA=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['sha256'])")
HMMER=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['hmmer'].split()[2])")
$M signalp --fasta "$FQ" --workdir "$WORKDIR" --results "$WORKDIR/raw/signalp/prediction_results.txt" --signalp-version "$(paste -sd';' "$WORKDIR/raw/signalp/version.txt")"
$M tm --fasta "$FQ" --workdir "$WORKDIR" --table "$WORKDIR/raw/tmhmm/tmhmm.tsv"
$M pfam --fasta "$FQ" --workdir "$WORKDIR" --domtbl "$WORKDIR/raw/pfam/domtbl.txt" --family-table "$PROJ_ROOT/data/sorting_hat/family_table.tsv" \
  --pfam-release 38.2 --pfam-sha256 "$SHA" --hmmer-version "$HMMER" --sp-module step1_rule@R0 --tm-module tm
$M repeat02 --fasta "$FQ" --workdir "$WORKDIR" --table "$WORKDIR/raw/repeats/repeat02.tsv" --script "$PROJ_ROOT/analysis/cocci_repeats/02_repeat_profile.py"
$M repeat14 --fasta "$FQ" --workdir "$WORKDIR" --table "$WORKDIR/raw/repeats/repeat14.tsv" --script "$PROJ_ROOT/analysis/cocci_repeats/14_repeat_detect_general.py"
$M allergen --fasta "$FQ" --workdir "$WORKDIR" --blast "$WORKDIR/raw/allergen/blast.tsv" --allergen-fasta "$S/iuis_fungal_allergens.faa" \
  --blast-version "$(head -1 "$WORKDIR/raw/allergen/version.txt")" --evalue 1 --seg no --max-target-seqs 200
if [ "${STATUS:-0}" = 1 ]; then
  cellsurface_sorting_hat_calibrate phasec --workdir "$WORKDIR" --metrics "$PROJ_ROOT/_workdir/step1_compare/phasec/metrics.json" \
    --clusters "$PROJ_ROOT/_workdir/step1_compare/phasec/clusters.tsv.gz" --eval-table "$PROJ_ROOT/_workdir/step1_compare/phasec/eval_table.tsv.gz" \
    --phasec-signalp-module signalp/6-gpu --phasec-signalp-mode fast --set-species "$PROJ_ROOT/data/sorting_hat/phasec_set_species.tsv" \
    --names-dmp "$T/names.dmp" --nodes-dmp "$T/nodes.dmp"
fi
# RS needs the lookup modules first (antigen, expression, cys); run them before this script and it picks them up.
cellsurface_sorting_hat --fasta "$FQ" --taxon "$TAXON" --taxdump "$T/nodes.dmp" --workdir "$WORKDIR" --out "$WORKDIR/out"
echo "done $PROT"
