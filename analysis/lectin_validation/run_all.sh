#!/bin/bash -l
# Rebuild the lectin validation from the fetched data. Login node, no SLURM.
# Usage: run_all.sh ROOT   (ROOT = absolute path of analysis/lectin_validation)
# Network steps (fetch_*.py, select_assumed.py) are NOT run here. Their outputs are in ROOT/data.
# No BASH_SOURCE is used. Paths come from the argument.
set -euo pipefail
R=${1:?absolute path of analysis/lectin_validation}
SPEC=$R/../lectin_class_spec
module load hmmer/3.4
module load MMseqs2/17-b804f
PY=/usr/bin/python3.12
mkdir -p "$R/out" "$R/results" "$R/controls"
$PY "$R/build_truth.py" --data "$R/data" --curated "$R/curated_positives.tsv" \
  --unilectin-tsv "$SPEC/unilectin_fungal_sequences.tsv" --unilectin-faa "$SPEC/unilectin_fungal_sequences.faa" \
  --shared "$R/../../data/controls/hard-negatives-shared" --out-dir "$R/out"
rm -rf "$R/out/cluster" "$R/out/tmp"
bash "$R/run_cluster_seed.sh" "$R/out/cluster" "$R/out/sequences.faa" "$R/data/seeds.faa.gz" "$R/out/tmp"
$PY "$R/finalize_truth.py" "$R/out/truth_all.tsv" "$R/out/cluster/clu_cluster.tsv" "$R/out/cluster/seed_hits.m8" \
  "$R/out/truth_final.tsv" "$R/out/clusters.tsv"
$PY "$R/scan_models.py" "$R/out/sequences.faa" "$R/out/scan" --cpu 4
$PY "$R/measure.py" --truth "$R/out/truth_final.tsv" --hits "$R/out/scan.hits.tsv" --out-dir "$R/results"
$PY "$R/scan_misses.py" "$R/results/misses.tsv" "$R/out/sequences.faa" "$R/out/scan.models.hmm" "$R/results" --tmpdir "$R/out/tmp"
rm -rf "$R/out/tmp"
$PY "$R/measure.py" --truth "$R/out/truth_final.tsv" --hits "$R/out/scan.hits.tsv" --out-dir "$R/results" --miss-dir "$R/results"
$PY "$R/swissprot_background.py" "$R/out/scan_swissprot.hits.tsv" "$R/data/raw/swissprot_fungi.tsv.gz" \
  "$R/out/truth_final.tsv" "$R/results/swissprot_background"
$PY "$R/margins.py" "$R/out/truth_final.tsv" "$R/out/scan.hits.tsv" "$R/out/scan.models.hmm" "$R/results/ga_margins.tsv"
$PY "$R/d4al26_vs_six.py" "$R/data/raw/uniprot_candidates.json.gz" "$R/../lectin_class_report/f27_six/six.faa" "$R/results/d4al26_vs_six.m8" "$R/out/tmp"
rm -rf "$R/out/tmp" "$R/out/chk"
$PY "$R/check_report_script.py" "$R/../lectin_class_report/lectin_report.py" "$R/out/truth_final.tsv" \
  "$R/out/sequences.faa" "$R/out/scan.hits.tsv" "$R/out/chk" "$R/results/report_script_check.tsv"
rm -rf "$R/out/chk"
$PY "$R/assemble_controls.py" "$R/out" "$R/controls" "$R/rejected_manual.tsv"
