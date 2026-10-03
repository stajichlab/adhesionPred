#!/usr/bin/bash
# Run all depth-vs-SOWgp statistics. Needs the mosdepth outputs from the Genotyping project
# (pipeline/11b_mosdepth_gene.sh, 11c_mosdepth_regions.sh). Run from this directory.
set -euo pipefail
PY=/usr/bin/python3.12
$PY 41_sowgp_status_depth.py     | tee sowgp_status_depth.log
$PY 42_sowgp_gene_depth_vs_units.py | tee sowgp_gene_depth_vs_units.log
$PY 40_sowgp_array_depth.py       | tee sowgp_array_depth.log
