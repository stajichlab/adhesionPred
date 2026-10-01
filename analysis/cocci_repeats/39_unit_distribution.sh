#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 2 --mem 16gb --time 1:00:00
#SBATCH -J unitdist -o logs/39_unit_distribution.%A.log
# Tables from the 37 search. See 39_unit_distribution.py. Run after 37_unit_search.sh.
# Submit from analysis/cocci_repeats:   sbatch 39_unit_distribution.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
/usr/bin/python3.12 39_unit_distribution.py
