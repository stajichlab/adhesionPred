#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 2 --mem 8gb --time 0:40:00
#SBATCH -J unitbad1 -o logs/38_unit_bad1_control.%A.log
# BAD1 positive control. See 38_unit_bad1_control.py. Run after 37_unit_search.sh.
# Submit from analysis/cocci_repeats:   sbatch 38_unit_bad1_control.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module load hmmer/3.4
/usr/bin/python3.12 38_unit_bad1_control.py
