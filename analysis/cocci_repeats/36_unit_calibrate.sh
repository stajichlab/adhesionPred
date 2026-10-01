#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 8 --mem 8gb --time 0:40:00
#SBATCH -J unitcal -o logs/36_unit_calibrate.%A.log
# Null floor and sensitivity of the SOWgp unit HMM. See 36_unit_calibrate.py.
# Submit from analysis/cocci_repeats:   sbatch 36_unit_calibrate.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module load hmmer/3.4
/usr/bin/python3.12 36_unit_calibrate.py
