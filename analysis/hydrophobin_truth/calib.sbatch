#!/bin/bash
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -J hyd_calib
#SBATCH -o /rhome/jstajich/projects/adhesionPred/analysis/hydrophobin_truth/calib.%j.log
set -euo pipefail
module load mmseqs2/17-b804f
cd /rhome/jstajich/projects/adhesionPred
/usr/bin/python3.12 analysis/hydrophobin_truth/make_calibration_inputs.py --dir analysis/hydrophobin_truth
