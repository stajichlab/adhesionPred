#!/bin/bash -l
# PLM embedding pilot on the reference proteomes: throughput, memory and GPU
# utilization per model, to size the NRP bulk run (issue #16). Submit from this
# directory after setup_env.sh:   sbatch run_pilot.sbatch
#SBATCH -p exfab
#SBATCH -A exfab
#SBATCH --gres=gpu:ada6000:1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH -J adh_plm_pilot
#SBATCH -o logs/pilot.%j.log
#SBATCH -e logs/pilot.%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export HF_HUB_OFFLINE=1
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv,noheader -l 10 \
  > logs/gpu_util.$SLURM_JOB_ID.csv &
SMI=$!
.venv/bin/python embed_pilot.py --genomes genomes.tsv --out out ${PILOT_ARGS:-}
kill $SMI
