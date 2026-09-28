#!/bin/bash -l
# Recovery job: the interactive SLURM allocation running this session's
# tools is itself a job (mem=16G, 4 cpus). 04_cluster_esmc300m.py's UMAP
# step was OOM-killed by that 16G memory cgroup (confirmed via dmesg:
# "Memory cgroup out of memory: Killed process ... oom_memcg=.../job_28009232").
# ESM-C's compute_umap_projection call operates on the raw 960-dim
# embeddings (not the PCA-reduced 15-dim array -- see cluster.py's
# signature), roughly double ESM2-classifier's 480-dim memory footprint for
# UMAP's nearest-neighbor graph, which is why ESM2's UMAP fit within 16G but
# ESM-C's did not. HDBSCAN's result (tables/cluster_labels_esmc300m.csv) is
# already on disk and intact -- this job only reruns the UMAP step via
# _umap_esmc300m_only.py, using the "highmem" partition with generous
# memory headroom.
#SBATCH -p highmem
#SBATCH --mem=100G
#SBATCH -c 4
#SBATCH --time=4:00:00
#SBATCH -J umap_esmc300m_recovery
#SBATCH -o logs/umap_esmc300m_recovery.%j.log
#SBATCH -e logs/umap_esmc300m_recovery.%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR not set}"

python analysis/embedding_clustering/_umap_esmc300m_only.py
