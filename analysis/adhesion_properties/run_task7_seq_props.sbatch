#!/bin/bash -l
#SBATCH -p highmem
#SBATCH --mem=100G
#SBATCH -c 4
#SBATCH --time=4:00:00
#SBATCH -J adhesion_seq_props
#SBATCH -o logs/task7_seq_props.%j.log
#SBATCH -e logs/task7_seq_props.%j.log

cd /bigdata/stajichlab/jstajich/projects/adhesionPred/.claude/worktrees/adhesion-protein-properties
/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin/python analysis/adhesion_properties/03_compute_sequence_properties.py
