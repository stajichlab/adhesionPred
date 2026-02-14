#!/usr/bin/bash -l
module load ncbi_datasets

datasets download genome taxon "Saccharomyces cerevisiae" --include protein
unzip ncbi_dataset.zip
cat ncbi_dataset/data/*/*.faa  > scer_proteins.faa

