#!/bin/bash -l
# Sub-threshold look: same models as scan_test_set.py, no --cut_ga, -E 10. Prints best sequence score per
# (protein, model) next to the model GA sequence cutoff. Usage: scan_below_ga.sh MODELS_HMM FASTA OUT_PREFIX
set -euo pipefail
module load hmmer/3.4
hmmsearch -E 10 --noali --tblout "$3.tbl" -o /dev/null "$1" "$2"
awk '/^NAME/{n=$2} /^ACC/{a=$2; sub(/\..*/,"",a)} /^GA /{gsub(";","",$2); print a"\t"n"\t"$2}' "$1" > "$3.ga.tsv"
