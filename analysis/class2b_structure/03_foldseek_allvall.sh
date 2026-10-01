#!/usr/bin/env bash
# All-vs-all structural comparison of the class 2b candidate set.
#
# Question: is class 2b ("small receptor-binding invasins") one structural class, or a
# residual bucket of proteins that only share "the repeat detector misses them"?
#
# Compares the pLDDT>=70 confident cores only (see 02_confidence_profile.py). Uses Foldseek's
# TM-align mode (--alignment-type 1) so the reported number is a TM-score, not a 3Di E-value.
# TM-score >= 0.5 is the conventional same-fold threshold.
#
# Usage (login node is fine, this is seconds of CPU on 10 small structures):
#   bash 03_foldseek_allvall.sh <workdir> <foldseek_binary>
set -euo pipefail

WORKDIR="${1:?usage: 03_foldseek_allvall.sh <workdir> <foldseek_binary>}"
FOLDSEEK="${2:?usage: 03_foldseek_allvall.sh <workdir> <foldseek_binary>}"

CORE="$WORKDIR/core_pdb"
OUT="$WORKDIR/foldseek"
mkdir -p "$OUT"
TMP="${SCRATCH:-$WORKDIR}/foldseek_tmp_allvall"
rm -rf "$TMP"; mkdir -p "$TMP"

# TM-align mode: the third-to-last column is the TM-score normalised by the query.
"$FOLDSEEK" easy-search "$CORE" "$CORE" "$OUT/allvall.tsv" "$TMP" \
    --alignment-type 1 \
    --tmscore-threshold 0.0 \
    -e 1000 \
    --max-seqs 1000 \
    --exhaustive-search 1 \
    --format-output "query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,qlen,tlen,lddt,prob" \
    -v 1

rm -rf "$TMP"
echo "wrote $OUT/allvall.tsv"
