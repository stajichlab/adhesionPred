#!/usr/bin/env bash
# Structural search of the class 2b confident cores against the experimental PDB.
#
# Two passes:
#   1. 3Di+AA mode (fast, gives an E-value) -- used to decide whether a hit exists at all.
#   2. TM-align mode (--alignment-type 1) -- used to report a TM-score for the top hits.
#
# CTS1 (endochitinase, GH18) is the positive control: it must return strong PDB hits. If it
# does not, the search is broken and nothing else here can be read.
#
# Foldseek is built with AVX2 and dies with "Illegal instruction" on some UCR HPCC login and
# `short` nodes. Run it under srun on epyc (same trap as MMseqs2, see docs/HANDOFF-HPCC.md).
#   srun -p epyc -c 8 --mem 16G -t 60 bash 04_foldseek_pdb.sh <workdir> <foldseek> <dbdir>
set -euo pipefail

WORKDIR="${1:?usage: 04_foldseek_pdb.sh <workdir> <foldseek_binary> <dbdir>}"
FOLDSEEK="${2:?usage: 04_foldseek_pdb.sh <workdir> <foldseek_binary> <dbdir>}"
DBDIR="${3:?usage: 04_foldseek_pdb.sh <workdir> <foldseek_binary> <dbdir>}"

CORE="$WORKDIR/core_pdb"
OUT="$WORKDIR/foldseek"
mkdir -p "$OUT" "$DBDIR"

if [ ! -s "$DBDIR/pdb" ]; then
    echo "downloading the Foldseek PDB database into $DBDIR"
    "$FOLDSEEK" databases PDB "$DBDIR/pdb" "${SCRATCH:-$DBDIR}/fsdb_tmp" --threads 8
fi

FMT="query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,qlen,tlen,lddt,prob,taxname"

TMP="${SCRATCH:?SCRATCH must be set inside a SLURM job}/foldseek_tmp_pdb"
rm -rf "$TMP"; mkdir -p "$TMP"
"$FOLDSEEK" easy-search "$CORE" "$DBDIR/pdb" "$OUT/vs_pdb_3di.tsv" "$TMP" \
    -e 10 --max-seqs 2000 --format-output "$FMT" -v 1
rm -rf "$TMP"

TMP="$SCRATCH/foldseek_tmp_pdb_tm"
rm -rf "$TMP"; mkdir -p "$TMP"
"$FOLDSEEK" easy-search "$CORE" "$DBDIR/pdb" "$OUT/vs_pdb_tmalign.tsv" "$TMP" \
    --alignment-type 1 --tmscore-threshold 0.3 -e 10 --max-seqs 2000 \
    --format-output "$FMT" -v 1
rm -rf "$TMP"

echo "wrote $OUT/vs_pdb_3di.tsv and $OUT/vs_pdb_tmalign.tsv"
