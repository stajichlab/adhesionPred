#!/bin/bash -l
# Foldseek search of the PRA3 full-length AlphaFold model (E9CRM7) against AlphaFold DB subsets.
# Queries: PRA3_full core (67 aa) and full model, plus two controls (Ag2_PRA core: CFEM; RodA core:
# hydrophobin). Each query runs in 3Di+AA mode and in TM-align mode.
# Environment (no script-relative paths; this runs under SLURM):
#   WORK     work directory of the full-length run (has core_pdb/ and full_pdb/)
#   FOLDSEEK foldseek binary
#   AFDB     prefix of the Foldseek database (for example .../afdb/afdb_swissprot)
#   TAG      short name of the database for the output file names
# Submit: sbatch -p highclock -c 16 --mem=<G> -t <time> --export=ALL,WORK=...,FOLDSEEK=...,AFDB=...,TAG=... 05_afdb_search.sh
#SBATCH -p highclock
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -J pra3_afdb
set -euo pipefail
: "${WORK:?}" "${FOLDSEEK:?}" "${AFDB:?}" "${TAG:?}"
TMP="${SCRATCH:?run as a SLURM job}/pra3_afdb_$TAG"
rm -rf "$TMP"; mkdir -p "$TMP/q_core" "$TMP/q_full" "$TMP/t1" "$TMP/t2" "$TMP/t3" "$TMP/t4" "$WORK/foldseek"
ln -s "$WORK/core_pdb/PRA3_full.pdb" "$TMP/q_core/PRA3_full_core.pdb"
ln -s "$WORK/core_pdb/Ag2_PRA.pdb" "$TMP/q_core/Ag2_PRA_core.pdb"
ln -s "$WORK/core_pdb/RodA.pdb" "$TMP/q_core/RodA_core.pdb"
ln -s "$WORK/full_pdb/PRA3_full.pdb" "$TMP/q_full/PRA3_full_model.pdb"
if [ -s "$WORK/core_pdb/ARB_05178.pdb" ]; then  # present when the work directory was built with candidates_d4ali1.tsv
  ln -s "$WORK/core_pdb/ARB_05178.pdb" "$TMP/q_core/ARB_05178_core.pdb"
  ln -s "$WORK/full_pdb/ARB_05178.pdb" "$TMP/q_full/ARB_05178_model.pdb"
fi
FMT="query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,qlen,tlen,lddt,prob,taxname"
T="${SLURM_CPUS_PER_TASK:-8}"
"$FOLDSEEK" easy-search "$TMP/q_core" "$AFDB" "$WORK/foldseek/vs_${TAG}_core_3di.tsv" "$TMP/t1" \
    -e 10 --max-seqs 2000 --threads "$T" --format-output "$FMT" -v 1
"$FOLDSEEK" easy-search "$TMP/q_core" "$AFDB" "$WORK/foldseek/vs_${TAG}_core_tmalign.tsv" "$TMP/t2" \
    --alignment-type 1 --tmscore-threshold 0.3 -e 10 --max-seqs 2000 --threads "$T" \
    --format-output "$FMT" -v 1
"$FOLDSEEK" easy-search "$TMP/q_full" "$AFDB" "$WORK/foldseek/vs_${TAG}_full_3di.tsv" "$TMP/t3" \
    -e 10 --max-seqs 2000 --threads "$T" --format-output "$FMT" -v 1
"$FOLDSEEK" easy-search "$TMP/q_full" "$AFDB" "$WORK/foldseek/vs_${TAG}_full_tmalign.tsv" "$TMP/t4" \
    --alignment-type 1 --tmscore-threshold 0.3 -e 10 --max-seqs 2000 --threads "$T" \
    --format-output "$FMT" -v 1
echo "done $TAG"
