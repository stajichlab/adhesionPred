#!/usr/bin/env bash
# Targeted pairwise TM-align against experimental reference structures.
#
# Why this step exists: the whole-PDB Foldseek search (04) returns NO hit to 4Y7S (the CFEM
# protein Csa2) for Ag2/PRA or PRA2, even in TM-align mode -- while the direct all-vs-all (03)
# scores Ag2/PRA against the Csa2 AlphaFold model at TM 0.81 / LDDT 0.77. The 3Di prefilter
# drops the true neighbour of a ~64-residue query. This step re-runs the comparison with
# --exhaustive-search 1, which skips the prefilter, against the experimental structures.
#
#   srun/sbatch on epyc (AVX2): bash 06_targeted_pairwise.sh <workdir> <foldseek>
set -euo pipefail

WORKDIR="${1:?usage: 06_targeted_pairwise.sh <workdir> <foldseek_binary>}"
FOLDSEEK="${2:?usage: 06_targeted_pairwise.sh <workdir> <foldseek_binary>}"

REF="$WORKDIR/ref_pdb"
OUT="$WORKDIR/foldseek"
mkdir -p "$REF" "$OUT"

# Experimental references, one per fold hypothesis under test.
#   4Y7S  CFEM protein Csa2 (Candida albicans)          -- the CFEM fold
#   5FID  elicitor MoHrip2 (Magnaporthe oryzae)         -- CalA's top PDB hit
#   1LL7  endochitinase CTS1 (Coccidioides immitis)     -- pipeline positive control
#   6GCJ  hydrophobin RodA (Aspergillus fumigatus)      -- class 2c control
for pdbid in 4y7s 5fid 1ll7 6gcj; do
    if [ ! -s "$REF/${pdbid}.pdb" ]; then
        curl -sL "https://files.rcsb.org/download/${pdbid^^}.pdb" -o "$REF/${pdbid}.pdb"
    fi
done

TMP="${SCRATCH:-$WORKDIR}/foldseek_tmp_ref"
rm -rf "$TMP"; mkdir -p "$TMP"
"$FOLDSEEK" easy-search "$WORKDIR/core_pdb" "$REF" "$OUT/vs_ref_tmalign.tsv" "$TMP" \
    --alignment-type 1 --tmscore-threshold 0.0 -e inf --exhaustive-search 1 \
    --format-output "query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,qlen,tlen,lddt,prob" \
    -v 1
rm -rf "$TMP"
echo "wrote $OUT/vs_ref_tmalign.tsv"
