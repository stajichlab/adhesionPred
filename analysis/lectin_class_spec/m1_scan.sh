#!/usr/bin/env bash
# Scan FASTA files with the five M1 models (and the M2, M3 and owned models for context)
# at the Pfam GA cutoff.
# Usage: bash m1_scan.sh OUTDIR FASTA [FASTA ...]
# Needs: hmmer/3.4 on PATH (module load hmmer/3.4). Pfam 38.2 path below.
# No BASH_SOURCE is used. Login node only for small inputs (here under 100k proteins).
set -euo pipefail
if ! command -v hmmsearch > /dev/null 2>&1; then
  source /usr/share/Modules/init/bash 2>/dev/null || true
  module load hmmer/3.4
fi
OUT=${1:?outdir}; shift
PFAM=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm
MODELS="PF07938:M1 PF21691:M1 PF00337:M1 PF07367:M1 PF21595:M1 PF18647:M2 PF22803:M2 PF01419:M2 PF00652:M2 PF14200:M2 PF13385:M2 PF08881:M2 PF01476:M2 PF03388:M3 PF00262:M3 PF11721:M3 PF07915:M3 PF02157:M3 PF07691:own PF10528:own"
mkdir -p "$OUT"
# map accession (no version) to model name once
awk '$1=="NAME"{n=$2} $1=="ACC"{split($2,a,"."); print a[1]"\t"n}' "$PFAM" > "$OUT/acc2name.tsv"
: > "$OUT/models.txt"
: > "$OUT/models.hmm"
for m in $MODELS; do
  acc=${m%%:*}
  name=$(awk -v a="$acc" '$1==a{print $2}' "$OUT/acc2name.tsv")
  hmmfetch "$PFAM" "$name" >> "$OUT/models.hmm"
  echo "$m" >> "$OUT/models.txt"
done
for f in "$@"; do
  b=$(basename "$f"); b=${b%.gz}; b=${b%.faa}; b=${b%.fasta}
  hmmsearch --cut_ga --cpu 8 --noali --domtblout "$OUT/$b.domtbl" "$OUT/models.hmm" "$f" > /dev/null
done
