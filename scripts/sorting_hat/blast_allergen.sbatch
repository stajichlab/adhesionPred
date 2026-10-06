#!/bin/bash -l
# BLASTP of one proteome against the WHO/IUIS fungal allergen sequences.
# Output: $WORKDIR/raw/allergen/blast.tsv (columns: qseqid sseqid pident length qlen slen bitscore evalue),
#         blast.tsv.done (completion marker: blast version and query count) and version.txt.
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=...,ALLERGEN_FASTA=... blast_allergen.sbatch
# Then:   cellsurface_sorting_hat_module allergen --blast $WORKDIR/raw/allergen/blast.tsv --allergen-fasta $ALLERGEN_FASTA \
#           --blast-version "<version.txt>" --max-target-seqs 200 --evalue 1 --seg no ...
# The --max-target-seqs, --evalue and --seg values of the module command must equal the blastp options below.
# -max_target_seqs 200 is larger than the database (111 sequences: 116 IUIS rows have a sequence, 5 are
# skipped), so no hit is lost to the early-stop behaviour of -max_target_seqs. The cap must exceed the
# database size: the leave-species-out recall test depends on every hit being kept. The job counts the
# database sequences and stops before blastp if the count is not below the cap. -seg no is the blastp default, stated here.
# An empty blast.tsv is refused by the module command, because a finished run with no hit looks like a failed run.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH -J csh_blast
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
: "${ALLERGEN_FASTA:?export ALLERGEN_FASTA (made by build_allergen_fasta) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_blast.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/allergen"
mkdir -p "$TMP" "$OUT"
rm -f "$OUT/blast.tsv.done"  # a marker from an earlier run must not outlive a failed run
CAP=200  # equals -max_target_seqs below
NDB="$(grep -c '^>' "$ALLERGEN_FASTA" || true)"
if [ "$NDB" -ge "$CAP" ]; then
  echo "FATAL: $ALLERGEN_FASTA has $NDB sequences; the -max_target_seqs cap ($CAP) must exceed the database size or hits are lost silently" >&2
  exit 1
fi
module load ncbi-blast/2.14.0+
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
makeblastdb -in "$ALLERGEN_FASTA" -dbtype prot -out "$TMP/allergens" -logfile /dev/null
if ! blastp -query "$TMP/in.fasta" -db "$TMP/allergens" -evalue 1 -max_target_seqs 200 -seg no \
  -outfmt "6 qseqid sseqid pident length qlen slen bitscore evalue" -num_threads "${SLURM_CPUS_ON_NODE:-2}" \
  -out "$TMP/blast.tsv"; then
  echo "FATAL: blastp exited with an error" >&2
  exit 1
fi
[ -e "$TMP/blast.tsv" ] || { echo "FATAL: blastp wrote no output file" >&2; exit 1; }
[ -s "$TMP/blast.tsv" ] || echo "WARNING: blast.tsv is empty (no hit); the module command will refuse it" >&2
VERSION="$(blastp -version | head -1)"
NQ="$(grep -c '^>' "$TMP/in.fasta" || true)"
printf 'blast=%s\nqueries=%s\nhit_lines=%s\n' "$VERSION" "$NQ" "$(wc -l < "$TMP/blast.tsv")" > "$TMP/blast.tsv.done"
cp "$TMP/blast.tsv" "$OUT/.tmp.blast.tsv" && mv "$OUT/.tmp.blast.tsv" "$OUT/blast.tsv"
echo "$VERSION" > "$OUT/.tmp.version.txt" && mv "$OUT/.tmp.version.txt" "$OUT/version.txt"
cp "$TMP/blast.tsv.done" "$OUT/.tmp.blast.tsv.done" && mv "$OUT/.tmp.blast.tsv.done" "$OUT/blast.tsv.done"
