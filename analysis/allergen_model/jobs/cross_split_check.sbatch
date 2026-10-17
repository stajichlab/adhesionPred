#!/bin/bash -l
# T1.3: inter-split identity check (a measurement; the job removes nothing).
# Every protein is searched against every other protein with MMseqs2. Only pairs whose two IDs come from different splits are kept.
# Split is the ID prefix before the first underscore: train_, val_, test_. The job fails if fewer than two splits are present.
# Identity is fident (fraction, 0-1), the scale of T_S = 0.4.
# --max-seqs 2000: the default 300 is filled by within-split hits on a set of ~11,000 proteins, which can push crossing hits out.
# Each crossing pair appears in both directions; only query < target is written, so each pair is counted once.
# Input:  $SPLIT_FASTA (all split members, prefixed train_/val_/test_), $ARM
# Output: $OUTDIR/cross_split_pairs.tsv.zst  (query, target, fident, alnlen)
#         $OUTDIR/cross_split_max.txt        (maximum fident across splits; 0 if there are no crossing pairs)
#         $OUTDIR/cross_split.done           (written last; split sizes, counts, input hash)
# Submit: sbatch --export=ALL,SPLIT_FASTA=...,OUTDIR=...,ARM=... cross_split_check.sbatch
# Partition: stajichlab, 16 cpus, 32 GB, 4 h.
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -J ap_xsplit
set -euo pipefail
: "${SPLIT_FASTA:?export SPLIT_FASTA (all split members, prefixed by split) before sbatch}"
: "${OUTDIR:?export OUTDIR before sbatch}"
: "${ARM:?export ARM before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/ap_xsplit.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$OUTDIR"
rm -f "$OUTDIR/cross_split.done"
command -v zstd >/dev/null || { echo "FATAL: zstd is not on PATH" >&2; exit 1; }
MMSEQS_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs
[ -x "$MMSEQS_BIN" ] || { echo "FATAL: MMseqs2 binary not found at $MMSEQS_BIN" >&2; exit 1; }  # non-AVX2 build: the AVX2 build crashes with SIGILL on some nodes (c03)
[ -s "$SPLIT_FASTA" ] || { echo "FATAL: $SPLIT_FASTA is missing or empty" >&2; exit 1; }
grep '^>' "$SPLIT_FASTA" | sed 's/^>//; s/_.*//' | sort | uniq -c > "$TMP/split_counts.txt"
NSPLIT="$(wc -l < "$TMP/split_counts.txt")"
[ "$NSPLIT" -ge 2 ] || { echo "FATAL: fewer than two split prefixes in $SPLIT_FASTA; the check would be vacuous" >&2; exit 1; }
for s in train val test; do
  grep -q " $s\$" "$TMP/split_counts.txt" || { echo "FATAL: split prefix '$s' missing from $SPLIT_FASTA" >&2; exit 1; }
done
"$MMSEQS_BIN" createdb "$SPLIT_FASTA" "$TMP/db" > "$TMP/createdb.log" 2>&1
"$MMSEQS_BIN" search "$TMP/db" "$TMP/db" "$TMP/res" "$TMP/tmp" \
  --min-seq-id 0 -c 0 --cov-mode 0 -s 7.5 -e 1000 --max-seqs 2000 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/search.log" 2>&1
"$MMSEQS_BIN" convertalis "$TMP/db" "$TMP/db" "$TMP/res" "$TMP/all.tsv" \
  --format-output "query,target,fident,alnlen" --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/convert.log" 2>&1
[ -e "$TMP/all.tsv" ] || { echo "FATAL: convertalis wrote no file" >&2; exit 1; }
# Crossing pairs: different split prefixes; keep query < target so each unordered pair is written once.
awk -F'\t' 'BEGIN{OFS="\t"} {split($1,a,"_"); split($2,b,"_"); if (a[1] != b[1] && $1 < $2) print $1, $2, $3, $4}' "$TMP/all.tsv" > "$TMP/cross.tsv"
MAX="$(cut -f3 "$TMP/cross.tsv" | sort -g | tail -1 || true)"
MAX="${MAX:-0}"
{ printf 'query\ttarget\tfident\talnlen\n'; cat "$TMP/cross.tsv"; } | zstd -q -T0 -f -o "$TMP/cross_split_pairs.tsv.zst"
cp "$TMP/cross_split_pairs.tsv.zst" "$OUTDIR/.tmp.cross_split_pairs.tsv.zst" && mv "$OUTDIR/.tmp.cross_split_pairs.tsv.zst" "$OUTDIR/cross_split_pairs.tsv.zst"
echo "$MAX" > "$OUTDIR/.tmp.cross_split_max.txt" && mv "$OUTDIR/.tmp.cross_split_max.txt" "$OUTDIR/cross_split_max.txt"
{
  echo "mmseqs=$("$MMSEQS_BIN" version)"
  echo "command=search -s 7.5 -e 1000 --max-seqs 2000 --min-seq-id 0 -c 0"
  echo "arm=$ARM"
  echo "split_counts:"; sed 's/^/  /' "$TMP/split_counts.txt"
  echo "max_fident_cross_split=$MAX"
  echo "crossing_pairs_unordered=$(wc -l < "$TMP/cross.tsv")"
  echo "input_sha256=$(sha256sum "$SPLIT_FASTA" | cut -d' ' -f1)"
  echo "date=$(date +%F)"
} > "$TMP/cross_split.done"
cp "$TMP/cross_split.done" "$OUTDIR/.tmp.cross_split.done" && mv "$OUTDIR/.tmp.cross_split.done" "$OUTDIR/cross_split.done"
