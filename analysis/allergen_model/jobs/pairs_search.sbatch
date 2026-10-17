#!/bin/bash -l
# T0.5: pairwise identity of every weak-negative candidate against every positive (MMseqs2 search).
# This table is the only input to the build gate (T0.6) and to pairing (T1.2).
# Identity is fident (fraction, 0-1), the scale the spec thresholds use (T_C 0.6 and 0.7). pident (percent) is not used.
# No coverage or alignment-length condition is applied here. The condition is an open owner decision (plan O7); the table
# carries qcov, tcov and alnlen so any condition can be applied later without rerunning the search.
# Input:  $QUERY (candidate negatives FASTA, N_ IDs), $TARGET (positives FASTA, P_ IDs), $ARM
# Output: $OUTDIR/pairs_candidates_vs_positives.tsv.zst   (query, target, fident, alnlen, qcov, tcov, evalue, bits)
#         $OUTDIR/pairs.done  (written last; counts and input hashes)
# Submit: sbatch --export=ALL,QUERY=...,TARGET=...,OUTDIR=...,ARM=shared pairs_search.sbatch
# Partition: stajichlab, 16 cpus, 32 GB, 6 h. Final time set by the T0.4 pilot.
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH -J ap_pairs
set -euo pipefail
: "${QUERY:?export QUERY (candidate negatives FASTA) before sbatch}"
: "${TARGET:?export TARGET (positives FASTA) before sbatch}"
: "${OUTDIR:?export OUTDIR before sbatch}"
: "${ARM:?export ARM (arm label) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/ap_pairs.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$OUTDIR"
rm -f "$OUTDIR/pairs.done"
command -v zstd >/dev/null || { echo "FATAL: zstd is not on PATH; load it before submitting, not after the search" >&2; exit 1; }
MMSEQS_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs
[ -x "$MMSEQS_BIN" ] || { echo "FATAL: MMseqs2 binary not found at $MMSEQS_BIN" >&2; exit 1; }  # non-AVX2 build: the AVX2 build crashes with SIGILL on some nodes (c03)
NQ="$(grep -c '^>N_' "$QUERY" || true)"
NT="$(grep -c '^>P_' "$TARGET" || true)"
[ "$NQ" -gt 0 ] || { echo "FATAL: $QUERY has no N_ candidates" >&2; exit 1; }
[ "$NT" -gt 0 ] || { echo "FATAL: $TARGET has no P_ positives" >&2; exit 1; }
"$MMSEQS_BIN" createdb "$QUERY" "$TMP/qdb" > "$TMP/q.log" 2>&1
"$MMSEQS_BIN" createdb "$TARGET" "$TMP/tdb" > "$TMP/t.log" 2>&1
# -e 1000 and --max-seqs keep every candidate-positive hit the prefilter returns for all positives.
"$MMSEQS_BIN" search "$TMP/qdb" "$TMP/tdb" "$TMP/res" "$TMP/tmp" \
  --min-seq-id 0 -c 0 --cov-mode 0 -s 7.5 -e 1000 --max-seqs 2000 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/search.log" 2>&1
"$MMSEQS_BIN" convertalis "$TMP/qdb" "$TMP/tdb" "$TMP/res" "$TMP/pairs.tsv" \
  --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits" --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/convert.log" 2>&1
[ -e "$TMP/pairs.tsv" ] || { echo "FATAL: convertalis wrote no file" >&2; exit 1; }
[ -s "$TMP/pairs.tsv" ] || echo "WARNING: no candidate-positive hits; the build gate will report zero pairable candidates" >&2
{ printf 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits\n'; cat "$TMP/pairs.tsv"; } > "$TMP/pairs_h.tsv"
zstd -q -T0 -f "$TMP/pairs_h.tsv" -o "$TMP/pairs.tsv.zst"
cp "$TMP/pairs.tsv.zst" "$OUTDIR/.tmp.pairs.tsv.zst" && mv "$OUTDIR/.tmp.pairs.tsv.zst" "$OUTDIR/pairs_candidates_vs_positives.tsv.zst"
{
  echo "mmseqs=$("$MMSEQS_BIN" version)"
  echo "command=search -s 7.5 -e 1000 --max-seqs 2000 --min-seq-id 0 -c 0"
  echo "arm=$ARM"
  echo "queries_candidates=$NQ targets_positives=$NT"
  echo "query_sha256=$(sha256sum "$QUERY" | cut -d' ' -f1)"
  echo "target_sha256=$(sha256sum "$TARGET" | cut -d' ' -f1)"
  echo "pair_rows=$(wc -l < "$TMP/pairs.tsv")"
  echo "date=$(date +%F)"
} > "$TMP/pairs.done"
cp "$TMP/pairs.done" "$OUTDIR/.tmp.pairs.done" && mv "$OUTDIR/.tmp.pairs.done" "$OUTDIR/pairs.done"
