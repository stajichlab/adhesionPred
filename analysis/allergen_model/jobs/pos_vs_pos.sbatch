#!/bin/bash -l
# B6 input: all-vs-all identity among the positives (MMseqs2 search, fident fraction), so each positive's nearest
# training positive can be read. Input: $POSITIVES (P_ IDs), $OUTDIR. Output: pos_vs_pos.tsv.zst, pos_vs_pos.done.
# Partition: stajichlab, 8 cpus, 16 GB, 1 h.
#SBATCH -p stajichlab
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -J ap_posvpos
set -euo pipefail
: "${POSITIVES:?export POSITIVES (positives FASTA, P_ IDs) before sbatch}"
: "${OUTDIR:?export OUTDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/ap_posvpos.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$OUTDIR"
rm -f "$OUTDIR/pos_vs_pos.done"
command -v zstd >/dev/null || { echo "FATAL: zstd is not on PATH" >&2; exit 1; }
MMSEQS_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs
[ -x "$MMSEQS_BIN" ] || { echo "FATAL: MMseqs2 binary not found" >&2; exit 1; }
NP="$(grep -c '^>P_' "$POSITIVES" || true)"
[ "$NP" -gt 1 ] || { echo "FATAL: $POSITIVES has fewer than two P_ positives" >&2; exit 1; }
"$MMSEQS_BIN" createdb "$POSITIVES" "$TMP/db" > "$TMP/createdb.log" 2>&1
"$MMSEQS_BIN" search "$TMP/db" "$TMP/db" "$TMP/res" "$TMP/tmp" \
  --min-seq-id 0 -c 0 --cov-mode 0 -s 7.5 -e 1000 --max-seqs 2000 --threads "${SLURM_CPUS_ON_NODE:-8}" > "$TMP/search.log" 2>&1
"$MMSEQS_BIN" convertalis "$TMP/db" "$TMP/db" "$TMP/res" "$TMP/pairs.tsv" \
  --format-output "query,target,fident,alnlen,qcov,tcov" --threads "${SLURM_CPUS_ON_NODE:-8}" > "$TMP/convert.log" 2>&1
[ -e "$TMP/pairs.tsv" ] || { echo "FATAL: convertalis wrote no file" >&2; exit 1; }
{ printf 'query\ttarget\tfident\talnlen\tqcov\ttcov\n'; cat "$TMP/pairs.tsv"; } | zstd -q -T0 -f -o "$TMP/pos_vs_pos.tsv.zst"
cp "$TMP/pos_vs_pos.tsv.zst" "$OUTDIR/.tmp.pos.zst" && mv "$OUTDIR/.tmp.pos.zst" "$OUTDIR/pos_vs_pos.tsv.zst"
{
  echo "mmseqs=$("$MMSEQS_BIN" version)"
  echo "positives=$NP"
  echo "input_sha256=$(sha256sum "$POSITIVES" | cut -d' ' -f1)"
  echo "date=$(date +%F)"
} > "$TMP/pos_vs_pos.done"
cp "$TMP/pos_vs_pos.done" "$OUTDIR/.tmp.done" && mv "$OUTDIR/.tmp.done" "$OUTDIR/pos_vs_pos.done"
