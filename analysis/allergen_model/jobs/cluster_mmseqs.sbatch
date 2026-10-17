#!/bin/bash -l
# T0.4: cluster the allergen-predictor positives and weak-negative candidates together.
# MMseqs2 easy-cluster, --min-seq-id 0.3 -c 0.5 --cov-mode 0 (COMMON-RULES rule 4), -s 7.5 (higher sensitivity than the default 4.0).
# Input:  $SEQS (FASTA; IDs prefixed P_ (positive) or N_ (candidate) so the source is recoverable), $ARM (arm label, e.g. shared)
# Output: $OUTDIR/clusters.tsv        accession<TAB>cluster_id   (COMMON-RULES schema)
#         $OUTDIR/cluster_params.txt  parameters and input hash
#         $OUTDIR/cluster.done        written last
# Submit: sbatch --export=ALL,SEQS=...,OUTDIR=...,ARM=shared cluster_mmseqs.sbatch
# Partition: stajichlab, 16 cpus, 32 GB, 6 h. Timed pilot on 1,000 candidates first (plan T0.4).
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH -J ap_cluster
set -euo pipefail
: "${SEQS:?export SEQS (positives + candidates FASTA) before sbatch}"
: "${OUTDIR:?export OUTDIR before sbatch}"
: "${ARM:?export ARM (arm label) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/ap_cluster.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$OUTDIR"
rm -f "$OUTDIR/cluster.done"   # a marker from an earlier run must not outlive a failed run
MMSEQS_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs
[ -x "$MMSEQS_BIN" ] || { echo "FATAL: MMseqs2 binary not found at $MMSEQS_BIN" >&2; exit 1; }  # non-AVX2 build: the AVX2 build crashes with SIGILL on some nodes (c03)
[ -s "$SEQS" ] || { echo "FATAL: $SEQS is missing or empty" >&2; exit 1; }
case "$SEQS" in *.gz) zcat "$SEQS" > "$TMP/in.faa" ;; *) cp "$SEQS" "$TMP/in.faa" ;; esac
NP="$(grep -c '^>P_' "$TMP/in.faa" || true)"
NN="$(grep -c '^>N_' "$TMP/in.faa" || true)"
NSEQ="$(grep -c '^>' "$TMP/in.faa" || true)"
[ "$NP" -gt 0 ] || { echo "FATAL: no P_ positives in $SEQS; the ID prefix scheme is wrong" >&2; exit 1; }
[ "$NN" -gt 0 ] || { echo "FATAL: no N_ candidates in $SEQS; the ID prefix scheme is wrong" >&2; exit 1; }
[ "$NSEQ" -eq $((NP + NN)) ] || { echo "FATAL: $NSEQ records but $((NP + NN)) prefixed; some IDs lack P_ or N_" >&2; exit 1; }
"$MMSEQS_BIN" easy-cluster "$TMP/in.faa" "$TMP/clu" "$TMP/tmp" \
  --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/cluster.log" 2>&1
[ -s "$TMP/clu_cluster.tsv" ] || { echo "FATAL: easy-cluster wrote no cluster table" >&2; exit 1; }
# easy-cluster writes representative<TAB>member. The COMMON-RULES schema is accession<TAB>cluster_id, with the representative as cluster_id.
awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1}' "$TMP/clu_cluster.tsv" > "$TMP/clusters.tsv"
NCL="$(cut -f2 "$TMP/clusters.tsv" | sort -u | wc -l)"
NROW="$(wc -l < "$TMP/clusters.tsv")"
[ "$NROW" -eq "$NSEQ" ] || { echo "FATAL: $NROW cluster rows for $NSEQ sequences" >&2; exit 1; }
INHASH="$(sha256sum "$SEQS" | cut -d' ' -f1)"
COMMIT="$(git -C /bigdata/stajichlab/jstajich/projects/adhesionPred rev-parse --short HEAD 2>/dev/null || echo unknown)"
{
  echo "mmseqs=$("$MMSEQS_BIN" version)"
  echo "command=easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5"
  echo "arm=$ARM"
  echo "input=$SEQS"
  echo "input_sha256=$INHASH"
  echo "sequences=$NSEQ positives=$NP candidates=$NN clusters=$NCL"
  echo "repo_commit=$COMMIT"
  echo "date=$(date +%F)"
} > "$TMP/cluster_params.txt"
cp "$TMP/clusters.tsv" "$OUTDIR/.tmp.clusters.tsv" && mv "$OUTDIR/.tmp.clusters.tsv" "$OUTDIR/clusters.tsv"
cp "$TMP/cluster_params.txt" "$OUTDIR/.tmp.cluster_params.txt" && mv "$OUTDIR/.tmp.cluster_params.txt" "$OUTDIR/cluster_params.txt"
cp "$TMP/cluster_params.txt" "$OUTDIR/.tmp.cluster.done" && mv "$OUTDIR/.tmp.cluster.done" "$OUTDIR/cluster.done"
