#!/bin/bash -l
# T0.4 timed pilot: runs the clustering and the candidate-vs-positive search on a random sample of candidates, times each
# step, and writes a projection to the full candidate set. The projection sets the time limits of cluster_mmseqs.sbatch and
# pairs_search.sbatch (plan T0.4). The pilot does not write the production tables and does not feed the build gate.
# Input:  $CANDIDATES (full candidate FASTA, N_ IDs), $POSITIVES (positives FASTA, P_ IDs), $N_SAMPLE (default 1000), $OUTDIR
# Output: $OUTDIR/pilot.tsv  (step, seconds, sequences, pair_rows, max_rss_kb or NA)
#         $OUTDIR/pilot_projection.txt (projected seconds for the full candidate set, linear in candidates)
#         $OUTDIR/pilot.done (written last)
# Submit: sbatch --export=ALL,CANDIDATES=...,POSITIVES=...,OUTDIR=... pilot_mmseqs.sbatch
# Partition: stajichlab, 16 cpus, 32 GB, 1 h (the pilot itself is small; the 1 h limit is an upper bound).
#SBATCH -p stajichlab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH -J ap_pilot
set -euo pipefail
: "${CANDIDATES:?export CANDIDATES (full candidate FASTA) before sbatch}"
: "${POSITIVES:?export POSITIVES (positives FASTA) before sbatch}"
: "${OUTDIR:?export OUTDIR before sbatch}"
N_SAMPLE="${N_SAMPLE:-1000}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/ap_pilot.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$OUTDIR"
rm -f "$OUTDIR/pilot.done"
command -v zstd >/dev/null || { echo "FATAL: zstd is not on PATH" >&2; exit 1; }
MMSEQS_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs
[ -x "$MMSEQS_BIN" ] || { echo "FATAL: MMseqs2 binary not found at $MMSEQS_BIN" >&2; exit 1; }  # non-AVX2 build: the AVX2 build crashes with SIGILL on some nodes (c03)
[ -s "$CANDIDATES" ] || { echo "FATAL: $CANDIDATES is missing or empty" >&2; exit 1; }
[ -s "$POSITIVES" ] || { echo "FATAL: $POSITIVES is missing or empty" >&2; exit 1; }
NCAND="$(grep -c '^>N_' "$CANDIDATES" || true)"
[ "$NCAND" -gt 0 ] || { echo "FATAL: $CANDIDATES has no N_ candidates" >&2; exit 1; }
NS="$N_SAMPLE"
[ "$NS" -le "$NCAND" ] || NS="$NCAND"

# Random sample of candidates with a fixed seed, so the pilot is repeatable. Each record gets a random key;
# the records with the smallest keys are kept, and the sample is written back in FASTA form.
# Each record is put on one line (newlines become \001) so that sort works on whole records.
# sort writes to a file first: piping sort into head ends sort with SIGPIPE, which pipefail reports as a failure (exit 141).
awk 'BEGIN{RS=">"; ORS=""; srand(20261010)} NR>1 {gsub(/\n/, "\001"); printf "%.10f\t>%s\n", rand(), $0}' "$CANDIDATES" > "$TMP/keyed.txt"
sort -k1,1g "$TMP/keyed.txt" > "$TMP/sorted.txt"
head -n "$NS" "$TMP/sorted.txt" | cut -f2- | tr '\001' '\n' | sed '/^$/d' > "$TMP/sample.faa"
SAMPLE_N="$(grep -c '^>' "$TMP/sample.faa" || true)"
[ "$SAMPLE_N" -gt 0 ] || { echo "FATAL: the sample is empty" >&2; exit 1; }

TIMER=""
command -v /usr/bin/time >/dev/null && TIMER="/usr/bin/time -f %M"
now() { date +%s.%N; }

# Step 1: clustering of the sample plus all positives.
cat "$POSITIVES" "$TMP/sample.faa" > "$TMP/cluster_in.faa"
t0="$(now)"
if [ -n "$TIMER" ]; then
  $TIMER -o "$TMP/rss1.txt" "$MMSEQS_BIN" easy-cluster "$TMP/cluster_in.faa" "$TMP/clu" "$TMP/tmp1" \
    --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/c.log" 2>&1
  RSS1="$(tail -1 "$TMP/rss1.txt")"
else
  "$MMSEQS_BIN" easy-cluster "$TMP/cluster_in.faa" "$TMP/clu" "$TMP/tmp1" \
    --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/c.log" 2>&1
  RSS1="NA"
fi
t1="$(now)"
NCLU_SEQS="$(grep -c '^>' "$TMP/cluster_in.faa" || true)"

# Step 2: candidate-vs-positive search for the sample.
"$MMSEQS_BIN" createdb "$TMP/sample.faa" "$TMP/qdb" > /dev/null 2>&1
"$MMSEQS_BIN" createdb "$POSITIVES" "$TMP/tdb" > /dev/null 2>&1
t2="$(now)"
if [ -n "$TIMER" ]; then
  $TIMER -o "$TMP/rss2.txt" "$MMSEQS_BIN" search "$TMP/qdb" "$TMP/tdb" "$TMP/res" "$TMP/tmp2" \
    --min-seq-id 0 -c 0 --cov-mode 0 -s 7.5 -e 1000 --max-seqs 2000 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/s.log" 2>&1
  RSS2="$(tail -1 "$TMP/rss2.txt")"
else
  "$MMSEQS_BIN" search "$TMP/qdb" "$TMP/tdb" "$TMP/res" "$TMP/tmp2" \
    --min-seq-id 0 -c 0 --cov-mode 0 -s 7.5 -e 1000 --max-seqs 2000 --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/s.log" 2>&1
  RSS2="NA"
fi
"$MMSEQS_BIN" convertalis "$TMP/qdb" "$TMP/tdb" "$TMP/res" "$TMP/pairs.tsv" \
  --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits" --threads "${SLURM_CPUS_ON_NODE:-16}" > "$TMP/cv.log" 2>&1
t3="$(now)"
PAIR_ROWS="$(wc -l < "$TMP/pairs.tsv" 2>/dev/null || echo 0)"

CLU_S="$(echo "$t1 - $t0" | bc -l)"
SRCH_S="$(echo "$t3 - $t2" | bc -l)"
{
  printf 'step\tseconds\tsequences\tpair_rows\tmax_rss_kb\n'
  printf 'cluster\t%.1f\t%s\tNA\t%s\n' "$CLU_S" "$NCLU_SEQS" "$RSS1"
  printf 'search\t%.1f\t%s\t%s\t%s\n' "$SRCH_S" "$SAMPLE_N" "$PAIR_ROWS" "$RSS2"
} > "$TMP/pilot.tsv"

# Projection: the search scales with candidates (each candidate is one query). Clustering is reported as measured for the
# sample size, because MMseqs2 clustering does not scale strictly linearly; the full-run limit is set from it with a 2x margin.
FULL="$NCAND"
PROJ_SRCH="$(echo "$SRCH_S * $FULL / $SAMPLE_N" | bc -l)"
PROJ_CLU_2X="$(echo "$CLU_S * 2 * $FULL / $SAMPLE_N" | bc -l)"
{
  echo "sample_candidates=$SAMPLE_N full_candidates=$FULL"
  printf 'projected_search_seconds=%.0f\n' "$PROJ_SRCH"
  printf 'projected_cluster_seconds_2x_margin=%.0f\n' "$PROJ_CLU_2X"
  echo "suggested_pairs_time_limit_h=$(echo "scale=2; ($PROJ_SRCH * 1.5) / 3600" | bc -l)"
  echo "suggested_cluster_time_limit_h=$(echo "scale=2; ($PROJ_CLU_2X * 1.5) / 3600" | bc -l)"
  echo "note=linear projection; cluster scaling is not linear; verify with the full run"
  echo "mmseqs=$("$MMSEQS_BIN" version)"
  echo "date=$(date +%F)"
} > "$TMP/pilot_projection.txt"

cp "$TMP/pilot.tsv" "$OUTDIR/.tmp.pilot.tsv" && mv "$OUTDIR/.tmp.pilot.tsv" "$OUTDIR/pilot.tsv"
cp "$TMP/pilot_projection.txt" "$OUTDIR/.tmp.pilot_projection.txt" && mv "$OUTDIR/.tmp.pilot_projection.txt" "$OUTDIR/pilot_projection.txt"
cp "$TMP/pilot_projection.txt" "$OUTDIR/.tmp.pilot.done" && mv "$OUTDIR/.tmp.pilot.done" "$OUTDIR/pilot.done"
