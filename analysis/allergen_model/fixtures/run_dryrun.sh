#!/bin/bash -l
# Dry run of the MMseqs2 job scripts on the T1.6 fixture, outside SLURM.
# The scripts keep their #SBATCH headers, which bash ignores. SCRATCH is set to a local directory.
# Usage: bash run_dryrun.sh <workdir>   (default: fixtures/dryrun_out, created here)
set -euo pipefail
FIX="/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/fixtures"
JOBS="/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/jobs"
WORK="${1:-$FIX/dryrun_out}"
rm -rf "$WORK"
mkdir -p "$WORK/scratch" "$WORK/out/cluster" "$WORK/out/pairs" "$WORK/out/xsplit"
export SCRATCH="$WORK/scratch"

# Clustering input: positives and candidates together. The job checks the P_ and N_ prefixes.
cat "$FIX/data/positives.faa" "$FIX/data/candidates.faa" > "$WORK/cluster_in.faa"
SEQS="$WORK/cluster_in.faa" OUTDIR="$WORK/out/cluster" ARM=shared bash -l "$JOBS/cluster_mmseqs.sbatch"

# Pairs: candidates against positives.
QUERY="$FIX/data/candidates.faa" TARGET="$FIX/data/positives.faa" OUTDIR="$WORK/out/pairs" ARM=shared \
  bash -l "$JOBS/pairs_search.sbatch"

# Cross-split check on the split FASTA.
SPLIT_FASTA="$FIX/data/split.faa" OUTDIR="$WORK/out/xsplit" ARM=shared bash -l "$JOBS/cross_split_check.sbatch"

echo "dry run finished: $WORK"
