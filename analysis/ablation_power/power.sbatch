#!/bin/bash -l
# Ablation power calibration: how large an improvement could adhesion_v1 actually detect?
#
# Runs the simulation at both declared training pools. CPU only, no GPU, no network.
#
# Why this is a job rather than a login-node command: the serial version took over 15 minutes for
# ONE training pool at 200 replicates, and 200 replicates is not enough. A power estimate near 0.8
# from 200 draws has a standard error of about 2.8%, and the minimum detectable effect inherits
# that noise. Every (shape, effect, replicate) is independent, so the work parallelises almost
# perfectly and the replicate count can go up instead.
#
# NOT `python3 -I`: isolated mode implies -s, which drops user site-packages, and numpy, scipy and
# scikit-learn live in ~/.local/lib/python3.12 on this cluster. The run records their versions.
#
# Inputs:  $PROJ_ROOT (a clone holding this commit), $BENCH (the frozen benchmark), $OUT (must NOT exist).
# Outputs: $OUT/{all,scored}/{baseline.json,type_i.tsv,power.tsv,summary.md}, provenance.txt
#
# Submit (after the go-ahead; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,BENCH=/path/benchmark,OUT=/path/new_out \
#     -o /path/logs/%x.%j.log -e /path/logs/%x.%j.log analysis/ablation_power/power.sbatch
#SBATCH -p exfab
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J ablation_power
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (a clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${BENCH:?export BENCH (the frozen benchmark directory) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
PINNED="analysis/ablation_power"
PY="${PY:-/usr/bin/python3.12}"
T="${SLURM_CPUS_PER_TASK:-16}"
REPLICATES="${REPLICATES:-1000}"
PERMUTATIONS="${PERMUTATIONS:-1000}"

fail() { echo "FATAL: $*" >&2; exit 1; }

git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$PINNED" || fail "$PINNED changed since $CODE_COMMIT"
[ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$PINNED")" ] || fail "uncommitted changes in $PINNED"
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -d "$BENCH" ] || fail "$BENCH is not a directory"
(cd "$BENCH" && sha256sum -c --quiet SHA256SUMS) || fail "$BENCH does not match its SHA256SUMS"
"$PY" -c "import numpy, scipy, sklearn" || fail "$PY cannot import numpy/scipy/sklearn (do NOT add -I)"

cd "$PROJ_ROOT"
mkdir -p "$OUT"

# Both declared training pools. The choice was an undeclared assumption in the original plan and it
# changes the leakage drop counts, so neither is treated as the default truth.
for POOL in all scored; do
  echo "=== training pool: $POOL ==="
  "$PY" analysis/ablation_power/run_power.py \
    --benchmark "$BENCH" \
    --out "$OUT/$POOL" \
    --train-pool "$POOL" \
    --replicates "$REPLICATES" \
    --permutations "$PERMUTATIONS" \
    --workers "$T" \
    --seed 0
done

{
  echo "job: ablation_power (power calibration for the 2b harness)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "node: $(hostname); cores: $T"
  echo "python: $("$PY" --version 2>&1)"
  echo "packages: $("$PY" -c 'import numpy,scipy,sklearn;print("numpy",numpy.__version__,"scipy",scipy.__version__,"sklearn",sklearn.__version__)')"
  echo "benchmark: $BENCH sha256(SHA256SUMS) $(sha256sum < "$BENCH/SHA256SUMS" | cut -d' ' -f1)"
  echo "replicates: $REPLICATES; permutations: $PERMUTATIONS; seed 0; alpha 0.05"
  echo "training pools: all, scored"
} > "$OUT/provenance.txt"

echo "wrote $OUT"
