#!/usr/bin/env bash
# Reproduce the class 2b structural analysis end to end.
#
#   bash analysis/class2b_structure/run.sh            # submit the Foldseek jobs
#
# Must be run from the repository root. Foldseek is built for AVX2 and dies with
# "Illegal instruction" on the older UCR HPCC nodes (c01 has no AVX2). Steps 3 and 4
# therefore go to the `epyc` partition through sbatch; steps 1, 2 and 5 are plain CPU/network
# work and run wherever you are.
set -euo pipefail

REPO="${REPO:-$PWD}"
DIR="$REPO/analysis/class2b_structure"
WORK="$REPO/_workdir/class2b_structure"
DB="$REPO/_workdir/foldseek_db"
FOLDSEEK="${FOLDSEEK:-$REPO/_workdir/bin/foldseek/bin/foldseek}"

if [ ! -x "$FOLDSEEK" ]; then
    echo "Foldseek is not installed. Install the static binary:"
    echo "  mkdir -p $REPO/_workdir/bin && cd $REPO/_workdir/bin"
    echo "  curl -sL https://github.com/steineggerlab/foldseek/releases/download/9-427df8a/foldseek-linux-avx2.tar.gz | tar xz"
    exit 1
fi

# 1-2: fetch AlphaFold DB models and extract the pLDDT>=70 confident cores.
/usr/bin/python3.12 "$DIR/01_fetch_structures.py" --outdir "$WORK"
/usr/bin/python3.12 "$DIR/02_confidence_profile.py" --workdir "$WORK"

# 3-4: Foldseek. Needs an AVX2 node.
mkdir -p "$WORK/logs"
JOB=$(sbatch --parsable -p epyc -c 8 --mem 16G -t 1:30:00 -J class2b_foldseek \
      -o "$WORK/logs/foldseek_%j.out" -e "$WORK/logs/foldseek_%j.err" \
      --wrap "bash '$DIR/03_foldseek_allvall.sh' '$WORK' '$FOLDSEEK' && \
              bash '$DIR/04_foldseek_pdb.sh' '$WORK' '$FOLDSEEK' '$DB'")
echo "submitted Foldseek job $JOB; logs in $WORK/logs/"
echo "  bash '$DIR/06_targeted_pairwise.sh' (submit to epyc too; run.sh does 03-04 only)"
echo "when it finishes:"
echo "  /usr/bin/python3.12 $DIR/05_summarize.py --workdir $WORK"
