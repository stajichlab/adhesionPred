#!/bin/bash -l
# Align tips.faa (mafft), trim (trimAl -automated1), tree (IQ-TREE 2, ModelFinder, 1000 UFBoot).
# Environment: WORK (has tips.faa). Outputs in $WORK/tree/. Temp files in node-local $SCRATCH.
# Submit: sbatch -p highclock -c 8 --mem=16G -t 2:00:00 --export=ALL,WORK=... 04_tree.sh
#SBATCH -p highclock
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -J pf28404_tree
set -euo pipefail
: "${WORK:?export WORK}"
module load mafft/7.505 trimal/1.4.1 iqtree/2.2.2.6
OUT="$WORK/tree"; mkdir -p "$OUT"
T="${SLURM_CPUS_PER_TASK:-4}"
mafft --auto --thread "$T" "$WORK/tips.faa" > "$OUT/tips.aln.faa" 2> "$OUT/mafft.log"
trimal -in "$OUT/tips.aln.faa" -out "$OUT/tips.trim.faa" -automated1 -htmlout "$OUT/trimal.html" > "$OUT/trimal.log" 2>&1
python3 - "$OUT/tips.aln.faa" "$OUT/tips.trim.faa" <<'PY' | tee "$OUT/columns.txt"
import sys
def width(p):
    seq, first = 0, True
    for ln in open(p):
        if ln.startswith(">"):
            if not first:
                break
            first = False
        else:
            seq += len(ln.strip())
    return seq
print(f"alignment columns: {width(sys.argv[1])} -> trimmed: {width(sys.argv[2])}")
PY
iqtree2 -s "$OUT/tips.trim.faa" -m MFP -B 1000 -T "$T" --seed 20261002 --prefix "$OUT/pf28404" -redo > "$OUT/iqtree.stdout" 2>&1
echo done
