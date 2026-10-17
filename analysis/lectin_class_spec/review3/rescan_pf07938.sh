#!/usr/bin/env bash
# Review 3: independent re-run of hmmsearch --cut_ga with PF07938 (Pfam 38.2, HMMER 3.4)
# on a random sample of the 653 discovery proteomes (seeded shuf), login node, 2 CPUs.
# Usage: rescan_pf07938.sh WORK_DIR N_SAMPLE
# WORK_DIR must hold pf07938.hmm (hmmfetch Pfam-A.hmm Fungal_lectin). Output: WORK_DIR/sample.txt, WORK_DIR/tbl/*.tbl
set -euo pipefail
D=${1:?work dir}
N=${2:-25}
SEL=/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/fungi_scan/selection_discovery.tsv
H=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch
mkdir -p "$D/tbl"
tail -n +2 "$SEL" | cut -f2,3 | shuf --random-source=<(yes 20261010) -n "$N" > "$D/sample.txt"
while IFS=$'\t' read -r name fa; do
  "$H" --cut_ga --cpu 2 --noali --tblout "$D/tbl/$name.tbl" "$D/pf07938.hmm" "$fa" > /dev/null
done < "$D/sample.txt"
echo DONE
