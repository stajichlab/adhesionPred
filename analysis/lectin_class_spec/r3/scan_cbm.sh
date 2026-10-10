#!/usr/bin/env bash
# Revision 3: whole-proteome scan of the 653 discovery proteomes with two chitin-binding models
# (CBM_14 PF01607, Chitin_bind_1 PF00187), Pfam 38.2, hmmsearch --cut_ga, login node, 2 CPUs.
# Usage: scan_cbm.sh WORK_DIR   (WORK_DIR is an empty scratch directory outside git)
set -euo pipefail
D=${1:?work dir}
SEL=/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/fungi_scan/selection_discovery.tsv
PF=/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm
HMMER=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin
mkdir -p "$D/scan"
if [ ! -s "$D/m.hmm" ]; then
  printf 'PF01607\nPF00187\n' > "$D/acc.txt"
  "$HMMER/hmmfetch" -f "$PF" "$D/acc.txt" > "$D/m.hmm" 2>/dev/null || {
    # hmmfetch needs an index for accession fetch; fall back to key fetch by name
    printf 'CBM_14\nChitin_bind_1\n' > "$D/names.txt"
    "$HMMER/hmmfetch" -f "$PF" "$D/names.txt" > "$D/m.hmm"
  }
fi
tail -n +2 "$SEL" | cut -f3 > "$D/fastas.txt"
while read -r fa; do
  b=$(basename "$fa" .proteins.fa)
  if [ ! -s "$D/scan/$b.tbl" ]; then
    "$HMMER/hmmsearch" --cut_ga --cpu 2 --noali --tblout "$D/scan/$b.tbl.tmp" "$D/m.hmm" "$fa" > /dev/null
    mv "$D/scan/$b.tbl.tmp" "$D/scan/$b.tbl"
  fi
done < "$D/fastas.txt"
echo DONE > "$D/scan.done"
