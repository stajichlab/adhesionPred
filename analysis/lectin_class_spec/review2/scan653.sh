#!/usr/bin/env bash
# Review 2: whole-proteome (ungated) scan of the 653 discovery proteomes with the five rare lectin models.
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/adhesionPred/.claude/worktrees/agent-a487e4e8175609091/_review2_scratch
H=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch
mkdir -p "$D/scan653"
while read -r fa; do
  b=$(basename "$fa" .proteins.fa)
  if [ ! -s "$D/scan653/$b.tbl" ]; then
    "$H" --cut_ga --cpu 2 --noali --tblout "$D/scan653/$b.tbl.tmp" "$D/m.hmm" "$fa" > /dev/null
    mv "$D/scan653/$b.tbl.tmp" "$D/scan653/$b.tbl"
  fi
done < "$D/fastas.txt"
echo DONE > "$D/scan653.done"
