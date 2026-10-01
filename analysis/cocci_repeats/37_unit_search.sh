#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 16 --mem 48gb --time 2:00:00
#SBATCH -J unitsrch -o logs/37_unit_search.%A.log
# Is the SOWgp repeat unit found in any other protein? Three searches over 6,313 proteomes
# (5,813 Fungi_5k + 7 long-read Coccidioides + 493 pangenome), 16 shards in parallel in one job.
#
#   A. unit HMM      hmmsearch of sowgp_unit.hmm (34) on every proteome. All domain hits
#                    with E <= 1000 at a fixed database size -Z 6e7 (~ Fungi_5k protein
#                    count), so scores and E-values are comparable across proteomes.
#   B. anchor        30_anchor_family_search.py search, motifs KKYGDC and PTDCYGDC, every hit
#                    aligned to SOWgp58. Membership is decided by the alignment, not the motif.
#   C. architecture  14_repeat_detect_general.py on the Onygenales proteomes (Fungi_5k) only:
#                    periodicity + composition, for the "Pro/Cys-rich repeat, no unit
#                    homology" level.
#
# Needs sowgp_unit.hmm (34_sowgp_unit_hmm.py) and unit_genomes.tsv (35_unit_genomes.py).
# Submit from analysis/cocci_repeats:   sbatch 37_unit_search.sh
#
# NOTE: no $(dirname "${BASH_SOURCE[0]}") -- wrong under sbatch. Use $SLURM_SUBMIT_DIR.

set -euo pipefail
DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$DIR"
mkdir -p logs
module load hmmer/3.4

WORK="${SCRATCH:?SCRATCH is not set}/unitsrch.$$"
mkdir -p "$WORK"
NS=16
REF=sowgp_seed.fa
REFID=SOWgp58_Cocci_immitis_published
[ -f "$REF" ] || gunzip -kc "$REF.gz" > "$REF"

# column 2 = path, column 8 = order, column 1 = source (header row skipped)
tail -n +2 unit_genomes.tsv | cut -f2 > "$WORK/all.txt"
tail -n +2 unit_genomes.tsv | awk -F'\t' '($1=="fungi5k" && $8=="Onygenales") || $1=="cocci_longread" {print $2}' > "$WORK/ony.txt"
wc -l "$WORK/all.txt" "$WORK/ony.txt"
split -n r/$NS -d "$WORK/all.txt" "$WORK/A."
split -n r/$NS -d "$WORK/ony.txt" "$WORK/C."

# ---- A + B per shard
PIDS=()
for S in "$WORK"/A.??; do
  (
    : > "$S.dom"
    while read -r F; do
      L=$(basename "$F"); L=${L%.proteins.fa}; L=${L%.fasta}; L=${L%.fa}
      hmmsearch --cpu 1 -Z 60000000 --domZ 60000000 -E 1000 --domE 1000 --noali \
          --domtblout "$S.one" -o /dev/null sowgp_unit.hmm "$F"
      grep -v '^#' "$S.one" | awk -v L="$L" '{print L"\t"$1"\t"$3"\t"$6"\t"$10"\t"$11"\t"$13"\t"$14"\t"$16"\t"$17"\t"$18"\t"$19}' >> "$S.dom" || true
    done < "$S"
    /usr/bin/python3.12 30_anchor_family_search.py search \
        --motif KKYGDC PTDCYGDC --reference "$REF" --reference-id "$REFID" --period 47 \
        --no-gff --proteomes $(cat "$S") --out "$S.anchor.tsv" > "$S.anchor.log" 2>&1
  ) &
  PIDS+=($!)
done
for P in "${PIDS[@]}"; do wait "$P"; done
echo "phase A/B done"

# ---- C per shard (Onygenales architecture scan)
PIDS=()
for S in "$WORK"/C.??; do
  /usr/bin/python3.12 14_repeat_detect_general.py $(cat "$S") --out "$S.arch.tsv" \
      > "$S.arch.log" 2>&1 &
  PIDS+=($!)
done
for P in "${PIDS[@]}"; do wait "$P"; done
echo "phase C done"

# ---- merge and compress
{ printf 'proteome\tprotein\tlength\tfull_score\tdom_n\tdom_of\tdom_ievalue\tdom_score\thmm_from\thmm_to\tali_from\tali_to\n'
  cat "$WORK"/A.??.dom; } | gzip > unit_hmm_hits.tsv.gz
head -1 "$WORK/A.00.anchor.tsv" > "$WORK/anchor.all"
for S in "$WORK"/A.??.anchor.tsv; do tail -n +2 "$S" >> "$WORK/anchor.all"; done
gzip -c "$WORK/anchor.all" > unit_anchor_hits.tsv.gz
head -1 "$WORK/C.00.arch.tsv" > "$WORK/arch.all"
for S in "$WORK"/C.??.arch.tsv; do tail -n +2 "$S" >> "$WORK/arch.all"; done
gzip -c "$WORK/arch.all" > unit_architecture_onygenales.tsv.gz
tail -n 2 "$WORK"/A.??.anchor.log | grep -E "proteomes with|calls" || true
zcat unit_hmm_hits.tsv.gz | wc -l; zcat unit_anchor_hits.tsv.gz | wc -l; zcat unit_architecture_onygenales.tsv.gz | wc -l
rm -rf "$WORK"
