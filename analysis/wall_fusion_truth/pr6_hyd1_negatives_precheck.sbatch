#!/bin/bash -l
# PR 6 hyd1_like pre-freeze check (DRAFT, NOT SUBMITTED; nothing of this has run on HPCC). Modeled on pr6_ecm33_negatives_precheck.sbatch.
# Which of the 161 shared hard negatives are hyd1_like family members (phmmer of A0ACE4KFX7, E <= 1e-5: the freeze's labelled_rule) or
# hydrophobin-class proteins (the eight hydrophobin-class Pfam models at GA)? The same two searches over the four kept proteomes give the proteome
# proteins that the freeze would call `labelled` (phmmer) and the `other hydrophobin` proteins (panel). A third, descriptive search says how many
# curated hydrophobins (the committed truth proteins and the 8 other Kiss hydrophobin positives) the Kiss protein is related to at all.
# It builds no model, scores no held-out unit, writes no call, status or tier, and does not read or change any freeze file.
#
# Needs the owner's go-ahead BEFORE it is submitted (D21): the clone update and this job. A small job (minutes), still a SLURM job.
# Inputs: $IN/kiss_sequences.faa (+ .sha256; already staged, A0ACE4KFX7 is taken from it), $PROJ_ROOT (a clone at CODE_COMMIT: shared negatives,
#   truth proteins, scripts), $PROTEOME_DIR (the four kept proteomes; protein counts are checked against inputs.tsv), Pfam-A (PFAM_HMM, release PFAM_RELEASE).
# Modules: hmmer/3.4. NOT RUN: the exact options below are unconfirmed on HPCC.
#
# Submit (after the go-ahead):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_hyd1_negatives_precheck.sbatch
#   (the log directory must exist; $OUT must NOT exist.)
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH -J pr6_hyd1_negcheck
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
QUERY=A0ACE4KFX7
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
DATA_ROOT="${DATA_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
PROTEOME_DIR="${PROTEOME_DIR:-$DATA_ROOT/_workdir/sorting_hat}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_hyd1_negcheck.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
WF=analysis/wall_fusion_truth
NEG=data/controls/hard-negatives-shared
TRUTH=analysis/hydrophobin_truth/positives_t2t3.faa
PROTEOMES="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"
# the model NAMES of the eight hydrophobin-class models of the panel (data/sorting_hat/family_table.tsv: PF01185 PF06766 PF22354 PF28987 PF29802 PF29785 PF29465 PF12296)
PANEL_NAMES="Hydrophobin Hydrophobin_2 Eas DewD Hyd1F Hydrophobin_D Hydrophobin_like HsbA"
# the other Kiss hydrophobin positives (controls.tsv annotation 'hydrophobin'; the descriptive relatedness set)
KISS_HYD="A0A0D1BZL3 Q7S3P5 P79073 I1RDW3 A0ACE4JU10 P41746 E9QT94 Q04571"
fail() { echo "FATAL: $*" >&2; exit 1; }
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$NEG" analysis/hydrophobin_truth/positives_t2t3.faa scripts/sorting_hat/refresh_3_10/inputs.tsv \
    || fail "inputs changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$NEG")" ] || fail "uncommitted changes in $WF or $NEG"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }
pick() { awk -v ids="$2" 'BEGIN { n = split(ids, a, " "); for (i = 1; i <= n; i++) want[a[i]] = 1 } /^>/ { keep = (substr($1, 2) in want) } keep' "$1"; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -f "$IN/kiss_sequences.faa" ] && sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa is missing or differs from its .sha256"
for f in sequences.faa controls.tsv; do [ -s "$PROJ_ROOT/$NEG/$f" ] || fail "$NEG/$f is missing or empty"; done
[ -s "$PROJ_ROOT/$TRUTH" ] || fail "$TRUTH is missing"
mkdir -p "$TMP"
grep '^>' "$PROJ_ROOT/$NEG/sequences.faa" | sed 's/^>//; s/[[:space:]].*//' | sort > "$TMP/fa_ids"
tail -n +2 "$PROJ_ROOT/$NEG/controls.tsv" | cut -f1 | sort > "$TMP/ctl_ids"
cmp -s "$TMP/fa_ids" "$TMP/ctl_ids" || fail "sequence names of $NEG/sequences.faa differ from the accessions of controls.tsv"
[ "$(wc -l < "$TMP/fa_ids")" -eq 161 ] || fail "$NEG does not hold 161 sequences"
pick "$IN/kiss_sequences.faa" "$QUERY" > "$TMP/query.faa"
[ "$(grep -c '^>' "$TMP/query.faa")" -eq 1 ] || fail "$QUERY is not in kiss_sequences.faa exactly once"
pick "$IN/kiss_sequences.faa" "$KISS_HYD" > "$TMP/kiss_hyd.faa"
[ "$(grep -c '^>' "$TMP/kiss_hyd.faa")" -eq 8 ] || fail "the 8 other Kiss hydrophobin positives are not all in kiss_sequences.faa"
cat "$PROJ_ROOT/$TRUTH" "$TMP/kiss_hyd.faa" > "$TMP/related_targets.faa"
for p in $PROTEOMES; do
  f="$PROTEOME_DIR/$p.faa"; [ -s "$f" ] || fail "$f is missing"
  want=$(awk -F'\t' -v p="$p" '$1==p {print $3}' "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/inputs.tsv")
  [ "$want" = "$(grep -c '^>' "$f")" ] || fail "$f does not have the $want proteins of inputs.tsv"
done

module load hmmer/3.4
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
tr ' ' '\n' <<< "$PANEL_NAMES" > "$TMP/panel.names"
hmmfetch -f "$RESOLVED" "$TMP/panel.names" > "$TMP/panel.hmm"
[ "$(grep -c '^ACC ' "$TMP/panel.hmm")" -eq 8 ] || fail "panel.hmm does not hold 8 models"

ok_tbl() { tail -n 3 "$1" | grep -qx '# \[ok\]' || fail "$1 has no '# [ok]' trailer"; }
NEGF="$PROJ_ROOT/$NEG/sequences.faa"
phmmer -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/phmmer_E1e-5.tbl" "$TMP/query.faa" "$NEGF"; ok_tbl "$TMP/phmmer_E1e-5.tbl"
phmmer -E 1e-3 --domE 1e-3 --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/phmmer_E1e-3.tbl" "$TMP/query.faa" "$NEGF"; ok_tbl "$TMP/phmmer_E1e-3.tbl"
hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/panel_ga.tbl" "$TMP/panel.hmm" "$NEGF"; ok_tbl "$TMP/panel_ga.tbl"
phmmer -E 1e-3 --domE 1e-3 --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/related_E1e-3.tbl" "$TMP/query.faa" "$TMP/related_targets.faa"; ok_tbl "$TMP/related_E1e-3.tbl"
printf 'proteome\tid\tmodel_or_query\tscore\n' > "$TMP/proteome_labelled_phmmer.tsv"; cp "$TMP/proteome_labelled_phmmer.tsv" "$TMP/proteome_panel_ga.tsv"
for p in $PROTEOMES; do
  phmmer -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/l_$p.tbl" "$TMP/query.faa" "$PROTEOME_DIR/$p.faa"; ok_tbl "$TMP/l_$p.tbl"
  hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/g_$p.tbl" "$TMP/panel.hmm" "$PROTEOME_DIR/$p.faa"; ok_tbl "$TMP/g_$p.tbl"
  awk -v p="$p" 'BEGIN{OFS="\t"} !/^#/ {print p, $1, $3, $6}' "$TMP/l_$p.tbl" >> "$TMP/proteome_labelled_phmmer.tsv"
  awk -v p="$p" 'BEGIN{OFS="\t"} !/^#/ {print p, $1, $3, $6}' "$TMP/g_$p.tbl" >> "$TMP/proteome_panel_ga.tsv"
done

/usr/bin/python3.12 -I "$PROJ_ROOT/$WF/pr6_hyd1_negatives_precheck.py" summarize \
  --strict "$TMP/phmmer_E1e-5.tbl" --weak "$TMP/phmmer_E1e-3.tbl" --panel "$TMP/panel_ga.tbl" --related "$TMP/related_E1e-3.tbl" \
  --controls "$PROJ_ROOT/$NEG/controls.tsv" --negatives-fasta "$NEGF" --out "$TMP/tables"
{
  echo "job: pr6_hyd1_negcheck (pre-freeze check of the hyd1_like hard negatives and proteomes)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(phmmer -h | sed -n 2p)"
  echo "pfam_resolved: $RESOLVED"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "input negatives sequences.faa sha256 $(sha256sum < "$NEGF" | cut -d' ' -f1)"
  echo "input truth proteins sha256 $(sha256sum < "$PROJ_ROOT/$TRUTH" | cut -d' ' -f1)"
  echo "query: $QUERY; panel models: $PANEL_NAMES"
} > "$TMP/tables/provenance.txt"
cp "$TMP"/phmmer_E1e-5.tbl "$TMP"/phmmer_E1e-3.tbl "$TMP"/panel_ga.tbl "$TMP"/related_E1e-3.tbl "$TMP"/proteome_labelled_phmmer.tsv "$TMP"/proteome_panel_ga.tsv "$TMP/tables/"
mkdir -p "$(dirname "$OUT")"
mv "$TMP/tables" "$OUT"
echo "wrote $OUT"
