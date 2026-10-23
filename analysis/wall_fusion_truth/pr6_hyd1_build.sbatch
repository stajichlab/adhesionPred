#!/bin/bash -l
# PR 6 for hyd1_like (DRAFT, NOT SUBMITTED; nothing of this has run on HPCC). Modeled on pr6_ecm33_build.sbatch, trimmed. SLURM array: task i is the i-th
# (unit, mode) pair of `pr6_hyd1_build.py tasks` (clade units and the Kiss unit in both independence modes, cluster units in the headline mode, then `all`).
# Per task: manifest -> alignment -> HMM -> false-hit searches -> cutoff -> ONLY THEN the held-out searches (refused if cutoff.json is missing).
# It writes no call, no status, no tier and does not touch data/sorting_hat. A custom model is report-only (decision 6).
#
# THE FREEZE COMES FIRST: $PROJ_ROOT/analysis/wall_fusion_truth/pr6_hyd1_freeze.json must say "frozen": true. Only the owner turns
# pr6_hyd1_freeze_PROPOSED.json into it (the freeze commit). `tasks` and `manifest` refuse an unfrozen file, so a missing freeze stops the array at task 0.
# Needs the owner's go-ahead for EACH step (D21): clone update, pre-check job, extract job, freeze commit, this job.
#
# Inputs: $EXTRACT (output of pr6_hyd1_extract.sbatch), $PROTEOME_DIR/<name>.faa (the four kept proteomes), Pfam-A (PFAM_HMM, release PFAM_RELEASE),
#   data/controls/hard-negatives-shared (161 shared negatives), data/sorting_hat/hydrophobin_relaxed.hmm (committed relaxed panel comparator).
# Modules: hmmer/3.4 mafft/7.505 seqkit/2.4.0 (MMseqs2 ran in the extract job). UNCONFIRMED: `esl-reformat` and jackhmmer taking the seed as written;
# run time of L-INS-i on <= 30 short sequences (minutes expected). Resources are a guess.
#
# Submit (after the go-ahead). N = number of tasks - 1:
#   N=$(( $(/usr/bin/python3.12 -I $PROJ_ROOT/analysis/wall_fusion_truth/pr6_hyd1_build.py tasks --units $EXTRACT/members/units.tsv \
#        --units-cluster $EXTRACT/units_cluster.tsv --freeze $PROJ_ROOT/analysis/wall_fusion_truth/pr6_hyd1_freeze.json | wc -l) - 1 ))
#   sbatch --array=0-$N%4 --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,EXTRACT=/path/to/extract_out,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o /path/to/logs/%x.%A_%a.log -e /path/to/logs/%x.%A_%a.log analysis/wall_fusion_truth/pr6_hyd1_build.sbatch
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH -J pr6_hyd1_build
set -euo pipefail
: "${PROJ_ROOT:?}" "${CODE_COMMIT:?}" "${EXTRACT:?}" "${OUT:?}" "${PFAM_RELEASE:?}" "${SLURM_ARRAY_TASK_ID:?this is an array job}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
DATA_ROOT="${DATA_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
PROTEOME_DIR="${PROTEOME_DIR:-$DATA_ROOT/_workdir/sorting_hat}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_hyd1_build.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-8}"
WF=analysis/wall_fusion_truth
NEG=data/controls/hard-negatives-shared
FREEZE="$PROJ_ROOT/$WF/pr6_hyd1_freeze.json"
PROTEOMES="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"
PANEL_NAMES="Hydrophobin Hydrophobin_2 Eas DewD Hyd1F Hydrophobin_D Hydrophobin_like HsbA"
QUERY=A0ACE4KFX7
PY=/usr/bin/python3.12
BUILD="$PY -I $PROJ_ROOT/$WF/pr6_hyd1_build.py"
fail() { echo "FATAL: $*" >&2; exit 1; }
jget() { $PY -I -c "import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])" "$FREEZE" "$1"; }

git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$NEG" data/sorting_hat/hydrophobin_relaxed.hmm scripts/sorting_hat/refresh_3_10/inputs.tsv \
  || fail "reviewed inputs changed since $CODE_COMMIT"
[ -f "$FREEZE" ] || fail "$FREEZE is missing: the freeze commit comes before any build"
grep -q '"frozen": true' "$FREEZE" || fail "$FREEZE is not frozen"
[ "$(sha256sum < "$EXTRACT/hyd1_members.faa" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$EXTRACT/hyd1_members.faa.sha256")" ] || fail "hyd1_members.faa differs from its .sha256"
FZ_READING="$($PY -I -c "import json,sys; print(json.load(open(sys.argv[1]))['training_rule']['reading'])" "$FREEZE")"
grep -q "reading $FZ_READING\$" "$EXTRACT/provenance.txt" || fail "the extract job used another training-rule reading than the freeze ($FZ_READING): rerun it with READING=$FZ_READING"
MEM="$EXTRACT/members/members.tsv"; MFA="$EXTRACT/hyd1_members.faa"; CLU="$EXTRACT/clu_cluster.tsv"

mkdir -p "$TMP/w" "$OUT"; W="$TMP/w"
$BUILD tasks --units "$EXTRACT/members/units.tsv" --units-cluster "$EXTRACT/units_cluster.tsv" --freeze "$FREEZE" > "$W/tasks.list"
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$W/tasks.list")"
[ -n "$LINE" ] || { echo "task $SLURM_ARRAY_TASK_ID has no (unit, mode) pair; nothing to do"; exit 0; }
UNIT="$(printf '%s' "$LINE" | cut -f1)"; MODE="$(printf '%s' "$LINE" | cut -f2)"
UDIR="$OUT/${UNIT//:/_}/$MODE"
[ ! -e "$UDIR" ] || fail "$UDIR exists; nothing is overwritten"
IS_KISS=no; if [ "$UNIT" != all ] && awk -F'\t' -v u="$UNIT" '$1==u && $NF=="kiss_single" {f=1} END{exit !f}' "$EXTRACT/members/units.tsv" "$EXTRACT/units_cluster.tsv"; then IS_KISS=yes; fi
echo "task $SLURM_ARRAY_TASK_ID: unit $UNIT, mode $MODE, kiss unit: $IS_KISS"

for p in $PROTEOMES; do
  f="$PROTEOME_DIR/$p.faa"; [ -s "$f" ] || fail "$f is missing"
  want=$(awk -F'\t' -v p="$p" '$1==p {print $3}' "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/inputs.tsv")
  [ "$want" = "$(grep -c '^>' "$f")" ] || fail "$f does not have the $want proteins of inputs.tsv"
  printf '%s\t%s\n' "$p" "$want" >> "$W/sizes.body"
done
{ printf 'proteome\tn_proteins\n'; cat "$W/sizes.body"; } > "$W/sizes.tsv"
module load hmmer/3.4 mafft/7.505 seqkit/2.4.0
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac

# 1. Manifest (refuses an unfrozen freeze). Exit 3 = not buildable: recorded, no model. A cluster unit has one column (cluster30 = none).
set +e
$BUILD manifest --members "$MEM" --unit "$UNIT" --clusters "$CLU" --mode "$MODE" --freeze "$FREEZE" --out "$W/unit"
RC=$?
set -e
if [ "$RC" -eq 3 ]; then mkdir -p "$UDIR"; cp "$W/unit/manifest.json" "$UDIR/"; echo "$UNIT ($MODE): not_buildable; no model"; exit 0
elif [ "$RC" -ne 0 ]; then fail "manifest failed (exit $RC)"; fi

# 2. Alignment (L-INS-i), the manifest check (no held-out id, no excluded neighbour) and the freeze's occupancy check; a failure stops the unit.
seqkit grep -f "$W/unit/train_ids.txt" "$MFA" > "$W/train.faa"
[ "$(grep -c '^>' "$W/train.faa")" -eq "$(wc -l < "$W/unit/train_ids.txt")" ] || fail "train.faa does not hold every training id"
mafft --localpair --maxiterate 1000 --quiet --thread "$THREADS" "$W/train.faa" > "$W/train.aln.fa"
$BUILD check-alignment --manifest "$W/unit/manifest.json" --alignment "$W/train.aln.fa" --freeze "$FREEZE" || fail "alignment check failed for $UNIT"
hmmbuild --amino --informat afa --cpu "$THREADS" -n "hyd1_like_${UNIT//:/_}" "$W/model.hmm" "$W/train.aln.fa" > "$W/hmmbuild.log"

# 3. The false-hit side: the four kept proteomes and the shared negatives. `labelled` = phmmer of the Kiss protein at E 1e-5 (the freeze's labelled_rule),
#    plus, when the freeze says proteome_other_hydrophobin = ignore, the proteome proteins hit by the hydrophobin-class panel at GA.
FLOOR="$(jget cutoff_floor)"; Z="$(jget evalue_db_size)"; OTHER="$(jget proteome_other_hydrophobin)"
awk -F'\t' 'NR>1 {print $1 "\t" $6}' "$PROJ_ROOT/$NEG/controls.tsv" > "$W/ng.body"; { printf 'accession\tgroup\n'; cat "$W/ng.body"; } > "$W/neg_groups.tsv"
tr ' ' '\n' <<< "$PANEL_NAMES" > "$W/panel.names"; hmmfetch -f "$RESOLVED" "$W/panel.names" > "$W/panel.hmm"
[ "$(grep -c '^ACC ' "$W/panel.hmm")" -eq 8 ] || fail "panel.hmm does not hold 8 models"
awk -v q="$QUERY" '/^>/ { keep = (substr($1, 2) == q) } keep' "$MFA" > "$W/query.faa"; [ -s "$W/query.faa" ] || fail "$QUERY is not in hyd1_members.faa"
: > "$W/labelled.txt"; printf 'proteome\tid\tscore\n' > "$W/proteome_hits.tsv"
for p in $PROTEOMES; do
  hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/c_$p.tbl" "$W/model.hmm" "$PROTEOME_DIR/$p.faa"
  awk -v p="$p" 'BEGIN{OFS="\t"} !/^#/ {print p, $1, $6}' "$W/c_$p.tbl" >> "$W/proteome_hits.tsv"
  phmmer -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$W/l_$p.tbl" "$W/query.faa" "$PROTEOME_DIR/$p.faa"
  awk '!/^#/ {print $1}' "$W/l_$p.tbl" >> "$W/labelled.txt"
  if [ "$OTHER" = ignore ]; then
    hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$W/g_$p.tbl" "$W/panel.hmm" "$PROTEOME_DIR/$p.faa"; awk '!/^#/ {print $1}' "$W/g_$p.tbl" >> "$W/labelled.txt"
  fi
done
sort -u "$W/labelled.txt" -o "$W/labelled.txt"
hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/c_neg.tbl" "$W/model.hmm" "$PROJ_ROOT/$NEG/sequences.faa"
{ printf 'accession\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/c_neg.tbl"; } > "$W/neg_hits.tsv"

# 4. The cutoff: no positive score is an input; refuses an unfrozen freeze.
$BUILD cutoff --freeze "$FREEZE" --proteome-hits "$W/proteome_hits.tsv" --labelled "$W/labelled.txt" --neg-groups "$W/neg_groups.tsv" \
  --neg-hits "$W/neg_hits.tsv" --sizes "$W/sizes.tsv" --out "$W/cutoff.json"

# 5. ONLY NOW the held-out members. Comparators independent of the training alignment: the Pfam panel at GA, the committed relaxed panel, phmmer of the
#    single Kiss protein; the ceiling: jackhmmer -N 3 from the single Kiss protein. -Z fixes the E-value database size to UniRef50 (the freeze), because the
#    held-out file is tiny and an unset -Z would make every E-value look better than in the search that found the members. In the Kiss unit the held-out member
#    IS the seed, so the two seed comparators are not run (the report says seed_in_unit).
[ -s "$W/cutoff.json" ] || fail "no cutoff: the held-out step does not start"
if [ "$UNIT" != all ]; then   # the `all` model has no held-out unit
seqkit grep -f "$W/unit/heldout_ids.txt" "$MFA" > "$W/heldout.faa"
hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_custom.tbl" "$W/model.hmm" "$W/heldout.faa"
hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_ga.tbl" "$W/panel.hmm" "$W/heldout.faa"
hmmsearch --nobias --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_rel.tbl" "$PROJ_ROOT/data/sorting_hat/hydrophobin_relaxed.hmm" "$W/heldout.faa"
{ printf 'id\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/h_custom.tbl"; } > "$W/hits_custom.tsv"
{ printf 'id\n'; awk '!/^#/ {print $1}' "$W/h_ga.tbl" | sort -u; } > "$W/hits_pfam_ga.tsv"
{ printf 'id\n'; awk '!/^#/ {print $1}' "$W/h_rel.tbl" | sort -u; } > "$W/hits_relaxed_panel.tsv"
if [ "$IS_KISS" = no ]; then
  phmmer -Z "$Z" -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_seed.tbl" "$W/query.faa" "$W/heldout.faa"
  jackhmmer -Z "$Z" -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_jack.tbl" "$W/query.faa" "$W/heldout.faa"
  { printf 'id\tevalue\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $5}' "$W/h_seed.tbl"; } > "$W/hits_phmmer_seed.tsv"
  { printf 'id\tevalue\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $5}' "$W/h_jack.tbl"; } > "$W/hits_jackhmmer_seed.tsv"
fi
fi
shopt -s nullglob

# 6. Provenance and one atomic move into $OUT/<unit>/<mode>.
{
  echo "job: pr6_hyd1_build (unit $UNIT, mode $MODE); commit $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned $CODE_COMMIT)"
  echo "freeze sha256: $(sha256sum < "$FREEZE" | cut -d' ' -f1)"
  echo "hmmer: $(hmmsearch -h | sed -n 2p); mafft: $(mafft --version 2>&1 | head -1); pfam_resolved: $RESOLVED"
  echo "hyd1_members.faa sha256: $(sha256sum < "$MFA" | cut -d' ' -f1); model.hmm sha256: $(sha256sum < "$W/model.hmm" | cut -d' ' -f1)"
} > "$W/provenance.txt"
mkdir -p "$(dirname "$UDIR")" "$UDIR.tmp.$$"
cp "$W/unit/manifest.json" "$W/unit/train_ids.txt" "$W/unit/heldout_ids.txt" "$W/train.aln.fa" "$W/model.hmm" "$W/hmmbuild.log" "$W/cutoff.json" \
   "$W/neg_hits.tsv" "$W/proteome_hits.tsv" "$W/provenance.txt" "$W"/hits_*.tsv "$UDIR.tmp.$$/"
mv "$UDIR.tmp.$$" "$UDIR"
echo "wrote $UDIR"
