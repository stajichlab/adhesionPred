#!/bin/bash -l
# PR 6 steps 6.1 and 6.2 for ecm33 (DRAFT, NOT SUBMITTED; nothing of PR 6 beyond the gap search has run on HPCC): per held-out unit, build a
# custom Ecm33 HMM from the OTHER clades, set its cutoff from the false-hit side only, then score the held-out members and the three
# comparators. A SLURM array: task i is the i-th (unit, mode) pair: for every held-out unit (status `unit` of units.tsv) the HEADLINE independence mode (cluster30) and the second, DESCRIPTIVE one (none), then the `all` model (no held-out unit, built once). 2 x units + 1 tasks (23 for the 11 units of the 2026-10-22 members).
# It writes no call, no status and no tier and does not touch data/sorting_hat. A custom model is report-only (decision 6).
# Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md (sections 4.3 to 4.8).
#
# THE FREEZE COMES FIRST. $FREEZE (analysis/wall_fusion_truth/pr6_ecm33_freeze.json in the clone) must exist and say "frozen": true before
# this job starts: the proposal is analysis/wall_fusion_truth/pr6_ecm33_freeze_PROPOSED.json and only the owner turns it into the freeze
# commit. The order inside each task is fixed: manifest -> alignment -> HMM -> false-hit searches -> cutoffs -> ONLY THEN the held-out
# searches (the held-out step refuses to start if cutoff.json is missing). No held-out score ever reaches a cutoff.
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update, the extract job and this job. Heavy work is a
# SLURM job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $PROJ_ROOT (a clone at CODE_COMMIT)  scripts, the freeze file, the committed control sets (shared negatives, Kiss controls),
#                                        scripts/sorting_hat/refresh_3_10/inputs.tsv.
#   $EXTRACT                             the OUTPUT directory of pr6_ecm33_extract.sbatch (members/, ecm33_members.faa and its .sha256).
#   $PROTEOME_DIR/<name>.faa             the four kept proteomes (default the OLD dir of scripts/sorting_hat/refresh_3_10/lib.sh); the
#                                        protein count of each must equal n_proteins of inputs.tsv. NOT RE-CHECKED ON HPCC BY THE AUTHOR.
#   Pfam: PFAM_HMM (default the `current` link) must resolve to a directory named *Pfam$PFAM_RELEASE.
# Modules: hmmer/3.4, mafft/7.505, seqkit/2.4.0, MMseqs2/17-b804f (all listed by `module avail` on HPCC on 2026-10-22, except that the exact
# options below were NOT run). UNCONFIRMED: `esl-reformat` (Easel, bundled with HMMER) to turn the aligned FASTA into Stockholm for the
# jackhmmer comparator; whether `jackhmmer` takes the alignment as its query as written; the run time of L-INS-i on up to 480 sequences
# (it is slow: the request below is a guess).
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script and the freeze):
#   N=$(( 2 * $(awk -F'\t' 'NR>1 && $7=="unit"' $EXTRACT/members/units.tsv | wc -l) ))   # 2 modes per unit, then `all`: the array is 0..N (N+1 tasks)
#   sbatch --array=0-$N%4 --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,EXTRACT=/path/to/extract_out,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o /path/to/logs/%x.%A_%a.log -e /path/to/logs/%x.%A_%a.log analysis/wall_fusion_truth/pr6_ecm33_build.sbatch
#   (the directory of the log must exist; $OUT is created by the first task and each unit writes only to $OUT/<unit>, never overwriting.)
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -J pr6_ecm33_build
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT before sbatch}"
: "${EXTRACT:?export EXTRACT (the output of pr6_ecm33_extract.sbatch) before sbatch}"
: "${OUT:?export OUT before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
: "${SLURM_ARRAY_TASK_ID:?this is an array job (sbatch --array=0-N)}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
DATA_ROOT="${DATA_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
PROTEOME_DIR="${PROTEOME_DIR:-$DATA_ROOT/_workdir/sorting_hat}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_ecm33_build.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-8}"
WF=analysis/wall_fusion_truth
NEG=data/controls/hard-negatives-shared
CTRL=data/controls/wall-fusion-kiss2019
FREEZE="$PROJ_ROOT/$WF/pr6_ecm33_freeze.json"
PROTEOMES="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"   # the frozen kept set; checked against the freeze file below
PY=/usr/bin/python3.12
BUILD="$PY -I $PROJ_ROOT/$WF/pr6_ecm33_build.py"

fail() { echo "FATAL: $*" >&2; exit 1; }

check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$NEG" "$CTRL" scripts/sorting_hat/refresh_3_10/inputs.tsv \
    || fail "$WF, $NEG, $CTRL or inputs.tsv changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$NEG" "$CTRL")" ] || fail "uncommitted changes in $WF, $NEG or $CTRL"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

check_pin
[ -f "$FREEZE" ] || fail "$FREEZE is missing: the freeze commit comes before any build"
grep -q '"frozen": true' "$FREEZE" || fail "$FREEZE is not frozen"
[ -f "$EXTRACT/ecm33_members.faa" ] && [ -f "$EXTRACT/ecm33_members.faa.sha256" ] || fail "$EXTRACT is not an extract output"
[ "$(sha256sum < "$EXTRACT/ecm33_members.faa" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$EXTRACT/ecm33_members.faa.sha256")" ] \
  || fail "$EXTRACT/ecm33_members.faa differs from its .sha256"
MEM="$EXTRACT/members/members.tsv"; UNITS="$EXTRACT/members/units.tsv"

# The tasks: for each held-out unit the headline mode then the second mode (both read from the freeze), then `all` (built once, mode none:
# with no held-out unit the two modes are the same).
HEAD_MODE="$($PY -I -c "import json,sys; print(json.load(open(sys.argv[1]))['independence_mode'])" "$FREEZE")"
SECOND_MODE="$($PY -I -c "import json,sys; print(json.load(open(sys.argv[1]))['independence_second_column'])" "$FREEZE")"
mkdir -p "$TMP"
awk -F'\t' -v a="$HEAD_MODE" -v b="$SECOND_MODE" 'NR>1 && $7=="unit" {print $1 "\t" a; print $1 "\t" b}' "$UNITS" > "$TMP/tasks.list"
printf 'all\tnone\n' >> "$TMP/tasks.list"
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TMP/tasks.list")"
[ -n "$LINE" ] || { echo "task $SLURM_ARRAY_TASK_ID has no (unit, mode) pair (the array is larger than the task list); nothing to do"; exit 0; }
UNIT="$(printf '%s' "$LINE" | cut -f1)"; MODE="$(printf '%s' "$LINE" | cut -f2)"
UDIR="$OUT/$UNIT/$MODE"
[ ! -e "$UDIR" ] || fail "$UDIR exists; nothing is overwritten"
mkdir -p "$TMP/w" "$OUT"
W="$TMP/w"
echo "task $SLURM_ARRAY_TASK_ID: unit $UNIT, independence mode $MODE (headline $HEAD_MODE, second column $SECOND_MODE)"

# Proteomes: the frozen kept set and the counts of inputs.tsv.
for p in $PROTEOMES; do
  f="$PROTEOME_DIR/$p.faa"; [ -s "$f" ] || fail "$f is missing"
  want=$(awk -F'\t' -v p="$p" '$1==p {print $3}' "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/inputs.tsv")
  got=$(grep -c '^>' "$f"); [ "$want" = "$got" ] || fail "$f has $got proteins, inputs.tsv says $want"
  printf '%s\t%s\n' "$p" "$got" >> "$W/sizes.body"
done
{ printf 'proteome\tn_proteins\n'; cat "$W/sizes.body"; } > "$W/sizes.tsv"

module load hmmer/3.4 mafft/7.505 seqkit/2.4.0 MMseqs2/17-b804f
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac

# 1. Independence clusters over the members (30 percent identity, 50 percent coverage: the Task 5.3 rule), then the manifest.
mmseqs easy-cluster "$EXTRACT/ecm33_members.faa" "$W/clu" "$W/clu_tmp" --min-seq-id 0.3 -c 0.5 --threads "$THREADS" -v 1
# The HEADLINE manifest of the unit is always made (its count of excluded neighbours tells whether the second build would differ).
if [ "$UNIT" != "all" ]; then
  $BUILD manifest --members "$MEM" --unit "$UNIT" --clusters "$W/clu_cluster.tsv" --mode "$HEAD_MODE" --out "$W/head" || [ "$?" -eq 3 ] || fail "headline manifest failed"
fi
set +e
$BUILD manifest --members "$MEM" --unit "$UNIT" --clusters "$W/clu_cluster.tsv" --mode "$MODE" --out "$W/unit"
RC=$?
set -e
if [ "$RC" -eq 3 ]; then
  mkdir -p "$UDIR"; cp "$W/unit/manifest.json" "$UDIR/"; echo "$UNIT ($MODE): not_buildable (see manifest.json); no model is built"; exit 0
elif [ "$RC" -ne 0 ]; then fail "manifest failed (exit $RC)"; fi
# The second column is not rebuilt when the headline rule removed no neighbour: the two training sets are then the same.
if [ "$UNIT" != "all" ] && [ "$MODE" != "$HEAD_MODE" ]; then
  N_EXCL="$($PY -I -c "import json,sys; print(json.load(open(sys.argv[1]))['n_excluded_neighbours'])" "$W/head/manifest.json")"
  if [ "$N_EXCL" = "0" ]; then
    mkdir -p "$UDIR"; cp "$W/unit/manifest.json" "$UDIR/"; echo "the headline rule removed no neighbour of $UNIT: same training set; not rebuilt" > "$UDIR/SAME_AS_HEADLINE"
    echo "$UNIT ($MODE): same training set as the headline ($HEAD_MODE); not rebuilt"; exit 0
  fi
fi

# 2. Alignment: L-INS-i (the precedent's options), then the manifest check (the held-out unit and its 30 percent neighbours are absent)
#    and the pre-registered occupancy check of the freeze file. A failure stops the unit; no other aligner is tried.
seqkit grep -f "$W/unit/train_ids.txt" "$EXTRACT/ecm33_members.faa" > "$W/train.faa"
[ "$(grep -c '^>' "$W/train.faa")" -eq "$(wc -l < "$W/unit/train_ids.txt")" ] || fail "train.faa does not hold every training id"
mafft --localpair --maxiterate 1000 --quiet --thread "$THREADS" "$W/train.faa" > "$W/train.aln.fa"
$BUILD check-alignment --manifest "$W/unit/manifest.json" --alignment "$W/train.aln.fa" --freeze "$FREEZE" || fail "alignment check failed for $UNIT"
hmmbuild --amino --informat afa --cpu "$THREADS" -n "ecm33_$UNIT" "$W/model.hmm" "$W/train.aln.fa" > "$W/hmmbuild.log"

# 3. The false-hit side: the four kept proteomes and the shared hard negatives. FLOOR is the freeze's cutoff_floor.
FLOOR="$($PY -I -c "import json,sys; print(json.load(open(sys.argv[1]))['cutoff_floor'])" "$FREEZE")"
awk -F'\t' 'NR>1 {print $1 "\t" $6}' "$PROJ_ROOT/$NEG/controls.tsv" > "$W/neg_groups.body"
{ printf 'accession\tgroup\n'; cat "$W/neg_groups.body"; } > "$W/neg_groups.tsv"
echo Ecm33 > "$W/pf.names"   # the model name of PF12454 (analysis/wall_fusion_report/wall_fusion_models.tsv); hmmfetch -f looks up by NAME
hmmfetch -f "$RESOLVED" "$W/pf.names" > "$W/pf12454.hmm"
grep -q '^ACC *PF12454' "$W/pf12454.hmm" || fail "the fetched model is not PF12454"
: > "$W/labelled.txt"
awk '/^>/ {keep = ($1 == ">Q5AGC4")} keep' "$EXTRACT/ecm33_members.faa" > "$W/q5agc4.faa"
[ -s "$W/q5agc4.faa" ] || fail "Q5AGC4 is not in ecm33_members.faa"
{ printf 'proteome\tid\tscore\n'; } > "$W/proteome_hits_custom.tsv"; cp "$W/proteome_hits_custom.tsv" "$W/proteome_hits_relaxed.tsv"
for p in $PROTEOMES; do
  hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/c_$p.tbl" "$W/model.hmm" "$PROTEOME_DIR/$p.faa"
  hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/r_$p.tbl" "$W/pf12454.hmm" "$PROTEOME_DIR/$p.faa"
  awk -v p="$p" 'BEGIN{OFS="\t"} !/^#/ {print p, $1, $6}' "$W/c_$p.tbl" >> "$W/proteome_hits_custom.tsv"
  awk -v p="$p" 'BEGIN{OFS="\t"} !/^#/ {print p, $1, $6}' "$W/r_$p.tbl" >> "$W/proteome_hits_relaxed.tsv"
  # labelled family members of the proteome: the same definition of the family as the members (phmmer of Q5AGC4, E 1e-5)
  phmmer -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$W/l_$p.tbl" "$W/q5agc4.faa" "$PROTEOME_DIR/$p.faa"
  awk '!/^#/ {print $1}' "$W/l_$p.tbl" >> "$W/labelled.txt"
done
sort -u "$W/labelled.txt" -o "$W/labelled.txt"
cat "$PROJ_ROOT/$NEG/sequences.faa" > "$W/negatives.faa"
hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/c_neg.tbl" "$W/model.hmm" "$W/negatives.faa"
hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/r_neg.tbl" "$W/pf12454.hmm" "$W/negatives.faa"
{ printf 'accession\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/c_neg.tbl"; } > "$W/neg_hits_custom.tsv"
{ printf 'accession\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/r_neg.tbl"; } > "$W/neg_hits_relaxed.tsv"

# 4. The cutoffs (custom and relaxed Pfam): no positive score is an input. Both refuse a missing or unfrozen freeze file.
$BUILD cutoff --freeze "$FREEZE" --proteome-hits "$W/proteome_hits_custom.tsv" --labelled "$W/labelled.txt" \
  --neg-groups "$W/neg_groups.tsv" --neg-hits "$W/neg_hits_custom.tsv" --sizes "$W/sizes.tsv" --out "$W/cutoff.json"
$BUILD cutoff --freeze "$FREEZE" --proteome-hits "$W/proteome_hits_relaxed.tsv" --labelled "$W/labelled.txt" \
  --neg-groups "$W/neg_groups.tsv" --neg-hits "$W/neg_hits_relaxed.tsv" --sizes "$W/sizes.tsv" --out "$W/relaxed_cutoff.json"

# 5. ONLY NOW the held-out members (not for `all`).
if [ "$UNIT" != "all" ]; then
  [ -s "$W/cutoff.json" ] && [ -s "$W/relaxed_cutoff.json" ] || fail "no cutoff: the held-out step does not start"
  seqkit grep -f "$W/unit/heldout_ids.txt" "$EXTRACT/ecm33_members.faa" > "$W/heldout.faa"
  hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_custom.tbl" "$W/model.hmm" "$W/heldout.faa"
  hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_ga.tbl" "$W/pf12454.hmm" "$W/heldout.faa"
  hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_relaxed.tbl" "$W/pf12454.hmm" "$W/heldout.faa"
  { printf 'id\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/h_custom.tbl"; } > "$W/hits_custom.tsv"
  { printf 'id\n'; awk '!/^#/ {print $1}' "$W/h_ga.tbl"; } > "$W/hits_pfam_ga.tsv"
  { printf 'id\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/h_relaxed.tbl"; } > "$W/hits_pfam_relaxed.tsv"
  # the iteration comparator: jackhmmer seeded by the training alignment, over the held-out members, the proteomes and the negatives
  esl-reformat stockholm "$W/train.aln.fa" > "$W/train.sto"
  cat "$W/heldout.faa" "$W/negatives.faa" $(for p in $PROTEOMES; do echo "$PROTEOME_DIR/$p.faa"; done) > "$W/iter_target.faa"
  jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_iter.tbl" "$W/train.sto" "$W/iter_target.faa"
  { printf 'id\tevalue\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $5}' "$W/h_iter.tbl" | awk -F'\t' 'NR==FNR {h[$1]=1; next} ($1 in h)' "$W/unit/heldout_ids.txt" -; } > "$W/hits_iter.tsv"
fi

# 6. Provenance, then one atomic move into $OUT/<unit>.
cp "$W/neg_hits_custom.tsv" "$W/neg_hits.tsv"; cp "$W/proteome_hits_custom.tsv" "$W/proteome_hits.tsv"; cp "$W/neg_hits_relaxed.tsv" "$W/neg_hits_relaxed_out.tsv"
{
  echo "job: pr6_ecm33_build (unit $UNIT, independence mode $MODE; headline $HEAD_MODE)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "freeze sha256: $(sha256sum < "$FREEZE" | cut -d' ' -f1)"
  echo "hmmer: $(hmmsearch -h | sed -n 2p)"
  echo "mafft: $(mafft --version 2>&1 | head -1)"
  echo "mmseqs: $(mmseqs version)"
  echo "pfam_resolved: $RESOLVED"
  echo "ecm33_members.faa sha256: $(sha256sum < "$EXTRACT/ecm33_members.faa" | cut -d' ' -f1)"
  echo "model.hmm sha256: $(sha256sum < "$W/model.hmm" | cut -d' ' -f1)"
} > "$W/provenance.txt"
mkdir -p "$(dirname "$UDIR")" "$UDIR.tmp.$$"
cp "$W/unit/manifest.json" "$W/unit/train_ids.txt" "$W/unit/heldout_ids.txt" "$W/train.aln.fa" "$W/model.hmm" "$W/hmmbuild.log" \
   "$W/cutoff.json" "$W/relaxed_cutoff.json" "$W/neg_hits.tsv" "$W/proteome_hits.tsv" "$W/provenance.txt" "$UDIR.tmp.$$/"
cp "$W/neg_hits_relaxed_out.tsv" "$UDIR.tmp.$$/neg_hits_relaxed.tsv"
if [ "$UNIT" != "all" ]; then cp "$W"/hits_custom.tsv "$W"/hits_pfam_ga.tsv "$W"/hits_pfam_relaxed.tsv "$W"/hits_iter.tsv "$UDIR.tmp.$$/"; fi
mv "$UDIR.tmp.$$" "$UDIR"
echo "wrote $UDIR"
