#!/bin/bash -l
# Independent-set evaluation of the frozen custom Ecm33 'all' model, step 2 (DRAFT, NOT SUBMITTED; nothing was run on HPCC): every search on
# the test sequences with FIXED thresholds. The custom model and its cutoff are the frozen ones of the build; nothing is retuned and no
# held-out score reaches a cutoff. It writes no call, no status and no tier; the custom model is report-only.
# Design: analysis/wall_fusion_truth/pr6_ecm33_eval_freeze_PROPOSED.json (PROPOSALS until the owner signs pr6_ecm33_eval_freeze.json).
#
# Methods (one job; the test set is 187 sequences, the heavy part is the iterative searches over the background):
#   custom                    hmmsearch -T 0.0 with the committed all-model (sha256 checked); the report applies 29.6 bits
#   pfam_ga                   PF12454 hmmsearch --cut_ga
#   pfam_relaxed              PF12454 hmmsearch -T 0.0; the report applies the relaxed cutoff (15.0, from the build record)
#   iteration_alignment_seeded  jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 seeded by the build alignment (the build job's options plus -Z 70198728, owner decision 2026-10-23)
#   phmmer_seed_only          phmmer -E 1e-5 --domE 1e-5 -Z 70198728 of Q5AGC4 alone
#   jackhmmer_seed_only       jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 -Z 70198728 of Q5AGC4 alone
# The iterative and phmmer methods search the test sequences plus a background (the 161 shared negatives and the four kept proteomes, the
# build's iteration target); the training members are NOT in any target. Only test ids are written to the hit tables.
#
# THE FREEZE COMES FIRST ($PROJ_ROOT/analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json with "frozen": true). Each result file's sha256 and the
# freeze sha256 go to run.json; the report refuses a result directory without them.
# Needs the owner's go-ahead (D21), after pr6_ecm33_eval_extract.sbatch has run and its output has been read.
#
# Inputs:
#   $PROJ_ROOT    a clone at CODE_COMMIT (scripts, freeze, committed model data/custom_models/ecm33_all_2026-10-23, controls, inputs.tsv)
#   $EXTRACT      the OUTPUT directory of pr6_ecm33_eval_extract.sbatch (test_clusters.faa and extract.sha256)
#   $BUILD_ALL    the build's all/none directory on HPCC (train.aln.fa), default
#                 /bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/ecm33_build/all/none
#                 Its alignment ids must be the 481 training ids of the committed manifest; its sha256 is recorded, not pinned (UNCONFIRMED: no committed hash).
#   $PROTEOME_DIR the four kept proteomes (as the build job)
#   PFAM_HMM default /srv/projects/db/pfam/current/Pfam-A.hmm, must resolve into a directory named *Pfam$PFAM_RELEASE
# Modules: hmmer/3.4 (esl-reformat from Easel is used as in the build job).
# Resources: NOT MEASURED, a guess. The build's iteration step ran in the build job (8 cores, 12 h limit); here the target is the same size
# and three searches are iterative, so: batch partition (the short partition's limit is not known to fit), 8 cores, 16 GB, 6 hours.
#
# Submit (after the go-ahead; $OUT must NOT exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,EXTRACT=/path/to/extract_out,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_eval_search.sbatch
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH -J pr6_ecm33_eval_search
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT before sbatch}"
: "${EXTRACT:?export EXTRACT (the output of pr6_ecm33_eval_extract.sbatch) before sbatch}"
: "${OUT:?export OUT before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
DATA_ROOT="${DATA_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
BUILD_ALL="${BUILD_ALL:-$DATA_ROOT/_workdir/kiss_cluster_split/ecm33_build/all/none}"
PROTEOME_DIR="${PROTEOME_DIR:-$DATA_ROOT/_workdir/sorting_hat}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_ecm33_eval_search.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-8}"
WF=analysis/wall_fusion_truth
NEG=data/controls/hard-negatives-shared
CTRL=data/controls/wall-fusion-kiss2019
MODEL_DIR=data/custom_models/ecm33_all_2026-10-23
BLD=docs/reports/data/pr6_ecm33_build_2026-10-23
FREEZE="$PROJ_ROOT/$WF/pr6_ecm33_eval_freeze.json"
PROTEOMES="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"
PY=/usr/bin/python3.12
FLOOR=0.0

fail() { echo "FATAL: $*" >&2; exit 1; }
jget() { $PY -I -c "import json,sys; d=json.load(open(sys.argv[1])); [d:=d[k] for k in sys.argv[2].split('.')]; print(d)" "$FREEZE" "$1"; }

git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$NEG" "$CTRL" "$MODEL_DIR" "$BLD" scripts/sorting_hat/refresh_3_10/inputs.tsv \
  || fail "$WF, $NEG, $CTRL, $MODEL_DIR, $BLD or inputs.tsv changed since $CODE_COMMIT"
[ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$NEG" "$CTRL" "$MODEL_DIR" "$BLD")" ] || fail "uncommitted changes in a pinned directory"
echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
[ -f "$FREEZE" ] || fail "$FREEZE is missing: the freeze commit comes before any search"
grep -q '"frozen": true' "$FREEZE" || fail "$FREEZE is not frozen"
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -s "$EXTRACT/test_clusters.faa" ] || fail "$EXTRACT is not an extract output"
( cd "$EXTRACT" && sha256sum -c extract.sha256 ) || fail "the extract output differs from extract.sha256"
mkdir -p "$TMP"
W="$TMP"

# The frozen model and cutoffs: committed files, hashes and numbers checked against the freeze.
MODEL="$PROJ_ROOT/$MODEL_DIR/model.hmm"
[ "$(sha256sum < "$MODEL" | cut -d' ' -f1)" = "$(jget custom_model.sha256)" ] || fail "model.hmm is not the model of the freeze (sha256)"
[ "$(sha256sum < "$PROJ_ROOT/$MODEL_DIR/cutoff.json" | cut -d' ' -f1)" = "$(jget custom_model.cutoff_json_sha256)" ] || fail "cutoff.json differs from the freeze"
$PY -I -c "
import json,sys
fz=json.load(open(sys.argv[1])); c=json.load(open(sys.argv[2])); r=json.load(open(sys.argv[3]))
assert float(c['cutoff'])==float(fz['custom_model']['cutoff_bits']), 'custom cutoff differs from the freeze'
assert float(r['cutoff'])==float(fz['searches']['relaxed_pfam_cutoff_bits']), 'relaxed Pfam cutoff differs from the freeze'
" "$FREEZE" "$PROJ_ROOT/$MODEL_DIR/cutoff.json" "$PROJ_ROOT/$BLD/ecm33_build/all/none/relaxed_cutoff.json" || fail "cutoff check"
EVALDB="$(jget searches.evalue_db_size)"
[ "$EVALDB" = "70198728" ] || fail "evalue_db_size of the freeze is $EVALDB, not the UniRef50 sequence count 70198728"

# The alignment of the build for the alignment-seeded comparator: its ids must be exactly the committed training ids.
[ -s "$BUILD_ALL/train.aln.fa" ] || fail "$BUILD_ALL/train.aln.fa is missing"
$PY -I -c "
import json,sys
want={t['id'] for t in json.load(open(sys.argv[1]))['training']}
got={l[1:].split()[0] for l in open(sys.argv[2]) if l.startswith('>')}
assert got==want, f'alignment ids differ from the committed manifest: {len(got)} vs {len(want)}'
" "$PROJ_ROOT/$BLD/ecm33_build/all/none/manifest.json" "$BUILD_ALL/train.aln.fa" || fail "alignment ids check"

module load hmmer/3.4
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
echo Ecm33 > "$W/pf.names"   # the model name of PF12454; hmmfetch -f looks up by NAME (as the build job)
hmmfetch -f "$RESOLVED" "$W/pf.names" > "$W/pf12454.hmm"
grep -q '^ACC *PF12454' "$W/pf12454.hmm" || fail "the fetched model is not PF12454"

# Targets: the test sequences; and the test sequences plus the background (negatives and the four kept proteomes; never the training members).
TEST="$EXTRACT/test_clusters.faa"
cp "$PROJ_ROOT/$NEG/sequences.faa" "$W/negatives.faa"
for p in $PROTEOMES; do [ -s "$PROTEOME_DIR/$p.faa" ] || fail "$PROTEOME_DIR/$p.faa is missing"; done
cat "$TEST" "$W/negatives.faa" $(for p in $PROTEOMES; do echo "$PROTEOME_DIR/$p.faa"; done) > "$W/bg_target.faa"
awk '/^>/ {keep = ($1 == ">Q5AGC4")} keep' "$PROJ_ROOT/$CTRL/sequences.faa" > "$W/q5agc4.faa"
[ "$(grep -c '^>' "$W/q5agc4.faa")" -eq 1 ] || fail "Q5AGC4 is not exactly once in $CTRL/sequences.faa"
esl-reformat stockholm "$BUILD_ALL/train.aln.fa" > "$W/train.sto"
grep '^>' "$TEST" | sed 's/^>//' > "$W/test_ids.txt"

# Searches.
hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/custom.tbl" "$MODEL" "$TEST"
hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --tblout "$W/ga.tbl" "$W/pf12454.hmm" "$TEST"
hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/rel.tbl" "$W/pf12454.hmm" "$TEST"
jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 -Z "$EVALDB" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/iter.tbl" "$W/train.sto" "$W/bg_target.faa"
phmmer -E 1e-5 --domE 1e-5 -Z "$EVALDB" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/phmmer.tbl" "$W/q5agc4.faa" "$W/bg_target.faa"
jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 -Z "$EVALDB" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/jack.tbl" "$W/q5agc4.faa" "$W/bg_target.faa"

# Hit tables: test ids only. tblout columns: 1 target, 5 full-sequence E-value, 6 full-sequence score.
keep() { awk -F'\t' 'NR==FNR {t[$1]=1; next} ($1 in t)' "$W/test_ids.txt" -; }
mkdir -p "$OUT.tmp.$$"
O="$OUT.tmp.$$"
{ printf 'id\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/custom.tbl" | keep; } > "$O/hits_custom.tsv"
{ printf 'id\n'; awk '!/^#/ {print $1}' "$W/ga.tbl" | keep; } > "$O/hits_pfam_ga.tsv"
{ printf 'id\tscore\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $6}' "$W/rel.tbl" | keep; } > "$O/hits_pfam_relaxed.tsv"
{ printf 'id\tevalue\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $5}' "$W/iter.tbl" | keep; } > "$O/hits_iter_aln.tsv"
{ printf 'id\tevalue\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $5}' "$W/phmmer.tbl" | keep; } > "$O/hits_phmmer_seed.tsv"
{ printf 'id\tevalue\n'; awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $5}' "$W/jack.tbl" | keep; } > "$O/hits_jack_seed.tsv"

# run.json: the freeze sha256 and the sha256 of every result file (the report refuses without them).
FILES="hits_custom.tsv hits_pfam_ga.tsv hits_pfam_relaxed.tsv hits_iter_aln.tsv hits_phmmer_seed.tsv hits_jack_seed.tsv"
FSHA="$(sha256sum < "$FREEZE" | cut -d' ' -f1)"
ASHA="$(sha256sum < "$BUILD_ALL/train.aln.fa" | cut -d' ' -f1)"
( cd "$O" && $PY -I -c "
import hashlib,json,sys
files={n: hashlib.sha256(open(n,'rb').read()).hexdigest() for n in sys.argv[5:]}
json.dump({'freeze_sha256': sys.argv[1], 'commit': sys.argv[2], 'train_aln_sha256': sys.argv[3], 'hmmer': sys.argv[4], 'files': files}, open('run.json','w'), indent=1, sort_keys=True)
" "$FSHA" "$(git -C "$PROJ_ROOT" rev-parse HEAD)" "$ASHA" "$(hmmsearch -h | sed -n 2p)" $FILES )
{
  echo "job: pr6_ecm33_eval_search (independent-set evaluation, ecm33 all model)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "freeze sha256: $FSHA"
  echo "model.hmm sha256: $(sha256sum < "$MODEL" | cut -d' ' -f1)"
  echo "pfam_resolved: $RESOLVED"
  echo "train.aln.fa sha256: $ASHA (recorded, not pinned)"
  echo "test_clusters.faa sha256: $(sha256sum < "$TEST" | cut -d' ' -f1)"
  echo "evalue_db_size (-Z): $EVALDB"
} > "$O/provenance.txt"
mv "$O" "$OUT"
echo "wrote $OUT"
