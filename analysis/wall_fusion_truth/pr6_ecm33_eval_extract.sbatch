#!/bin/bash -l
# Independent-set evaluation of the frozen custom Ecm33 'all' model, step 1 (DRAFT, NOT SUBMITTED; nothing was run on HPCC): select the test
# clusters from the committed candidates.tsv, fetch the cluster REPRESENTATIVE sequences from UniRef50 2025_03 by id, and run the cluster30
# separation (MMseqs2) of the test sequences against the 481 training members. It builds no model, scores nothing, and writes no call, no
# status and no tier. Sequences stay in $OUT on HPCC; none is committed.
# Design: analysis/wall_fusion_truth/pr6_ecm33_eval_freeze_PROPOSED.json (every choice there is a PROPOSAL until the owner signs
# pr6_ecm33_eval_freeze.json); README_pr6_ecm33_eval.md.
#
# THE FREEZE COMES FIRST: $PROJ_ROOT/analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json must exist and say "frozen": true.
# Needs the owner's go-ahead for EACH step before it is submitted (D21). Heavy work is a SLURM job, never on the login node. $SCRATCH is the
# temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $PROJ_ROOT   a clone at CODE_COMMIT (candidates.tsv, members.tsv, the freeze, the scripts; all committed; the clone must hold the
#                merge of the PR of this script).
#   $UNIREF50    default /srv/projects/db/Uniprot/2025_03/uniref50.fasta; its size must be 28,781,193,344 bytes (the round 2 provenance).
#   $TRAIN_FAA   the staged ecm33_members.faa of the build's extract job (481 training members); its sha256 must be the freeze's
#                separation.training_members_faa_sha256. Default: /bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/ecm33_extract/ecm33_members.faa
# Tools: seqkit/2.4.0 (the UniRef50 file is streamed ONCE), MMseqs2/17-b804f (`mmseqs easy-cluster --min-seq-id 0.3 -c 0.5`, from the freeze).
# Resources: NOT MEASURED, a guess. batch partition (the file-system pass over 28.8 GB is the long part; the 2-hour limit of a short
# partition is not assured for it), 8 cores, 16 GB (187 ids and 668 sequences are tiny; seqkit buffers little), 3 hours.
#
# Submit (after the go-ahead; $OUT must NOT exist; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,OUT=/path/to/new_outdir \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_eval_extract.sbatch
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH -J pr6_ecm33_eval_extract
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
UNIREF50="${UNIREF50:-/srv/projects/db/Uniprot/2025_03/uniref50.fasta}"
UNIREF50_BYTES=28781193344
TRAIN_FAA="${TRAIN_FAA:-/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/ecm33_extract/ecm33_members.faa}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_ecm33_eval_extract.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-8}"
WF=analysis/wall_fusion_truth
CAND=docs/reports/data/ecm33_independent_set_candidates
BLD=docs/reports/data/pr6_ecm33_build_2026-10-23
FREEZE="$PROJ_ROOT/$WF/pr6_ecm33_eval_freeze.json"
PY=/usr/bin/python3.12

fail() { echo "FATAL: $*" >&2; exit 1; }
jget() { $PY -I -c "import json,sys; d=json.load(open(sys.argv[1])); [d:=d[k] for k in sys.argv[2].split('.')]; print(d)" "$FREEZE" "$1"; }

git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$CAND" "$BLD/members" || fail "$WF, $CAND or $BLD/members changed since $CODE_COMMIT"
[ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$CAND" "$BLD/members")" ] || fail "uncommitted changes in $WF, $CAND or $BLD/members"
echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
[ -f "$FREEZE" ] || fail "$FREEZE is missing: the freeze commit comes before any run"
grep -q '"frozen": true' "$FREEZE" || fail "$FREEZE is not frozen"
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -s "$UNIREF50" ] || fail "$UNIREF50 is missing or empty"
[ "$(stat -c %s "$UNIREF50")" = "$UNIREF50_BYTES" ] || fail "$UNIREF50 is not $UNIREF50_BYTES bytes: it is not the file of the round 2 record"
[ -s "$TRAIN_FAA" ] || fail "$TRAIN_FAA is missing"
[ "$(sha256sum < "$TRAIN_FAA" | cut -d' ' -f1)" = "$(jget separation.training_members_faa_sha256)" ] || fail "$TRAIN_FAA is not the training FASTA of the freeze (sha256)"
mkdir -p "$TMP"

# 1. The test set from the committed table: ids and metadata only.
$PY -I "$PROJ_ROOT/$WF/pr6_ecm33_eval_members.py" select --candidates "$PROJ_ROOT/$CAND/candidates.tsv" --out "$TMP/sel"
$PY -I "$PROJ_ROOT/$WF/pr6_ecm33_eval_members.py" purity --candidates "$PROJ_ROOT/$CAND/candidates.tsv" --members "$PROJ_ROOT/$BLD/members/members.tsv" --out "$TMP/purity"
N_IDS=$(wc -l < "$TMP/sel/test_cluster_ids.txt")
[ "$N_IDS" -eq "$(jget test_set.expected_clusters)" ] || fail "$N_IDS test clusters, the freeze says $(jget test_set.expected_clusters)"
N_PE=$(wc -l < "$TMP/sel/headline_pe_cluster_ids.txt")
[ "$N_PE" -eq "$(jget test_set.expected_headline_pe_clusters)" ] || fail "$N_PE clusters at the PANTHER E cut, the freeze says $(jget test_set.expected_headline_pe_clusters)"

# 2. One pass over UniRef50: the representatives whose first header token is exactly a cluster id (seqkit grep -f: exact ids, no regex).
module load seqkit/2.4.0 MMseqs2/17-b804f
seqkit grep -j "$THREADS" -f "$TMP/sel/test_cluster_ids.txt" "$UNIREF50" | awk '/^>/ { print $1; next } { print }' > "$TMP/test_clusters.faa"
N_FOUND=$(grep -c '^>' "$TMP/test_clusters.faa")
[ "$N_FOUND" -eq "$N_IDS" ] || fail "found $N_FOUND of $N_IDS representatives in $UNIREF50"
[ "$(grep '^>' "$TMP/test_clusters.faa" | sort | uniq -d | wc -l)" -eq 0 ] || fail "duplicate ids in the test FASTA"

# 3. The cluster30 separation: training members and test sequences clustered TOGETHER (the freeze's command).
cat "$TRAIN_FAA" "$TMP/test_clusters.faa" > "$TMP/all.faa"
[ "$(grep '^>' "$TMP/all.faa" | sort | uniq -d | wc -l)" -eq 0 ] || fail "an id is in both the training and the test FASTA"
mmseqs easy-cluster "$TMP/all.faa" "$TMP/clu" "$TMP/clu_tmp" --min-seq-id 0.3 -c 0.5 --threads "$THREADS" -v 1
grep '^>' "$TRAIN_FAA" | sed 's/^>//' > "$TMP/train_ids.txt"
$PY -I "$PROJ_ROOT/$WF/pr6_ecm33_eval_members.py" separate --clusters-tsv "$TMP/clu_cluster.tsv" --train-ids "$TMP/train_ids.txt" \
  --test-ids "$TMP/sel/test_cluster_ids.txt" --out "$TMP/separation.tsv"

# 4. Provenance and one atomic move.
mkdir -p "$OUT.tmp.$$"
cp "$TMP/test_clusters.faa" "$TMP/separation.tsv" "$OUT.tmp.$$/"
cp "$TMP/sel/test_clusters.tsv" "$TMP/sel/test_proteins.tsv" "$TMP/sel/test_cluster_ids.txt" "$TMP/sel/headline_pe_cluster_ids.txt" "$TMP/sel/counts.json" "$OUT.tmp.$$/"
cp "$TMP/purity/purity_members.tsv" "$TMP/purity/purity.json" "$OUT.tmp.$$/"
( cd "$OUT.tmp.$$" && sha256sum test_clusters.faa separation.tsv test_clusters.tsv > extract.sha256 )
{
  echo "job: pr6_ecm33_eval_extract (independent-set evaluation, ecm33 all model)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "freeze sha256: $(sha256sum < "$FREEZE" | cut -d' ' -f1)"
  echo "seqkit: $(seqkit version 2>&1 | head -1)"
  echo "mmseqs: $(mmseqs version)"
  echo "uniref50: $UNIREF50 ($UNIREF50_BYTES bytes)"
  echo "training faa: $TRAIN_FAA sha256 $(sha256sum < "$TRAIN_FAA" | cut -d' ' -f1)"
  echo "test clusters: $N_IDS (at the PANTHER E cut: $N_PE); separated: $(awk -F'\t' 'NR>1 && $5=="yes"' "$TMP/separation.tsv" | wc -l)"
} > "$OUT.tmp.$$/provenance.txt"
mv "$OUT.tmp.$$" "$OUT"
echo "wrote $OUT"
