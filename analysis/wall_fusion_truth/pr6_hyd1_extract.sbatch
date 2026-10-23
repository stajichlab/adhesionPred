#!/bin/bash -l
# PR 6 for hyd1_like (DRAFT, NOT SUBMITTED; nothing of this has run on HPCC). Modeled on pr6_ecm33_extract.sbatch: select the members of A0ACE4KFX7
# from the committed drift record (pr6_hyd1_members.py), extract their UniRef50 sequences, cluster the members at 30 percent identity / 50 percent
# coverage (MMseqs2), and list the leave-cluster-out units. It builds no model and writes no call, status or tier. No freeze is needed: nothing is scored.
# Needs the owner's go-ahead for EACH step before it is submitted (D21). Heavy work is a SLURM job. $SCRATCH is the temporary directory.
# Inputs: $PROJ_ROOT (a clone at CODE_COMMIT) with docs/reports/data/pr6_drift_check_2026-10-22/{drift_hits,hits_manifest}.tsv and the Kiss control set
#   (A0ACE4KFX7 is read from data/controls/wall-fusion-kiss2019/sequences.faa); $UNIREF50 (default the 2025_03 file, size checked).
# Modules: seqkit/2.4.0, MMseqs2/17-b804f. Resources NOT MEASURED (about 30 members; one pass over 28.8 GB by seqkit): a guess.
# Submit: sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,OUT=/path/to/new_outdir \
#   -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_hyd1_extract.sbatch    ($OUT must NOT exist)
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J pr6_hyd1_extract
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
UNIREF50="${UNIREF50:-/srv/projects/db/Uniprot/2025_03/uniref50.fasta}"
UNIREF50_BYTES=28781193344
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_hyd1_extract.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
WF=analysis/wall_fusion_truth
DRIFT=docs/reports/data/pr6_drift_check_2026-10-22
CTRL=data/controls/wall-fusion-kiss2019
QUERY=A0ACE4KFX7
PY=/usr/bin/python3.12
fail() { echo "FATAL: $*" >&2; exit 1; }
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$DRIFT" "$CTRL" || fail "$WF, $DRIFT or $CTRL changed since $CODE_COMMIT"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$DRIFT" "$CTRL")" ] || fail "uncommitted changes in $WF, $DRIFT or $CTRL"
}
check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
for f in drift_hits.tsv hits_manifest.tsv; do [ -s "$PROJ_ROOT/$DRIFT/$f" ] || fail "$DRIFT/$f is missing"; done
[ -s "$UNIREF50" ] && [ "$(stat -c %s "$UNIREF50")" = "$UNIREF50_BYTES" ] || fail "$UNIREF50 is missing or not the $UNIREF50_BYTES-byte file of the round 2 record"
mkdir -p "$TMP"

# 1. Members and clade units (stdlib; the committed record only). The reading of the training rule is the one in the freeze PROPOSAL/freeze: pass READING.
READING="${READING:-total}"
$PY -I "$PROJ_ROOT/$WF/pr6_hyd1_members.py" --drift-hits "$PROJ_ROOT/$DRIFT/drift_hits.tsv" --hits-manifest "$PROJ_ROOT/$DRIFT/hits_manifest.tsv" \
  --reading "$READING" --out "$TMP/members"
N_IDS=$(wc -l < "$TMP/members/member_ids.txt"); [ "$N_IDS" -gt 0 ] || fail "no member ids"

# 2. One pass over UniRef50 for the member ids (exact ids); the query from the Kiss control set.
module load seqkit/2.4.0 MMseqs2/17-b804f
seqkit grep -j "$THREADS" -f "$TMP/members/member_ids.txt" "$UNIREF50" > "$TMP/reps.faa"
[ "$(grep -c '^>' "$TMP/reps.faa")" -eq "$N_IDS" ] || fail "not every member id was found in $UNIREF50"
awk -v q="$QUERY" '/^>/ { keep = (substr($1, 2) == q) } keep' "$PROJ_ROOT/$CTRL/sequences.faa" > "$TMP/query.faa"
[ "$(grep -c '^>' "$TMP/query.faa")" -eq 1 ] || fail "$QUERY is not in $CTRL/sequences.faa exactly once"
{ cat "$TMP/query.faa" "$TMP/reps.faa"; } | awk '/^>/ { print $1; next } { print }' > "$TMP/hyd1_members.faa"
N_MEM=$(($(wc -l < "$TMP/members/members.tsv") - 1))
[ "$(grep -c '^>' "$TMP/hyd1_members.faa")" -eq "$N_MEM" ] || fail "hyd1_members.faa and members.tsv differ in size"
[ "$(grep '^>' "$TMP/hyd1_members.faa" | sort | uniq -d | wc -l)" -eq 0 ] || fail "duplicate ids"

# 3. Independence clusters over the members, then the leave-cluster-out units.
mmseqs easy-cluster "$TMP/hyd1_members.faa" "$TMP/clu" "$TMP/clu_tmp" --min-seq-id 0.3 -c 0.5 --threads "$THREADS" -v 1
$PY -I "$PROJ_ROOT/$WF/pr6_hyd1_build.py" cluster-units --members "$TMP/members/members.tsv" --clusters "$TMP/clu_cluster.tsv" --reading "$READING" --out "$TMP/units_cluster.tsv"

mkdir -p "$OUT.tmp.$$"
cp -r "$TMP/members" "$OUT.tmp.$$/members"
cp "$TMP/hyd1_members.faa" "$TMP/clu_cluster.tsv" "$TMP/units_cluster.tsv" "$OUT.tmp.$$/"
( cd "$OUT.tmp.$$" && sha256sum hyd1_members.faa | cut -d' ' -f1 > hyd1_members.faa.sha256 )
{
  echo "job: pr6_hyd1_extract; commit $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned $CODE_COMMIT)"
  echo "seqkit: $(seqkit version 2>&1 | head -1); mmseqs: $(mmseqs version)"
  echo "uniref50: $UNIREF50 ($UNIREF50_BYTES bytes; round 2 sha256 bc52bdd6c97ba2cc365b7e42ae7abb34865e0a1a82a2f95a74590aae040c81ff)"
  echo "drift_hits.tsv sha256 $(sha256sum < "$PROJ_ROOT/$DRIFT/drift_hits.tsv" | cut -d' ' -f1)"
  echo "members: $N_MEM (query $QUERY + $N_IDS UniRef50 representatives); reading $READING"
  echo "hyd1_members.faa sha256 $(cut -d' ' -f1 "$OUT.tmp.$$/hyd1_members.faa.sha256")"
} > "$OUT.tmp.$$/provenance.txt"
mv "$OUT.tmp.$$" "$OUT"
echo "wrote $OUT ($N_MEM members)"
