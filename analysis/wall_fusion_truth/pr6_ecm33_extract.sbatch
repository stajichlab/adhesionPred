#!/bin/bash -l
# PR 6 steps 6.1 and 6.2 for ecm33 (DRAFT, NOT SUBMITTED; nothing of PR 6 beyond the gap search has run on HPCC): select the family members
# of Q5AGC4 from the committed round 2 record (the conservative phmmer set, pr6_ecm33_members.py) and extract their sequences from
# UniRef50 2025_03 into a staged FASTA for the build job. It builds no model and writes no call, no status and no tier.
# Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md (decisions 2, 5, 7; sections 4.1 to 4.4).
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update and this job. Heavy work is a SLURM
# job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $PROJ_ROOT (a clone at CODE_COMMIT)  holds everything else, all committed: the round 2 record
#                                        docs/reports/data/pr6_gap_uniref50_2026-10-22/gap_homologs.tsv.gz and the Kiss control set
#                                        (Q5AGC4 is read from data/controls/wall-fusion-kiss2019/sequences.faa). The clone must contain
#                                        the merge of the PR that holds the round 2 record.
#   $UNIREF50                            the UniRef50 2025_03 FASTA (default /srv/projects/db/Uniprot/2025_03/uniref50.fasta). Its size
#                                        must equal the 28,781,193,344 bytes of the round 2 provenance (the 28.8 GB file is NOT hashed
#                                        again here; the sha256 bc52bdd6... of the round 2 record is the identity of the file).
# Tool: seqkit (module seqkit/2.4.0 is listed by `module avail seqkit` on HPCC). The job streams UniRef50 ONCE.
# Resources: NOT MEASURED. The request below (batch partition, 4 cores, 8 GB, 2 hours) is a guess: one pass over 28.8 GB by seqkit
# should take minutes to a few tens of minutes at the file-system read speed, and memory is small (480 ids). The values are defaults on the
# #SBATCH lines; override on the sbatch command line (-p, --time, -c, --mem) if needed.
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,OUT=/path/to/new_outdir \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_extract.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J pr6_ecm33_extract
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
UNIREF50="${UNIREF50:-/srv/projects/db/Uniprot/2025_03/uniref50.fasta}"
UNIREF50_BYTES=28781193344  # the round 2 provenance (uniprot_fasta_bytes)
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_ecm33_extract.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
WF=analysis/wall_fusion_truth
REC=docs/reports/data/pr6_gap_uniref50_2026-10-22
CTRL=data/controls/wall-fusion-kiss2019

fail() { echo "FATAL: $*" >&2; exit 1; }

check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$REC" "$CTRL" \
    || fail "$WF, $REC or $CTRL changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$REC" "$CTRL")" ] || fail "uncommitted changes in $WF, $REC or $CTRL"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}
one_seq() { awk -v id="$2" '/^>/ { keep = ($1 == ">" id) } keep' "$1"; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
HOM="$PROJ_ROOT/$REC/gap_homologs.tsv.gz"
[ -f "$HOM" ] || fail "$HOM is missing (the clone does not hold the merged round 2 record)"
[ -s "$UNIREF50" ] || fail "$UNIREF50 is missing or empty"
[ "$(stat -c %s "$UNIREF50")" = "$UNIREF50_BYTES" ] || fail "$UNIREF50 is not $UNIREF50_BYTES bytes: it is not the file of the round 2 record"
mkdir -p "$TMP"

# 1. Members, held-out units and ids (stdlib; the committed record only).
/usr/bin/python3.12 -I "$PROJ_ROOT/$WF/pr6_ecm33_members.py" --homologs "$HOM" --out "$TMP/members"
N_IDS=$(wc -l < "$TMP/members/member_ids.txt")
[ "$N_IDS" -gt 0 ] || fail "no member ids"

# 2. One pass over UniRef50: the sequences whose first header token is exactly a member id (seqkit grep -f: exact IDs, no regex).
module load seqkit/2.4.0
seqkit grep -j "$THREADS" -f "$TMP/members/member_ids.txt" "$UNIREF50" > "$TMP/reps.faa"
N_FOUND=$(grep -c '^>' "$TMP/reps.faa")
[ "$N_FOUND" -eq "$N_IDS" ] || fail "found $N_FOUND of $N_IDS member sequences in $UNIREF50"
# The query itself (Q5AGC4, a Kiss control with a public UniProt sequence): from the committed control set, under its accession.
one_seq "$PROJ_ROOT/$CTRL/sequences.faa" Q5AGC4 > "$TMP/query.faa"
[ -s "$TMP/query.faa" ] || fail "Q5AGC4 is not in $CTRL/sequences.faa"
# One header token per record (UniRef50_<id> or Q5AGC4); the rest of the UniRef50 header (cluster name, n=, Tax=) is dropped.
{ cat "$TMP/query.faa"; cat "$TMP/reps.faa"; } | awk '/^>/ { print $1; next } { print }' > "$TMP/ecm33_members.faa"
N_FA=$(grep -c '^>' "$TMP/ecm33_members.faa")
N_MEM=$(($(wc -l < "$TMP/members/members.tsv") - 1))
[ "$N_FA" -eq "$N_MEM" ] || fail "ecm33_members.faa has $N_FA sequences but members.tsv has $N_MEM rows"
# No duplicate id.
[ "$(grep '^>' "$TMP/ecm33_members.faa" | sort | uniq -d | wc -l)" -eq 0 ] || fail "duplicate ids in ecm33_members.faa"

mkdir -p "$OUT.tmp.$$"
cp -r "$TMP/members" "$OUT.tmp.$$/members"
cp "$TMP/ecm33_members.faa" "$OUT.tmp.$$/ecm33_members.faa"
( cd "$OUT.tmp.$$" && sha256sum ecm33_members.faa | cut -d' ' -f1 > ecm33_members.faa.sha256 )
{
  echo "job: pr6_ecm33_extract (PR 6 steps 6.1 and 6.2, ecm33)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "seqkit: $(seqkit version 2>&1 | head -1)"
  echo "uniref50: $UNIREF50 ($UNIREF50_BYTES bytes, the file of the round 2 record; sha256 bc52bdd6c97ba2cc365b7e42ae7abb34865e0a1a82a2f95a74590aae040c81ff)"
  echo "homologs: $REC/gap_homologs.tsv.gz sha256 $(sha256sum < "$HOM" | cut -d' ' -f1)"
  echo "members: $N_MEM (member ids fetched from UniRef50: $N_IDS; the query Q5AGC4 from the Kiss control set)"
  echo "ecm33_members.faa sha256 $(cut -d' ' -f1 "$OUT.tmp.$$/ecm33_members.faa.sha256")"
} > "$OUT.tmp.$$/provenance.txt"
mv "$OUT.tmp.$$" "$OUT"
echo "wrote $OUT ($N_MEM members)"
