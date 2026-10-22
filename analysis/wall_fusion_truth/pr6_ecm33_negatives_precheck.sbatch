#!/bin/bash -l
# PR 6 ecm33 pre-freeze check (DRAFT, NOT SUBMITTED; nothing of this has run on HPCC): which of the 161 shared hard negatives are
# members of the Ecm33 family? `phmmer` of the Kiss Ecm33 query Q5AGC4 against the committed shared negatives, at E <= 1e-5 (the
# definition of 'family member' of the freeze's labelled_rule) plus a separate descriptive pass at E <= 1e-3.
# The result is the complete list for `negatives_ignored` in the freeze (P38248 is listed; Q4WNS8 is a second known member; a member
# under another name cannot be excluded by reading names). It builds no model, scores no held-out clade and writes no call, status
# or tier. The freeze file is NOT read or changed by this job.
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update and this job. A tiny job (one query, 161
# sequences, seconds), but still a SLURM job, never on the login node. $SCRATCH is the temporary directory (nothing is written outside it and $OUT).
#
# Inputs:
#   $IN/kiss_sequences.faa (+ .sha256)   already staged at .../_workdir/kiss_cluster_split/in/ (no new staging); Q5AGC4 is taken from it
#   $PROJ_ROOT (a clone at CODE_COMMIT)  holds data/controls/hard-negatives-shared/{sequences.faa,controls.tsv} and the scripts
# Module: hmmer/3.4 (listed by `module avail` on HPCC 2026-10-22).
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_negatives_precheck.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --time=0:30:00
#SBATCH -J pr6_ecm33_negcheck
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
QUERY=Q5AGC4
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_ecm33_negcheck.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-2}"
WF=analysis/wall_fusion_truth
NEG=data/controls/hard-negatives-shared

fail() { echo "FATAL: $*" >&2; exit 1; }

# The code and the negative set must be the reviewed ones.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$NEG" \
    || fail "$WF or $NEG changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$NEG")" ] || fail "uncommitted changes in $WF or $NEG"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

# sidecar_ok FILE: the .sha256 next to FILE holds the sha256 of FILE
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -f "$IN/kiss_sequences.faa" ] && [ -f "$IN/kiss_sequences.faa.sha256" ] || fail "$IN/kiss_sequences.faa or its .sha256 is missing"
sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa differs from its .sha256"
for f in sequences.faa controls.tsv; do [ -s "$PROJ_ROOT/$NEG/$f" ] || fail "$PROJ_ROOT/$NEG/$f is missing or empty"; done
mkdir -p "$TMP"

# The negatives: every sequence name is an accession of controls.tsv, and the reverse (the ids stay exactly as in the file).
grep '^>' "$PROJ_ROOT/$NEG/sequences.faa" | sed 's/^>//; s/[[:space:]].*//' | sort > "$TMP/neg_fa_ids.txt"
tail -n +2 "$PROJ_ROOT/$NEG/controls.tsv" | cut -f1 | sort > "$TMP/neg_ctl_ids.txt"
cmp -s "$TMP/neg_fa_ids.txt" "$TMP/neg_ctl_ids.txt" || fail "the sequence names of $NEG/sequences.faa differ from the accessions of controls.tsv"
N_NEG=$(wc -l < "$TMP/neg_fa_ids.txt")
[ "$N_NEG" -eq 161 ] || fail "$NEG has $N_NEG sequences, not 161"

# The query: exact match on the first token of the header; it must exist exactly once.
awk -v q="$QUERY" '/^>/ { name = substr($1, 2); keep = (name == q) } keep { print }' "$IN/kiss_sequences.faa" > "$TMP/query.faa"
[ "$(grep -c '^>' "$TMP/query.faa")" -eq 1 ] || fail "$QUERY is not in $IN/kiss_sequences.faa exactly once"
echo "inputs ok: query $QUERY, $N_NEG shared negatives"

module load hmmer/3.4

# 1. The definition of the freeze: E <= 1e-5. 2. A separate DESCRIPTIVE pass at E <= 1e-3 (never part of the proposed list).
phmmer -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/phmmer_E1e-5.tbl" "$TMP/query.faa" "$PROJ_ROOT/$NEG/sequences.faa"
tail -n 3 "$TMP/phmmer_E1e-5.tbl" | grep -qx '# \[ok\]' || fail "phmmer_E1e-5.tbl has no '# [ok]' trailer"
phmmer -E 1e-3 --domE 1e-3 --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/phmmer_E1e-3.tbl" "$TMP/query.faa" "$PROJ_ROOT/$NEG/sequences.faa"
tail -n 3 "$TMP/phmmer_E1e-3.tbl" | grep -qx '# \[ok\]' || fail "phmmer_E1e-3.tbl has no '# [ok]' trailer"

# 3. The tables (standard library only). The summary script refuses a non-empty output directory.
/usr/bin/python3.12 -I "$PROJ_ROOT/$WF/pr6_ecm33_negatives_precheck.py" summarize \
  --strict "$TMP/phmmer_E1e-5.tbl" --weak "$TMP/phmmer_E1e-3.tbl" \
  --controls "$PROJ_ROOT/$NEG/controls.tsv" --negatives-fasta "$PROJ_ROOT/$NEG/sequences.faa" \
  --already-ignored P38248 --out "$TMP/tables"

{
  echo "job: pr6_ecm33_negcheck (pre-freeze check of the ecm33 hard negatives)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(phmmer -h | sed -n 2p)"
  echo "python: $(/usr/bin/python3.12 --version)"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "input shared negatives sequences.faa sha256 $(sha256sum < "$PROJ_ROOT/$NEG/sequences.faa" | cut -d' ' -f1)"
  echo "input shared negatives controls.tsv sha256 $(sha256sum < "$PROJ_ROOT/$NEG/controls.tsv" | cut -d' ' -f1)"
  echo "query: $QUERY"
  echo "settings: phmmer -E 1e-5 --domE 1e-5 (the definition); phmmer -E 1e-3 --domE 1e-3 (descriptive)"
} > "$TMP/tables/provenance.txt"
cp "$TMP/phmmer_E1e-5.tbl" "$TMP/phmmer_E1e-3.tbl" "$TMP/tables/"

mkdir -p "$(dirname "$OUT")"
mv "$TMP/tables" "$OUT"
echo "wrote $OUT"
