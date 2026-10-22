#!/bin/bash -l
# PR 6 drift check (DRAFT, NOT SUBMITTED; nothing of this has run on HPCC): for the queries whose outside homologs come
# only from jackhmmer iteration (aga1 P32323, ece1 Q07730 and, after the clade-bucket fix, hyd1_like A0ACE4KFX7) ask of every
# jackhmmer outside fungal hit whether a direct pairwise alignment to the query supports it or composition explains it.
# It builds no model and writes no call, no status and no tier. Design, freeze and classes: pr6_drift_check.py (docstring).
# Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md (decision 5, decision 7, 6.2); owner choice c of 2026-10-22.
#
# What it does (no jackhmmer is run again; the raw output holds the final iteration only):
#   1. python pr6_drift_check.py extract-ids   the jackhmmer outside fungal hits of the three queries, read from the COMMITTED
#                                              round 2 record (docs/reports/data/pr6_gap_uniref50_2026-10-22/)
#   2. seqkit grep -f ids.txt UniRef50         one streaming pass over the 28.8 GB FASTA pulls the sequences of those hits
#   3. python pr6_drift_check.py split         per query: the query, its hits, and the hits with residues SHUFFLED (seed 1)
#   4. phmmer -E 10 --domE 10 --domtblout      query vs its hits (real) and query vs the shuffled hits (null): seconds
#   5. python pr6_drift_check.py summarize     class of every hit and the training rule on SUPPORTED hits only
# The classifier thresholds are fixed in the docstring of pr6_drift_check.py BEFORE any result: PROPOSALS for the owner.
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update and this job. Heavy work is a SLURM
# job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $IN/kiss_sequences.faa (+ .sha256)   already staged at .../_workdir/kiss_cluster_split/in/ (no new staging)
#   $UNIPROT_FASTA                       /srv/projects/db/Uniprot/2025_03/uniref50.fasta (the file of round 2; its size is checked)
#   $PROJ_ROOT (a clone at CODE_COMMIT)  holds the round 2 record and the scripts (merge of the PR that holds this script and #114)
# Modules: seqkit/2.4.0 and hmmer/3.4 (both listed by `module avail` on HPCC 2026-10-22; `seqkit grep -f` behaviour on a
# 70-million-record file is not measured).
#
# Time: the one pass over 28.8 GB is the long step (grep over the same file took part of a 3 h 19 min job, but the share is not
# measured); the phmmer steps take seconds. The #SBATCH lines are defaults; override on the sbatch line if needed.
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch -p batch --time=4:00:00 -c 4 --mem=8G \
#     --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,UNIPROT_FASTA=/srv/projects/db/Uniprot/2025_03/uniref50.fasta \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_drift_check.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH -J pr6_drift_check
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${UNIPROT_FASTA:?export UNIPROT_FASTA (the UniRef50 FASTA of round 2; no path is assumed)}"
EXPECTED_UNIPROT_BYTES=28781193344   # round 2 provenance: uniref50.fasta 2025_03 (sha256 bc52bdd6c97ba2cc365b7e42ae7abb34865e0a1a82a2f95a74590aae040c81ff)
QUERIES="P32323 Q07730 A0ACE4KFX7"   # aga1, ece1, hyd1_like (the buildable_iterated_only queries under the fixed rule)
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_drift_check.$$"
trap 'rm -rf "$TMP" "$OUT".tmp.$$' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
WF=analysis/wall_fusion_truth
RECORD=docs/reports/data/pr6_gap_uniref50_2026-10-22

fail() { echo "FATAL: $*" >&2; exit 1; }

check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$RECORD" \
    || fail "$WF or $RECORD changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$RECORD")" ] || fail "uncommitted changes in $WF or $RECORD"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }
# one_seq FASTA ACCESSION: the sequence whose header's first token is exactly the accession (nothing if absent)
one_seq() { awk -v id="$2" '/^>/ { keep = ($1 == ">" id) } keep' "$1"; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -f "$PROJ_ROOT/$RECORD/gap_homologs.tsv.gz" ] || fail "$RECORD/gap_homologs.tsv.gz is missing (the clone does not hold the merged round 2 record)"
[ -f "$IN/kiss_sequences.faa" ] && [ -f "$IN/kiss_sequences.faa.sha256" ] || fail "$IN/kiss_sequences.faa or its .sha256 is missing"
sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa differs from its .sha256"
[ -s "$UNIPROT_FASTA" ] || fail "$UNIPROT_FASTA is missing or empty"
BYTES="$(stat -c %s "$UNIPROT_FASTA")"
[ "$BYTES" = "$EXPECTED_UNIPROT_BYTES" ] || fail "$UNIPROT_FASTA has $BYTES bytes, round 2 used $EXPECTED_UNIPROT_BYTES: it is not the same file"
mkdir -p "$TMP/ex" "$TMP/q" "$TMP/real" "$TMP/null"
PY="/usr/bin/python3.12 -I $PROJ_ROOT/$WF/pr6_drift_check.py"

# 1. The hits to check, from the committed record.
$PY extract-ids --record "$PROJ_ROOT/$RECORD" --out "$TMP/ex/ids_dir" --queries "$(echo "$QUERIES" | tr ' ' ',')"
N_IDS=$(wc -l < "$TMP/ex/ids_dir/ids.txt")
[ "$N_IDS" -gt 0 ] || fail "no ids to extract"

# 2. One streaming pass over UniRef50; every id must come back exactly once.
module load seqkit/2.4.0 hmmer/3.4
seqkit grep -j "$THREADS" -f "$TMP/ex/ids_dir/ids.txt" "$UNIPROT_FASTA" > "$TMP/hits_all.faa"
N_FA=$(grep -c '^>' "$TMP/hits_all.faa")
[ "$N_FA" = "$N_IDS" ] || fail "seqkit returned $N_FA sequences for $N_IDS ids (a hit is missing or duplicated)"
echo "extracted $N_FA hit sequences"

# 3. The queries, then the per-query real and shuffled hit sets.
for acc in $QUERIES; do
  one_seq "$IN/kiss_sequences.faa" "$acc" > "$TMP/q_$acc.faa"
  [ -s "$TMP/q_$acc.faa" ] || fail "$acc is not in kiss_sequences.faa"
done
cat $(for acc in $QUERIES; do echo "$TMP/q_$acc.faa"; done) > "$TMP/queries.faa"
$PY split --manifest "$TMP/ex/ids_dir/hits_manifest.tsv" --hits-faa "$TMP/hits_all.faa" --queries-faa "$TMP/queries.faa" \
  --out "$TMP/q" --seed 1

# 4. The pairwise searches: all pairs up to E 10 so that a weak or absent pair is visible; a missing trailer stops the job.
for acc in $QUERIES; do
  phmmer -E 10 --domE 10 --cpu "$THREADS" --noali -o /dev/null --domtblout "$TMP/real/$acc.dom" "$TMP/q/$acc.faa" "$TMP/q/$acc.hits.faa"
  phmmer -E 10 --domE 10 --cpu "$THREADS" --noali -o /dev/null --domtblout "$TMP/null/$acc.dom" "$TMP/q/$acc.faa" "$TMP/q/$acc.null.faa"
  for d in real null; do
    tail -n 3 "$TMP/$d/$acc.dom" | grep -qx '# \[ok\]' || fail "$d/$acc.dom has no '# [ok]' trailer"
  done
  echo "searched $acc"
done

# 5. The tables (standard library only; no sequence is written). The summary refuses a non-empty output directory.
$PY summarize --manifest "$TMP/ex/ids_dir/hits_manifest.tsv" --queries-faa "$TMP/queries.faa" --hits-faa "$TMP/hits_all.faa" \
  --real-dir "$TMP/real" --null-dir "$TMP/null" --out "$TMP/tables"

{
  echo "job: pr6_drift_check (PR 6, drift check of the iterated-only queries)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(phmmer -h | sed -n 2p)"
  echo "seqkit: $(seqkit version)"
  echo "python: $(/usr/bin/python3.12 --version)"
  echo "uniprot_fasta: $UNIPROT_FASTA"
  echo "uniprot_fasta_bytes: $BYTES (equals the round 2 file; sha256 of round 2: bc52bdd6c97ba2cc365b7e42ae7abb34865e0a1a82a2f95a74590aae040c81ff, not recomputed)"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "queries: $QUERIES"
  echo "hit_sequences_extracted: $N_FA"
  echo "settings: phmmer -E 10 --domE 10 (pairwise, real and shuffled); shuffle seed 1; thresholds in the docstring of pr6_drift_check.py"
} > "$TMP/tables/provenance.txt"
cp "$TMP/ex/ids_dir/hits_manifest.tsv" "$TMP/tables/hits_manifest.tsv"
cp -r "$TMP/tables" "$OUT.tmp.$$" && mv "$OUT.tmp.$$" "$OUT"
echo "wrote $OUT"
