#!/bin/bash -l
# PR 6 step 6.0 (DRAFT, NOT SUBMITTED; nothing of PR 6 has run on HPCC): scan the members of the 44 clusters that the 39-model
# panel misses against ALL of Pfam-A 38.2 with `hmmscan --cut_ga`, and write the per-cluster tables of
# analysis/wall_fusion_truth/pr6_pfam_scan_summary.py. It writes no call, no status and no tier.
# Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md, decision 1 and section 6 (step 6.0).
#
# Why hmmscan and not hmmsearch: the query set is 66 sequences and the target set is the whole of Pfam-A (about 20,000
# models). hmmscan is built for that shape: it reads the pressed binary database (Pfam-A.hmm.h3f, .h3i, .h3m, .h3p, which
# exist next to Pfam-A.hmm) and reports per sequence which models hit. hmmsearch would read all models as text for
# every run. `--cut_ga` is valid for hmmscan: a hit passes the GA threshold of the model it hits, so no model is
# skipped (a model without a GA line reports nothing under --cut_ga, which is the same rule as the Task 5.4 job).
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update and this job. Heavy work is a
# SLURM job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $IN/kiss_sequences.faa (+ .sha256)   the same file as in the cluster-split and baseline jobs (293 sequences);
#                                        already staged at .../_workdir/kiss_cluster_split/in/. No new staging.
#   $PROJ_ROOT (a clone at CODE_COMMIT)  holds the rest, all committed: the Task 5.4 record
#                                        docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv (the 44 clusters and
#                                        their member accessions), the cluster record, the groups table of the spec and
#                                        the two model lists. The clone must contain the merge of the PRs that hold these.
# Pfam: PFAM_HMM (default the `current` link) must resolve to a directory named *Pfam$PFAM_RELEASE; sha256 as in the
# Task 5.4 job (a different sha256 is recorded and reported, a release-name mismatch stops the job). Module: hmmer/3.4.
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_pfam_scan.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J pr6_pfam_scan
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
EARLIER_PFAM_SHA256=4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e  # the earlier provenance (Pfam 38.2)
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_pfam_scan.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
WF=analysis/wall_fusion_truth
BASE=docs/reports/data/kiss_baseline_2026-10-22
REC=docs/reports/data/kiss_cluster_split_2026-10-22
GROUPS_DIR=docs/superpowers/specs/2026-10-22-pr6-custom-hmms-data
MODELS=analysis/wall_fusion_report/wall_fusion_models.tsv
FAMILY=data/sorting_hat/family_table.tsv
N_MISSED_MEMBERS=66  # the members of the 44 clusters of clusters_missed.tsv (recomputed from the committed table)

fail() { echo "FATAL: $*" >&2; exit 1; }

# The code, the record tables and the model lists must be the reviewed ones.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$BASE" "$REC" "$GROUPS_DIR" "$MODELS" "$FAMILY" \
    || fail "$WF, a record table or a model list changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$BASE" "$REC" "$GROUPS_DIR" "$MODELS" "$FAMILY")" ] \
    || fail "uncommitted changes in $WF, a record table or a model list"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

# sidecar_ok FILE: the .sha256 next to FILE holds the sha256 of FILE
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
for f in "$PROJ_ROOT/$BASE/clusters_missed.tsv" "$PROJ_ROOT/$GROUPS_DIR/candidate_groups.tsv" "$PROJ_ROOT/$REC/members_split.tsv"; do
  [ -f "$f" ] || fail "$f is missing (the clone does not hold the merged record PRs)"
done
[ -f "$IN/kiss_sequences.faa" ] && [ -f "$IN/kiss_sequences.faa.sha256" ] || fail "$IN/kiss_sequences.faa or its .sha256 is missing"
sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa differs from its .sha256"
mkdir -p "$TMP"
# The sequences are the ones the cluster job used: every name is an accession of the committed members_split.tsv, and the reverse.
grep '^>' "$IN/kiss_sequences.faa" | sed 's/^>//; s/[[:space:]].*//' | sort > "$TMP/fa_ids.txt"
tail -n +2 "$PROJ_ROOT/$REC/members_split.tsv" | cut -f1 | sort > "$TMP/mem_ids.txt"
cmp -s "$TMP/fa_ids.txt" "$TMP/mem_ids.txt" || fail "the sequence names of kiss_sequences.faa differ from the accessions of members_split.tsv"

# The members of the 44 missed clusters (the table lists them comma-separated in the column `accessions`).
/usr/bin/python3.12 -I - "$PROJ_ROOT/$BASE/clusters_missed.tsv" "$TMP/missed_ids.txt" <<'PY'
import csv, sys
src, dst = sys.argv[1:3]
with open(src, encoding="utf-8-sig", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
ids = sorted({a for r in rows for a in r["accessions"].split(",")})
open(dst, "w").write("".join(i + "\n" for i in ids))
print(f"{len(rows)} missed clusters, {len(ids)} member sequences")
PY
N_IDS=$(wc -l < "$TMP/missed_ids.txt")
[ "$N_IDS" -eq "$N_MISSED_MEMBERS" ] || fail "clusters_missed.tsv lists $N_IDS members, not $N_MISSED_MEMBERS"
# Exact match on the first token of the header (a `#neg` copy is a different name and is never taken).
awk 'NR==FNR { want[$1] = 1; next }
     /^>/ { name = substr($1, 2); keep = (name in want); if (keep) seen[name] = 1 }
     keep { print }
     END { for (n in want) if (!(n in seen)) { print "MISSING " n > "/dev/stderr"; bad = 1 } exit bad }' \
  "$TMP/missed_ids.txt" "$IN/kiss_sequences.faa" > "$TMP/missed.faa" || fail "a member of a missed cluster is not in kiss_sequences.faa"
N_FA=$(grep -c '^>' "$TMP/missed.faa")
[ "$N_FA" -eq "$N_IDS" ] || fail "extracted $N_FA sequences for $N_IDS members"
echo "inputs ok: $N_FA sequences of the missed clusters"

module load hmmer/3.4

# Pfam: the resolved path must be a directory named *Pfam$PFAM_RELEASE (only this stops the job).
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
for ext in h3f h3i h3m h3p; do [ -f "$RESOLVED.$ext" ] || fail "$RESOLVED.$ext is missing: hmmscan needs the pressed database"; done
HMM_SHA="$(sha256sum "$RESOLVED" | cut -d' ' -f1)"
if [ "$HMM_SHA" = "$EARLIER_PFAM_SHA256" ]; then PFAM_SAME=yes; else
  PFAM_SAME=no
  echo "REPORT: sha256 of $RESOLVED is $HMM_SHA, not the $EARLIER_PFAM_SHA256 of the earlier jobs; recorded, the run goes on (release name $PFAM_RELEASE matches)."
fi

# 1. The GA scan: a row is a hit that passed the GA threshold of the model it hits.
hmmscan --cut_ga --cpu "$THREADS" --noali -o /dev/null --domtblout "$TMP/domtbl_ga.txt" "$PFAM_HMM" "$TMP/missed.faa"
tail -n 3 "$TMP/domtbl_ga.txt" | grep -qx '# \[ok\]' || fail "domtbl_ga.txt has no '# [ok]' trailer"
# 2. The weak-hit table (descriptive, clearly separate): E <= 1e-3, no GA cut-off. It changes no class.
hmmscan -E 1e-3 --domE 1e-3 --cpu "$THREADS" --noali -o /dev/null --domtblout "$TMP/domtbl_weak.txt" "$PFAM_HMM" "$TMP/missed.faa"
tail -n 3 "$TMP/domtbl_weak.txt" | grep -qx '# \[ok\]' || fail "domtbl_weak.txt has no '# [ok]' trailer"

# 3. The tables (standard library only). The summary script refuses a non-empty output directory.
/usr/bin/python3.12 -I "$PROJ_ROOT/$WF/pr6_pfam_scan_summary.py" \
  --domtbl-ga "$TMP/domtbl_ga.txt" --domtbl-weak "$TMP/domtbl_weak.txt" \
  --missed "$PROJ_ROOT/$BASE/clusters_missed.tsv" --groups "$PROJ_ROOT/$GROUPS_DIR/candidate_groups.tsv" \
  --wall-models "$PROJ_ROOT/$MODELS" --family-table "$PROJ_ROOT/$FAMILY" --out "$TMP/tables"

{
  echo "job: pr6_pfam_scan (PR 6 step 6.0)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(hmmscan -h | sed -n 2p)"
  echo "python: $(/usr/bin/python3.12 --version)"
  echo "pfam_hmm: $PFAM_HMM"
  echo "pfam_resolved: $RESOLVED"
  echo "pfam_release: $PFAM_RELEASE"
  echo "pfam_sha256: $HMM_SHA"
  echo "pfam_sha256_same_as_earlier_jobs: $PFAM_SAME (earlier: $EARLIER_PFAM_SHA256)"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "queries: $N_FA sequences of the missed clusters (sha256 of the extract $(sha256sum < "$TMP/missed.faa" | cut -d' ' -f1))"
  echo "settings: hmmscan --cut_ga (hits); hmmscan -E 1e-3 --domE 1e-3 without --cut_ga (weak hits, descriptive)"
} > "$TMP/provenance.txt"

# Write the output last, in a directory that did not exist. The extracted sequences stay in $TMP (they are not an output).
[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP"/tables/* "$OUT/"
cp "$TMP"/domtbl_ga.txt "$TMP"/domtbl_weak.txt "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
