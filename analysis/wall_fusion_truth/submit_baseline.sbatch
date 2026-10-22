#!/bin/bash -l
# Task 5.4 (DRAFT, NOT SUBMITTED; nothing of this plan has run on HPCC): search the 39 Pfam models (the 21 wall and fusion
# models of analysis/wall_fusion_report/wall_fusion_models.tsv and the 18 active family rows of data/sorting_hat/family_table.tsv)
# with `hmmsearch --cut_ga` against the Kiss control sequences and the shared negatives, then write the baseline recall tables
# (analysis/wall_fusion_truth/baseline_recall.py). It writes no call, no status and no tier. Spec D8, D9 step 1, D20; plan Task 5.4.
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone, the staging of the input, and this job.
# Heavy work is a SLURM job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $IN/kiss_sequences.faa (+ .sha256)   the same file as in the cluster-split job (293 sequences: 132 controls, 161 negatives);
#                                        already staged at .../_workdir/kiss_cluster_split/in/ for job 30229366
#   $PROJ_ROOT (a clone at CODE_COMMIT) holds the rest, all committed: the cluster tables of the run record
#                                        docs/reports/data/kiss_cluster_split_2026-10-22/{members_split,clusters_split}.tsv,
#                                        the control set data/controls/wall-fusion-kiss2019/ (both control files and
#                                        family_rep_sequences.faa for the 35 representatives) and the two model lists.
# Pfam: PFAM_HMM (default the `current` link) must resolve to a directory named *Pfam$PFAM_RELEASE (the check of
# scripts/sorting_hat/pfam_hmmsearch.sbatch). The sha256 of the resolved Pfam-A.hmm is recorded in provenance.json and compared
# with the one of the earlier jobs; a DIFFERENT sha256 is recorded and reported (a `REPORT:` line), it does not stop the job
# (only a release-name mismatch does). Module: hmmer/3.4 (the one of pfam_hmmsearch.sbatch).
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o $OUT.log -e $OUT.log analysis/wall_fusion_truth/submit_baseline.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J kiss_baseline
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
EARLIER_PFAM_SHA256=4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e  # the earlier provenance (Pfam 38.2)
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/kiss_baseline.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
SET=data/controls/wall-fusion-kiss2019
WF=analysis/wall_fusion_truth
REC=docs/reports/data/kiss_cluster_split_2026-10-22
MODELS=analysis/wall_fusion_report/wall_fusion_models.tsv
FAMILY=data/sorting_hat/family_table.tsv

fail() { echo "FATAL: $*" >&2; exit 1; }

# The code, the control set, the cluster record and the model lists must be the reviewed ones.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$SET" "$REC" "$MODELS" "$FAMILY" \
    || fail "$WF, $SET, $REC or a model list changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$SET" "$REC" "$MODELS" "$FAMILY")" ] \
    || fail "uncommitted changes in $WF, $SET, $REC or a model list"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

# sidecar_ok FILE: the .sha256 next to FILE holds the sha256 of FILE
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -f "$IN/kiss_sequences.faa" ] && [ -f "$IN/kiss_sequences.faa.sha256" ] || fail "$IN/kiss_sequences.faa or its .sha256 is missing"
sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa differs from its .sha256"
# The sequences are the ones the cluster job used: every one is a row of the committed members_split.tsv, and the reverse.
n_fa=$(grep -c '^>' "$IN/kiss_sequences.faa")
n_mem=$(($(wc -l < "$PROJ_ROOT/$REC/members_split.tsv") - 1))
[ "$n_fa" = "$n_mem" ] || fail "kiss_sequences.faa has $n_fa sequences but members_split.tsv has $n_mem rows"
mkdir -p "$TMP"
grep '^>' "$IN/kiss_sequences.faa" | sed 's/^>//; s/[[:space:]].*//' | sort > "$TMP/fa_ids.txt"
tail -n +2 "$PROJ_ROOT/$REC/members_split.tsv" | cut -f1 | sort > "$TMP/mem_ids.txt"
cmp -s "$TMP/fa_ids.txt" "$TMP/mem_ids.txt" || fail "the sequence names of kiss_sequences.faa differ from the accessions of members_split.tsv"
echo "inputs ok: $n_fa sequences"

module load hmmer/3.4

# Pfam: the resolved path must be a directory named *Pfam$PFAM_RELEASE (only this stops the job).
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
HMM_SHA="$(sha256sum "$RESOLVED" | cut -d' ' -f1)"
if [ "$HMM_SHA" = "$EARLIER_PFAM_SHA256" ]; then PFAM_SAME=yes; else
  PFAM_SAME=no
  echo "REPORT: sha256 of $RESOLVED is $HMM_SHA, not the $EARLIER_PFAM_SHA256 of the earlier jobs; recorded, the run goes on (release name $PFAM_RELEASE matches)."
fi

# The 39 models: ids and names from the two repository tables (standard library only).
/usr/bin/python3.12 -I - "$PROJ_ROOT/$MODELS" "$PROJ_ROOT/$FAMILY" "$TMP" <<'PY'
import csv, sys
wall, family, tmp = sys.argv[1:4]
rows = []
with open(wall, encoding="utf-8-sig", newline="") as fh:
    rows += [(r["id"], r["name"]) for r in csv.DictReader(fh, delimiter="\t")]
with open(family, encoding="utf-8-sig", newline="") as fh:
    rows += [(r["pfam_acc"], r["name"]) for r in csv.DictReader(fh, delimiter="\t") if r["active"] == "yes" and r["pfam_acc"]]
rows = sorted(set(rows))
accs = [a for a, _ in rows]
assert len(accs) == len(set(accs)), "an accession has two names"
open(f"{tmp}/keys.txt", "w").write("".join(n + "\n" for _, n in rows))
open(f"{tmp}/expected.tsv", "w").write("".join(f"{a}\t{n}\n" for a, n in rows))
PY
N_EXP=$(wc -l < "$TMP/expected.tsv")
[ "$N_EXP" -eq 39 ] || fail "the model lists give $N_EXP models, not 39"
# hmmfetch -f looks up by the model NAME in the pressed database, not by the accession
hmmfetch -f "$PFAM_HMM" "$TMP/keys.txt" > "$TMP/models.hmm"
N=$(grep -c '^NAME' "$TMP/models.hmm" || true)
[ "$N" -eq "$N_EXP" ] || fail "hmmfetch found $N of $N_EXP models"
awk '/^NAME/{n=$2} /^ACC/{split($2,a,"."); print a[1]"\t"n}' "$TMP/models.hmm" | sort > "$TMP/fetched.tsv"
sort "$TMP/expected.tsv" > "$TMP/expected.sorted.tsv"
diff -q "$TMP/fetched.tsv" "$TMP/expected.sorted.tsv" > /dev/null || { diff "$TMP/fetched.tsv" "$TMP/expected.sorted.tsv" >&2; fail "model names in the database differ from the model lists"; }

# A model without a GA line is REPORTED, not skipped: hmmsearch --cut_ga stops on such a model, so the models are split.
# The split models get an E-value cut-off, a separate table and a separate note; they enter no recall number.
awk -v ga="$TMP/with_ga.hmm" -v noga="$TMP/no_ga.hmm" -v list="$TMP/no_ga.txt" '
  { blk = blk $0 "\n" }
  /^ACC/ { split($2, a, "."); acc = a[1] }
  /^GA / { has = 1 }
  /^\/\// { if (has) printf "%s", blk > ga; else { printf "%s", blk > noga; print acc > list } blk = ""; has = 0; acc = "" }
' "$TMP/models.hmm"
touch "$TMP/with_ga.hmm" "$TMP/no_ga.hmm" "$TMP/no_ga.txt"
N_GA=$(grep -c '^NAME' "$TMP/with_ga.hmm" || true)
N_NOGA=$(grep -c '^NAME' "$TMP/no_ga.hmm" || true)
echo "models: $N_EXP fetched, $N_GA with a GA line, $N_NOGA without"
[ "$N_NOGA" -eq 0 ] || echo "REPORT: $N_NOGA model(s) have no GA line: $(tr '\n' ' ' < "$TMP/no_ga.txt")"

FAMREP="$PROJ_ROOT/$SET/family_rep_sequences.faa"
search_ga() { # search_ga FASTA DOMTBL
  hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --domtblout "$2" "$TMP/with_ga.hmm" "$1"
  tail -n 3 "$2" | grep -qx '# \[ok\]' || fail "$2 has no '# [ok]' trailer"
}
search_noga() { # search_noga FASTA DOMTBL
  hmmsearch -E 1e-5 --domE 1e-5 --cpu "$THREADS" --noali -o /dev/null --domtblout "$2" "$TMP/no_ga.hmm" "$1"
  tail -n 3 "$2" | grep -qx '# \[ok\]' || fail "$2 has no '# [ok]' trailer"
}
[ "$N_GA" -gt 0 ] || fail "no model has a GA line; nothing can be searched with --cut_ga"
search_ga "$IN/kiss_sequences.faa" "$TMP/domtbl_kiss.txt"
search_ga "$FAMREP" "$TMP/domtbl_famrep.txt"
EXTRA=()
if [ "$N_NOGA" -gt 0 ]; then
  search_noga "$IN/kiss_sequences.faa" "$TMP/domtbl_noga.txt"
  search_noga "$FAMREP" "$TMP/domtbl_noga_famrep.txt"
  EXTRA=(--models-no-ga "$TMP/no_ga.txt" --domtbl-noga "$TMP/domtbl_noga.txt" --domtbl-noga-famrep "$TMP/domtbl_noga_famrep.txt")
fi

# The tables (standard library only). baseline_recall.py refuses a non-empty output directory.
/usr/bin/python3.12 -I "$PROJ_ROOT/$WF/baseline_recall.py" \
  --domtbl "$TMP/domtbl_kiss.txt" --domtbl-famrep "$TMP/domtbl_famrep.txt" ${EXTRA[@]+"${EXTRA[@]}"} \
  --members "$PROJ_ROOT/$REC/members_split.tsv" --clusters "$PROJ_ROOT/$REC/clusters_split.tsv" \
  --controls "$PROJ_ROOT/$SET" --wall-models "$PROJ_ROOT/$MODELS" --family-table "$PROJ_ROOT/$FAMILY" \
  --out "$TMP/tables"

{
  echo "job: kiss_baseline (plan Task 5.4)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(hmmsearch -h | sed -n 2p)"
  echo "python: $(/usr/bin/python3.12 --version)"
  echo "pfam_hmm: $PFAM_HMM"
  echo "pfam_resolved: $RESOLVED"
  echo "pfam_release: $PFAM_RELEASE"
  echo "pfam_sha256: $HMM_SHA"
  echo "pfam_sha256_same_as_earlier_jobs: $PFAM_SAME (earlier: $EARLIER_PFAM_SHA256)"
  echo "models: $N_EXP (with GA $N_GA, without GA $N_NOGA)"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "input family_rep_sequences.faa sha256 $(sha256sum < "$FAMREP" | cut -d' ' -f1)"
  echo "settings: hmmsearch --cut_ga (models with a GA line); -E 1e-5 --domE 1e-5 only for a model without GA"
} > "$TMP/provenance.txt"
cp "$TMP/fetched.tsv" "$TMP/models_searched.tsv"

# Write the output last: tables, the raw domain tables and the provenance, all in a directory that did not exist.
[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP"/tables/* "$OUT/"
cp "$TMP"/domtbl_*.txt "$TMP/provenance.txt" "$TMP/models_searched.tsv" "$OUT/"
echo "wrote $OUT"
