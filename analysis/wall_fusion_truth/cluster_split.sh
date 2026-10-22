#!/bin/bash -l
# Task 5.3 (DRAFT, NOT SUBMITTED; nothing of this plan has run on HPCC): cluster the Kiss control sequences
# (and the hard negatives, only if the owner staged them) with MMseqs2 at 30% identity and 50% coverage, search
# them against the Pfam seed sequences of the wall, fusion and adhesion-family models, and label every cluster
# `seed-like` or `seed-independent` (analysis/wall_fusion_truth/seed_split.py). Spec D8, D20; plan Task 5.3.
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone, the staging of the inputs,
# and this job. It reads inputs only from $IN and writes only $OUT (an existing $OUT stops the job).
# Heavy work is a SLURM job, never on the login node. $SCRATCH is node-local; there is no /tmp and no BASH_SOURCE.
#
# Inputs in $IN (made on the owner's machine by build_kiss_fasta.py and fetch_seeds_models.py; each has a
# `.sha256` file next to it that this job checks):
#   kiss_sequences.faa   the sequences to cluster (header = accession)
#   members.tsv          accession, class, slice, ...  (one line per sequence)
#   seeds.faa.gz         Pfam seed sequences, header `<PFxxxxx>|<name>`
# Modules (from the repository, not from HPCC: could not be confirmed for today's HPCC): `MMseqs2/17-b804f`
# (docs/agent-tasks/COMMON-RULES.md rule 4; analysis/lectin_validation/run_cluster_seed.sh;
# analysis/calibration_truth/repeat_call_truth/build_truth.py). `hmmer/3.4` is NOT needed here (it is for Task 5.4).
# Python: /usr/bin/python3.12, standard library only (no virtual environment, no wheelhouse).
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds these scripts):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/inputs,OUT=/path/to/new_outdir \
#     -o $OUT.log -e $OUT.log analysis/wall_fusion_truth/cluster_split.sh
#   (the log directory must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J kiss_cluster_split
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the staged input directory) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/kiss_cluster_split.$$"
trap 'rm -rf "$TMP" "$OUT".tmp.$$' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
SET=data/controls/wall-fusion-kiss2019
WF=analysis/wall_fusion_truth

fail() { echo "FATAL: $*" >&2; exit 1; }

# The code and the control set must be the reviewed ones.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$SET" \
    || fail "$WF or $SET changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$SET")" ] || fail "uncommitted changes in $WF or $SET"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

# sidecar_ok FILE: the .sha256 next to FILE holds the sha256 of FILE
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
for f in kiss_sequences.faa members.tsv seeds.faa.gz; do
  [ -f "$IN/$f" ] || fail "$IN/$f is missing"
  [ -f "$IN/$f.sha256" ] || fail "$IN/$f.sha256 is missing"
  sidecar_ok "$IN/$f" || fail "$IN/$f differs from its .sha256"
done

# Every sequence of the FASTA is a member, and every non-negative member is in the committed control set.
n_fa=$(grep -c '^>' "$IN/kiss_sequences.faa")
n_mem=$(($(wc -l < "$IN/members.tsv") - 1))
[ "$n_fa" = "$n_mem" ] || fail "kiss_sequences.faa has $n_fa sequences but members.tsv has $n_mem rows"
n_pos=$(awk -F'\t' 'NR>1 && $2=="positive"' "$IN/members.tsv" | wc -l)
n_ctl=$(($(wc -l < "$PROJ_ROOT/$SET/controls.tsv") - 1 + $(wc -l < "$PROJ_ROOT/$SET/controls_local_sequence.tsv") - 1))
[ "$n_pos" = "$n_ctl" ] || fail "members.tsv has $n_pos positives but the committed control set has $n_ctl rows"
echo "inputs ok: $n_fa sequences, $n_pos positives"

mkdir -p "$TMP"
module load MMseqs2/17-b804f
mmseqs version > "$TMP/mmseqs_version.txt"
zcat "$IN/seeds.faa.gz" > "$TMP/seeds.faa"

# 1. Clusters: 30% identity, 50% coverage (COMMON-RULES rule 4; the call of run_cluster_seed.sh).
mmseqs easy-cluster "$IN/kiss_sequences.faa" "$TMP/clu" "$TMP/clu_tmp" --min-seq-id 0.3 -c 0.5 --threads "$THREADS" -v 1
# 2. Identity of every sequence to the nearest Pfam seed sequence (the leakage search of finalize_truth.py).
#    query = control sequence, target = seed sequence; identity at 0.3, target coverage 0.5, E <= 1e-3.
mmseqs easy-search "$IN/kiss_sequences.faa" "$TMP/seeds.faa" "$TMP/seed_hits.m8" "$TMP/search_tmp" \
  --min-seq-id 0.3 --cov-mode 1 -c 0.5 -s 7.5 -e 1e-3 --threads "$THREADS" -v 1 \
  --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"
# 3. The labels and the tables (standard library only).
/usr/bin/python3.12 -I "$PROJ_ROOT/$WF/seed_split.py" --members "$IN/members.tsv" --clusters "$TMP/clu_cluster.tsv" \
  --seed-hits "$TMP/seed_hits.m8" --controls "$PROJ_ROOT/$SET" --out "$TMP/split"

{
  echo "job: kiss_cluster_split (plan Task 5.3)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "mmseqs: $(cat "$TMP/mmseqs_version.txt")"
  echo "python: $(/usr/bin/python3.12 --version)"
  for f in kiss_sequences.faa members.tsv seeds.faa.gz; do echo "input $f sha256 $(sha256sum < "$IN/$f" | cut -d' ' -f1)"; done
  echo "settings: easy-cluster --min-seq-id 0.3 -c 0.5; easy-search --min-seq-id 0.3 --cov-mode 1 -c 0.5 -s 7.5 -e 1e-3"
} > "$TMP/provenance.txt"

# Copy under a temporary name, then rename; an existing $OUT is never replaced.
mkdir "$OUT.tmp.$$"
cp -a "$TMP/split/." "$OUT.tmp.$$/"
cp -a "$TMP/clu_cluster.tsv" "$TMP/seed_hits.m8" "$TMP/mmseqs_version.txt" "$TMP/provenance.txt" "$OUT.tmp.$$/"
mv "$OUT.tmp.$$" "$OUT"
echo "wrote $OUT"
