#!/bin/bash -l
# PR 6 gap search (DRAFT, NOT SUBMITTED; nothing of PR 6 has run on HPCC after step 6.0): the outside-homolog search of
# steps 6.1 and 6.2 for the clusters that have NO Pfam hit (the four rows of custom_targets.tsv of the step 6.0 record:
# hyd1_like, aga1, ece1, ecm33), and the Ecm33 check. It builds no model and writes no call, no status and no tier.
# Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md, decisions 3, 5 and 7.
#
# 1. The Ecm33 check FIRST. Q5AGC4 (C. albicans Ecm33) has no GA hit although Ecm33 (PF12454) is in the panel, and report 07
#    row 22 found no hit of S. cerevisiae Ecm33 (P38248) even at E 10. The job runs `hmmsearch --max -E 10` of the model PF12454
#    against both sequences and the summary prints the score against the GA thresholds of the model (read from the model file)
#    and the length of each sequence against the 406 aa of P38248. P38248 is read from the committed shared negatives
#    (data/controls/hard-negatives-shared/sequences.faa), so nothing is staged for it.
# 2. The outside-homolog search. For each of the four queries (members of the clusters, from the staged kiss_sequences.faa):
#    - `phmmer` (one pass): the conservative set; it cannot drift.
#    - `jackhmmer -N 3`: iteration can reach more distant homologs and can also drift into other families; so the summary
#      calls a target buildable only when the one-pass set meets the rule, and `buildable_iterated_only` when only the
#      iterated set does (a flag for the owner).
#    Both use the inclusion threshold --incE 1e-5 and --incdomE 1e-5, written here before any result is read; --tblout holds
#    the hits of the last iteration. The target is the UniProt set in $UNIPROT_FASTA (a REQUIRED argument: no UniProt path is
#    documented in the repository, so none is guessed). The spec names the reviewed fungal proteins as the first choice; the
#    set, its download date and its sha256 are recorded in provenance.txt.
#
# Target set: Swiss-Prot (UniProt headers, `OX=`) or UniRef50 (headers with `TaxID=` and `RepID=`; every hit is a cluster
# representative); pr6_gap_summary.py reads both. The partition, the time and the cores are defaults on the #SBATCH lines
# below: a longer run (UniRef50 is 28.8 GB) overrides them on the sbatch command line (`-p`, `--time`, `-c`, `--mem`);
# nothing in the script fixes them (THREADS follows SLURM_CPUS_PER_TASK). E-values depend on the size of the target set.
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update, the UniProt staging, and this job.
# Heavy work is a SLURM job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $IN/kiss_sequences.faa (+ .sha256)   already staged at .../_workdir/kiss_cluster_split/in/ (no new staging for it)
#   $UNIPROT_FASTA                       the UniProt protein FASTA to search (staged by the owner; see README_pr6_gap.md)
#   $TAXDUMP                             nodes.dmp and names.dmp (default: the taxdump of the sorting-hat jobs, see below)
#   $PROJ_ROOT (a clone at CODE_COMMIT)  holds the rest, all committed (step 6.0 record, shared negatives)
# Pfam: PFAM_HMM (default the `current` link) must resolve to a directory named *Pfam$PFAM_RELEASE (only this stops the job);
# the sha256 is recorded and compared with the earlier jobs. Module: hmmer/3.4.
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2,UNIPROT_FASTA=/path/to/uniprot.fasta \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_gap_search.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist.)
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -J pr6_gap_search
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
: "${UNIPROT_FASTA:?export UNIPROT_FASTA (the UniProt protein FASTA to search; no path is assumed)}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
DATA_ROOT="${DATA_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
TAXDUMP="${TAXDUMP:-$DATA_ROOT/_workdir/sorting_hat/taxdump}"  # TAXDUMP of scripts/sorting_hat/refresh_3_10/lib.sh
EARLIER_PFAM_SHA256=4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e  # the earlier provenance (Pfam 38.2)
INC_E=1e-5     # inclusion and reporting E-value of the searches (fixed before the results are read)
JACK_ITER=3    # at most 3 iterations
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_gap_search.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-8}"
WF=analysis/wall_fusion_truth
SCAN=docs/reports/data/pr6_pfam_scan_2026-10-22
NEG=data/controls/hard-negatives-shared

fail() { echo "FATAL: $*" >&2; exit 1; }

check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$SCAN" "$NEG" \
    || fail "$WF, $SCAN or $NEG changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$SCAN" "$NEG")" ] || fail "uncommitted changes in $WF, $SCAN or $NEG"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }
# one_seq FASTA ACCESSION: the sequence whose header's first token is exactly the accession (nothing if absent)
one_seq() { awk -v id="$2" '/^>/ { keep = ($1 == ">" id) } keep' "$1"; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
TARGETS="$PROJ_ROOT/$SCAN/custom_targets.tsv"
[ -f "$TARGETS" ] || fail "$TARGETS is missing (the clone does not hold the merged step 6.0 record)"
[ -f "$IN/kiss_sequences.faa" ] && [ -f "$IN/kiss_sequences.faa.sha256" ] || fail "$IN/kiss_sequences.faa or its .sha256 is missing"
sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa differs from its .sha256"
[ -s "$UNIPROT_FASTA" ] || fail "$UNIPROT_FASTA is missing or empty"
[ -f "$TAXDUMP/nodes.dmp" ] && [ -f "$TAXDUMP/names.dmp" ] || fail "$TAXDUMP/nodes.dmp or names.dmp is missing"
mkdir -p "$TMP/q" "$TMP/phmmer" "$TMP/jackhmmer"

# The queries: the accessions of the targets (the table lists members comma-separated in `accessions`).
tail -n +2 "$TARGETS" | cut -f7 | tr ',' '\n' | sort -u > "$TMP/queries.txt"
N_Q=$(wc -l < "$TMP/queries.txt")
[ "$N_Q" -eq 4 ] || fail "custom_targets.tsv lists $N_Q query accessions, not 4"
while read -r acc; do
  one_seq "$IN/kiss_sequences.faa" "$acc" > "$TMP/q/$acc.faa"
  [ -s "$TMP/q/$acc.faa" ] || fail "$acc is not in kiss_sequences.faa"
done < "$TMP/queries.txt"
# Ecm33 check inputs: Q5AGC4 (a target) and P38248 (a shared negative, committed).
one_seq "$PROJ_ROOT/$NEG/sequences.faa" P38248 > "$TMP/p38248.faa"
[ -s "$TMP/p38248.faa" ] || fail "P38248 is not in $NEG/sequences.faa"
cat "$TMP/q/Q5AGC4.faa" "$TMP/p38248.faa" > "$TMP/ecm33.faa"
awk '/^>/ { if (n) print n_id "\t" n; n_id = substr($1, 2); n = 0; next } { n += length($0) } END { print n_id "\t" n }' \
  "$TMP/ecm33.faa" > "$TMP/ecm33_lengths.tsv"
# One pass over the target set for the count (UniRef50 is 28.8 GB: do not read it twice for the count).
N_TARGET_SEQ="$(grep -c '^>' "$UNIPROT_FASTA")"
echo "inputs ok: 4 queries, target set $N_TARGET_SEQ sequences"

module load hmmer/3.4
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
HMM_SHA="$(sha256sum "$RESOLVED" | cut -d' ' -f1)"
if [ "$HMM_SHA" = "$EARLIER_PFAM_SHA256" ]; then PFAM_SAME=yes; else
  PFAM_SAME=no
  echo "REPORT: sha256 of $RESOLVED is $HMM_SHA, not the $EARLIER_PFAM_SHA256 of the earlier jobs; recorded, the run goes on."
fi

# 1. The Ecm33 check. hmmfetch -f looks up by the model NAME in the pressed database, not by the accession.
echo Ecm33 > "$TMP/keys.txt"
hmmfetch -f "$PFAM_HMM" "$TMP/keys.txt" > "$TMP/Ecm33.hmm"
grep -q '^ACC   PF12454' "$TMP/Ecm33.hmm" || fail "the fetched model Ecm33 is not PF12454"
hmmsearch --max -E 10 --domE 10 --cpu "$THREADS" --noali -o /dev/null --domtblout "$TMP/ecm33.domtbl" "$TMP/Ecm33.hmm" "$TMP/ecm33.faa"
tail -n 3 "$TMP/ecm33.domtbl" | grep -qx '# \[ok\]' || fail "ecm33.domtbl has no '# [ok]' trailer"

# 2. The searches, one query at a time (a missing trailer stops the job).
while read -r acc; do
  phmmer -E "$INC_E" --domE "$INC_E" --cpu "$THREADS" --noali -o /dev/null --tblout "$TMP/phmmer/$acc.tbl" \
    "$TMP/q/$acc.faa" "$UNIPROT_FASTA"
  jackhmmer -N "$JACK_ITER" --incE "$INC_E" --incdomE "$INC_E" -E "$INC_E" --domE "$INC_E" --cpu "$THREADS" --noali \
    -o /dev/null --tblout "$TMP/jackhmmer/$acc.tbl" "$TMP/q/$acc.faa" "$UNIPROT_FASTA"
  for m in phmmer jackhmmer; do
    tail -n 3 "$TMP/$m/$acc.tbl" | grep -qx '# \[ok\]' || fail "$m/$acc.tbl has no '# [ok]' trailer"
  done
  echo "searched $acc"
done < "$TMP/queries.txt"

# 3. The tables (standard library only). The summary script refuses a non-empty output directory.
#    If it fails after hours of searching, the raw tables are kept next to $OUT (the trap deletes $TMP) before the job stops.
if ! /usr/bin/python3.12 -I "$PROJ_ROOT/$WF/pr6_gap_summary.py" \
  --targets "$TARGETS" --phmmer-dir "$TMP/phmmer" --jackhmmer-dir "$TMP/jackhmmer" \
  --nodes "$TAXDUMP/nodes.dmp" --names "$TAXDUMP/names.dmp" \
  --ecm33-domtbl "$TMP/ecm33.domtbl" --ecm33-model "$TMP/Ecm33.hmm" --ecm33-lengths "$TMP/ecm33_lengths.tsv" \
  --out "$TMP/tables"; then
  SALVAGE="$OUT.raw_after_summary_failure"
  if [ ! -e "$SALVAGE" ]; then
    mkdir -p "$SALVAGE" && cp -r "$TMP/phmmer" "$TMP/jackhmmer" "$TMP/ecm33.domtbl" "$SALVAGE/" && echo "raw tables kept in $SALVAGE"
  fi
  fail "pr6_gap_summary.py failed"
fi

{
  echo "job: pr6_gap_search (PR 6 steps 6.1 and 6.2, the four clusters with no Pfam hit)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(phmmer -h | sed -n 2p)"
  echo "python: $(/usr/bin/python3.12 --version)"
  echo "pfam_hmm: $PFAM_HMM"
  echo "pfam_resolved: $RESOLVED"
  echo "pfam_release: $PFAM_RELEASE"
  echo "pfam_sha256: $HMM_SHA"
  echo "pfam_sha256_same_as_earlier_jobs: $PFAM_SAME (earlier: $EARLIER_PFAM_SHA256)"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "uniprot_fasta: $UNIPROT_FASTA"
  echo "uniprot_fasta_bytes: $(stat -c %s "$UNIPROT_FASTA")"
  echo "uniprot_fasta_sequences: $N_TARGET_SEQ"
  echo "uniprot_fasta_sha256: $(sha256sum < "$UNIPROT_FASTA" | cut -d' ' -f1)"
  echo "taxdump: $TAXDUMP"
  echo "settings: phmmer -E $INC_E; jackhmmer -N $JACK_ITER --incE $INC_E --incdomE $INC_E; Ecm33: hmmsearch --max -E 10"
} > "$TMP/provenance.txt"

# Write the output last, in a directory that did not exist. The query sequences stay in $TMP (they are not an output).
[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT/raw/phmmer" "$OUT/raw/jackhmmer"
cp "$TMP"/tables/* "$OUT/"
cp "$TMP/provenance.txt" "$TMP/ecm33.domtbl" "$OUT/"
cp "$TMP"/phmmer/*.tbl "$OUT/raw/phmmer/"
cp "$TMP"/jackhmmer/*.tbl "$OUT/raw/jackhmmer/"
echo "wrote $OUT"
