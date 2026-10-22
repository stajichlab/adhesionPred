#!/bin/bash -l
# PR 6 family specificity tests (DRAFT, NOT SUBMITTED; nothing of PR 6 beyond step 6.0 has run on HPCC): search the Pfam families
# that the step 6.0 scan names as candidates (union rule of the owner decision 1 of 2026-10-22: a hit covering >= 0.5 of the
# sequence OR a best hit with an E-value <= 1e-30; 26 clusters, 30 families) with `hmmsearch --cut_ga` over (i) the Kiss
# positives, (ii) the shared hard-negative set of 161 rows, (iii) the four kept proteomes (one per species, decision 9b), and
# write the report-only tables of analysis/wall_fusion_truth/pr6_family_specificity.py. It writes no call, no status and no tier,
# and it does not touch data/sorting_hat/family_measure.tsv or family_table.tsv (owner decision 2: nothing enters the active
# panel without a family_measure row AND the owner's sign-off; narrow_ok stays no).
#
# Needs the owner's go-ahead for EACH step before it is submitted (D21): the clone update and this job. Heavy work is a SLURM
# job, never on the login node. $SCRATCH is the temporary directory; there is no /tmp and no BASH_SOURCE.
#
# Inputs:
#   $IN/kiss_sequences.faa (+ .sha256)  the staged FASTA of the cluster-split job (293 sequences); only its 132 positives are used.
#                                       Already staged at .../_workdir/kiss_cluster_split/in/. No new staging.
#   $PROJ_ROOT (a clone at CODE_COMMIT) holds the rest, all committed: the step 6.0 record (scan_per_cluster.tsv), clusters_missed.tsv,
#                                       members_split.tsv, the shared negatives (data/controls/hard-negatives-shared/, FASTA, controls,
#                                       clusters), the committed candidate tables and scripts/sorting_hat/refresh_3_10/inputs.tsv.
#                                       The clone must contain the merges of the PRs that hold these.
#   $PROTEOME_DIR/<name>.faa            the four kept proteomes: Afum_Af293_UniProt, Calb_SC5314, Cimm_RS, Scer_S288C (read only). The
#                                       default is the `OLD` directory of scripts/sorting_hat/refresh_3_10/lib.sh. The protein count of
#                                       each file must equal n_proteins of inputs.tsv. NOT RE-CHECKED ON HPCC BY THE AUTHOR OF THIS DRAFT.
# The `#neg` copies exist only in the Kiss clustering FASTA. This job scores the negatives from the committed shared FASTA under their
# plain accession, so no `#neg` name is involved (the six accessions that are both a Kiss control and a shared negative are positives in
# one table and negatives in the other, as decided under D20).
# Pfam: PFAM_HMM (default the `current` link) must resolve to a directory named *Pfam$PFAM_RELEASE; sha256 as in the Task 5.4 job (a
# different sha256 is recorded and reported, a release-name mismatch stops the job). Module: hmmer/3.4.
#
# Submit (after the go-ahead; CODE_COMMIT is the merge commit of the PR that holds this script):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/to/in,OUT=/path/to/new_outdir,PFAM_RELEASE=38.2 \
#     -o /path/to/logs/%x.%j.log -e /path/to/logs/%x.%j.log analysis/wall_fusion_truth/pr6_family_specificity.sbatch
#   (the directory of the log must exist before sbatch; $OUT itself must NOT exist; optional: PROTEOME_DIR=/path.)
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J pr6_family_specificity
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the directory with kiss_sequences.faa) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm}"
PROTEOME_DIR="${PROTEOME_DIR:-/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat}"
EARLIER_PFAM_SHA256=4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e  # the earlier provenance (Pfam 38.2)
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/pr6_family_specificity.$$"
trap 'rm -rf "$TMP"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-4}"
WF=analysis/wall_fusion_truth
SCAN=docs/reports/data/pr6_pfam_scan_2026-10-22
CAND=docs/reports/data/pr6_family_candidates_2026-10-22
BASE=docs/reports/data/kiss_baseline_2026-10-22
REC=docs/reports/data/kiss_cluster_split_2026-10-22
NEG=data/controls/hard-negatives-shared
INPUTS=scripts/sorting_hat/refresh_3_10/inputs.tsv
FAMILY=data/sorting_hat/family_table.tsv
KEPT="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"  # one proteome per species (decision 9b)
N_POS_EXP=132   # the Kiss positives of members_split.tsv (recomputed from the committed table)
N_NEG_EXP=161   # the shared hard negatives (137 N1 and 24 N2)
N_FAM_EXP=30    # the candidate families (recomputed from the committed scan table)

fail() { echo "FATAL: $*" >&2; exit 1; }

# The code, the record tables and the candidate tables must be the reviewed ones.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$WF" "$SCAN" "$CAND" "$BASE" "$REC" "$NEG" "$INPUTS" "$FAMILY" \
    || fail "$WF, a record table or a candidate table changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$WF" "$SCAN" "$CAND" "$BASE" "$REC" "$NEG" "$INPUTS" "$FAMILY")" ] \
    || fail "uncommitted changes in $WF, a record table or a candidate table"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

# sidecar_ok FILE: the .sha256 next to FILE holds the sha256 of FILE
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
for f in "$PROJ_ROOT/$SCAN/scan_per_cluster.tsv" "$PROJ_ROOT/$BASE/clusters_missed.tsv" "$PROJ_ROOT/$REC/members_split.tsv" \
         "$PROJ_ROOT/$NEG/sequences.faa" "$PROJ_ROOT/$NEG/controls.tsv" "$PROJ_ROOT/$NEG/clusters.tsv" "$PROJ_ROOT/$INPUTS" \
         "$PROJ_ROOT/$CAND/candidate_families.tsv"; do
  [ -f "$f" ] || fail "$f is missing (the clone does not hold the merged record PRs)"
done
[ -f "$IN/kiss_sequences.faa" ] && [ -f "$IN/kiss_sequences.faa.sha256" ] || fail "$IN/kiss_sequences.faa or its .sha256 is missing"
sidecar_ok "$IN/kiss_sequences.faa" || fail "$IN/kiss_sequences.faa differs from its .sha256"
mkdir -p "$TMP"
PY=/usr/bin/python3.12

# 1. The candidate tables, made again from the committed scan table; they must equal the committed (reviewed) ones.
"$PY" -I "$PROJ_ROOT/$WF/pr6_family_specificity.py" candidates --scan "$PROJ_ROOT/$SCAN/scan_per_cluster.tsv" \
  --missed "$PROJ_ROOT/$BASE/clusters_missed.tsv" --family-table "$PROJ_ROOT/$FAMILY" --out "$TMP/cand"
for f in candidate_clusters.tsv candidate_families.tsv family_members.tsv families.txt; do
  cmp -s "$TMP/cand/$f" "$PROJ_ROOT/$CAND/$f" || { diff "$TMP/cand/$f" "$PROJ_ROOT/$CAND/$f" >&2 || true; fail "$f differs from the committed candidate table"; }
done
N_FAM=$(wc -l < "$TMP/cand/families.txt")
[ "$N_FAM" -eq "$N_FAM_EXP" ] || fail "$N_FAM candidate families, not $N_FAM_EXP"

# 2. The Kiss positives: exact match on the first token of the header (a `#neg` copy is a different name and is never taken).
tail -n +2 "$PROJ_ROOT/$REC/members_split.tsv" | awk -F'\t' '$3 == "positive" { print $1 }' | sort > "$TMP/pos_ids.txt"
N_IDS=$(wc -l < "$TMP/pos_ids.txt")
[ "$N_IDS" -eq "$N_POS_EXP" ] || fail "members_split.tsv lists $N_IDS positives, not $N_POS_EXP"
awk 'NR==FNR { want[$1] = 1; next }
     /^>/ { name = substr($1, 2); keep = (name in want); if (keep) seen[name] = 1 }
     keep { print }
     END { for (n in want) if (!(n in seen)) { print "MISSING " n > "/dev/stderr"; bad = 1 } exit bad }' \
  "$TMP/pos_ids.txt" "$IN/kiss_sequences.faa" > "$TMP/positives.faa" || fail "a positive of members_split.tsv is not in kiss_sequences.faa"
[ "$(grep -c '^>' "$TMP/positives.faa")" -eq "$N_POS_EXP" ] || fail "extracted a wrong number of positives"

# 3. The shared negatives (the committed FASTA): 161 sequences, named as the accessions of controls.tsv.
cp "$PROJ_ROOT/$NEG/sequences.faa" "$TMP/negatives.faa"
[ "$(grep -c '^>' "$TMP/negatives.faa")" -eq "$N_NEG_EXP" ] || fail "$NEG/sequences.faa does not hold $N_NEG_EXP sequences"
grep '^>' "$TMP/negatives.faa" | sed 's/^>//; s/[[:space:]].*//' | sort > "$TMP/neg_fa_ids.txt"
tail -n +2 "$PROJ_ROOT/$NEG/controls.tsv" | cut -f1 | sort > "$TMP/neg_ctl_ids.txt"
cmp -s "$TMP/neg_fa_ids.txt" "$TMP/neg_ctl_ids.txt" || fail "the names of $NEG/sequences.faa differ from the accessions of controls.tsv"

# 4. The kept proteomes: the files exist and the protein counts equal inputs.tsv (read only).
for p in $KEPT; do
  f="$PROTEOME_DIR/$p.faa"
  [ -f "$f" ] || fail "$f is missing"
  n_exp=$(awk -F'\t' -v p="$p" '$1 == p { print $3 }' "$PROJ_ROOT/$INPUTS")
  [ -n "$n_exp" ] || fail "$p is not in $INPUTS"
  n=$(grep -c '^>' "$f" || true)
  [ "$n" -eq "$n_exp" ] || fail "$f holds $n proteins, inputs.tsv says $n_exp"
done
echo "inputs ok: $N_FAM families, $N_POS_EXP positives, $N_NEG_EXP negatives, proteomes: $KEPT"

module load hmmer/3.4

# Pfam: the resolved path must be a directory named *Pfam$PFAM_RELEASE (only this stops the job).
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) fail "$RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" ;; esac
HMM_SHA="$(sha256sum "$RESOLVED" | cut -d' ' -f1)"
if [ "$HMM_SHA" = "$EARLIER_PFAM_SHA256" ]; then PFAM_SAME=yes; else
  PFAM_SAME=no
  echo "REPORT: sha256 of $RESOLVED is $HMM_SHA, not the $EARLIER_PFAM_SHA256 of the earlier jobs; recorded, the run goes on (release name $PFAM_RELEASE matches)."
fi

# 5. The candidate models, fetched by NAME (hmmfetch -f looks up by the model name, not by the accession) and checked.
hmmfetch -f "$PFAM_HMM" "$TMP/cand/families.txt" > "$TMP/models.hmm"
N=$(grep -c '^NAME' "$TMP/models.hmm" || true)
[ "$N" -eq "$N_FAM" ] || fail "hmmfetch found $N of $N_FAM models"
awk '/^NAME/{n=$2} /^ACC/{split($2,a,"."); print a[1]"\t"n}' "$TMP/models.hmm" | sort > "$TMP/fetched.tsv"
tail -n +2 "$TMP/cand/candidate_families.tsv" | cut -f1,2 | sort > "$TMP/expected.tsv"
diff -q "$TMP/fetched.tsv" "$TMP/expected.tsv" > /dev/null || { diff "$TMP/fetched.tsv" "$TMP/expected.tsv" >&2; fail "model names in the database differ from the candidate table"; }

# A model without a GA line is REPORTED, not skipped: hmmsearch --cut_ga stops on such a model, so the models are split.
# The split models are listed in no_ga_models.tsv and enter no rate (every Pfam-A family has a GA line; this is a guard).
awk -v ga="$TMP/with_ga.hmm" -v noga="$TMP/no_ga.hmm" -v list="$TMP/no_ga.txt" '
  { blk = blk $0 "\n" }
  /^ACC/ { split($2, a, "."); acc = a[1] }
  /^GA / { has = 1 }
  /^\/\// { if (has) printf "%s", blk > ga; else { printf "%s", blk > noga; print acc > list } blk = ""; has = 0; acc = "" }
' "$TMP/models.hmm"
touch "$TMP/with_ga.hmm" "$TMP/no_ga.hmm" "$TMP/no_ga.txt"
N_GA=$(grep -c '^NAME' "$TMP/with_ga.hmm" || true)
N_NOGA=$(grep -c '^NAME' "$TMP/no_ga.hmm" || true)
echo "models: $N_FAM fetched, $N_GA with a GA line, $N_NOGA without"
[ "$N_NOGA" -eq 0 ] || fail "$N_NOGA candidate model(s) have no GA line ($(tr '\n' ' ' < "$TMP/no_ga.txt")); a family without GA cannot be scored by --cut_ga: report to the owner"

# 6. hmmsearch --cut_ga of the candidate models over every FASTA.
search_ga() { # search_ga FASTA DOMTBL
  hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --domtblout "$2" "$TMP/with_ga.hmm" "$1"
  tail -n 3 "$2" | grep -qx '# \[ok\]' || fail "$2 has no '# [ok]' trailer"
}
search_ga "$TMP/positives.faa" "$TMP/domtbl_positives.txt"
search_ga "$TMP/negatives.faa" "$TMP/domtbl_negatives.txt"
PROT_ARGS=()
for p in $KEPT; do
  search_ga "$PROTEOME_DIR/$p.faa" "$TMP/domtbl_$p.txt"
  PROT_ARGS+=(--proteome "$p:$TMP/domtbl_$p.txt:$PROTEOME_DIR/$p.faa")
done

# 7. The tables (standard library only). The script refuses a non-empty output directory.
"$PY" -I "$PROJ_ROOT/$WF/pr6_family_specificity.py" summarize --candidates "$TMP/cand" \
  --members-split "$PROJ_ROOT/$REC/members_split.tsv" --neg-dir "$PROJ_ROOT/$NEG" \
  --domtbl-pos "$TMP/domtbl_positives.txt" --domtbl-neg "$TMP/domtbl_negatives.txt" \
  "${PROT_ARGS[@]}" --inputs "$PROJ_ROOT/$INPUTS" --out "$TMP/tables"

{
  echo "job: pr6_family_specificity (PR 6, owner decisions 1 and 2 of 2026-10-22)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "hmmer: $(hmmsearch -h | sed -n 2p)"
  echo "python: $($PY --version)"
  echo "pfam_hmm: $PFAM_HMM"
  echo "pfam_resolved: $RESOLVED"
  echo "pfam_release: $PFAM_RELEASE"
  echo "pfam_sha256: $HMM_SHA"
  echo "pfam_sha256_same_as_earlier_jobs: $PFAM_SAME (earlier: $EARLIER_PFAM_SHA256)"
  echo "input kiss_sequences.faa sha256 $(sha256sum < "$IN/kiss_sequences.faa" | cut -d' ' -f1)"
  echo "input negatives.faa sha256 $(sha256sum < "$TMP/negatives.faa" | cut -d' ' -f1)"
  for p in $KEPT; do echo "input $p.faa sha256 $(sha256sum < "$PROTEOME_DIR/$p.faa" | cut -d' ' -f1)"; done
  echo "families: $N_FAM (candidate rule: coverage >= 0.5 or E-value <= 1e-30, both inclusive)"
  echo "settings: hmmsearch --cut_ga with the candidate models over the Kiss positives, the shared negatives and the kept proteomes"
  echo "tier_changes: none; panel_changes: none; family_measure.tsv and family_table.tsv are not written"
} > "$TMP/provenance.txt"

# Write the output last, in a directory that did not exist. The FASTA files stay in $TMP (they are not outputs).
[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP"/tables/* "$OUT/"
cp "$TMP"/domtbl_*.txt "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
