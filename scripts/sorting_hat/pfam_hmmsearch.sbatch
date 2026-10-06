#!/bin/bash -l
# Search the Pfam models of the family table with hmmsearch --cut_ga (as analysis/cys_candidates/01_known_family_hmm.sh).
# Output: $WORKDIR/raw/pfam/domtbl.txt and provenance.json (path, resolved path, sha256 of Pfam-A.hmm).
# The domain table ends with "# [ok]" and its option line shows --cut_ga; the module command refuses a table without both.
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=...,FAMILY_TABLE=...,PFAM_RELEASE=38.2 pfam_hmmsearch.sbatch
# PFAM_RELEASE is stated by the person who submits; it is checked against the file name of the resolved path.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J csh_pfam
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
: "${FAMILY_TABLE:?export FAMILY_TABLE (data/sorting_hat/family_table.tsv) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_pfam.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/pfam"
mkdir -p "$TMP" "$OUT"
module load hmmer/3.4
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) echo "FATAL: $RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" >&2; exit 1 ;; esac
# hmmfetch -f looks up by the model NAME (column 2) in the pressed database, not by the accession
tail -n +2 "$FAMILY_TABLE" | cut -f2 | sort -u > "$TMP/keys.txt"
hmmfetch -f "$PFAM_HMM" "$TMP/keys.txt" > "$TMP/families.hmm"
N=$(grep -c '^NAME' "$TMP/families.hmm" || true)
[ "$N" -eq "$(wc -l < "$TMP/keys.txt")" ] || { echo "FATAL: hmmfetch found $N of $(wc -l < "$TMP/keys.txt") models" >&2; exit 1; }
awk '/^NAME/{n=$2} /^ACC/{split($2,a,"."); print a[1]"\t"n}' "$TMP/families.hmm" | sort > "$TMP/fetched.txt"
tail -n +2 "$FAMILY_TABLE" | cut -f1,2 | sort > "$TMP/table.txt"
diff -q "$TMP/fetched.txt" "$TMP/table.txt" > /dev/null || { echo "FATAL: model names in the database differ from the family table" >&2; diff "$TMP/fetched.txt" "$TMP/table.txt" >&2; exit 1; }
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
hmmsearch --cut_ga --cpu "${SLURM_CPUS_ON_NODE:-2}" --noali -o /dev/null --domtblout "$TMP/domtbl.txt" "$TMP/families.hmm" "$TMP/in.fasta"
tail -n 3 "$TMP/domtbl.txt" | grep -qx '# \[ok\]' || { echo "FATAL: domain table has no '# [ok]' trailer" >&2; exit 1; }
cp "$TMP/domtbl.txt" "$OUT/.tmp.domtbl.txt" && mv "$OUT/.tmp.domtbl.txt" "$OUT/domtbl.txt"
printf '{"pfam_hmm": "%s", "resolved": "%s", "release": "%s", "sha256": "%s", "hmmer": "%s"}\n' \
  "$PFAM_HMM" "$RESOLVED" "$PFAM_RELEASE" "$(sha256sum "$RESOLVED" | cut -d' ' -f1)" "$(hmmsearch -h | sed -n 2p)" > "$OUT/provenance.json"
