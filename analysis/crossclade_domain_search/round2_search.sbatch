#!/bin/bash -l
# Cross-clade domain search, round 2: structure search with real AlphaFold DB models (round2_inputs.py output).
# Needs the owner's go-ahead for EACH submission. CPU only.
#
#   1. Reseek v3 (-verysensitive, -stats superfamily), P-values up to 1 (cut-offs applied later). Reseek keeps a k-mer prefilter
#      even in -verysensitive mode (lowest threshold), so pairs that fail it have no row; the summary treats them as no hit;
#   2. Foldseek 3Di+AA on the AFDB models (CPU module), exhaustive, -e 10, --mask 0, same columns as round 1.
# The third method of the 2-of-3 rule is round 1's ProstT5 Foldseek (no re-run).
#
# Inputs:  $IN2 from round2_inputs.py (queries/, targets/ trimmed to pLDDT >= 50, models.tsv, SHA256SUMS); $RESEEK = reseek v3.01.
# Output:  $OUT/{reseek,foldseek_af}.tsv, provenance.txt. $OUT must NOT exist.
#
# Submit (after the go-ahead; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN2=/path/in2,RESEEK=/path/reseek,OUT=/path/new_out \
#     -o /path/logs/%x.%j.log -e /path/logs/%x.%j.log analysis/crossclade_domain_search/round2_search.sbatch
#SBATCH -p exfab
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH -J xclade_round2
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN2:?export IN2 (the round2_inputs.py output directory) before sbatch}"
: "${RESEEK:?export RESEEK (the reseek binary) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/xclade_round2.$$"
trap 'rm -rf "$TMP"' EXIT
PINNED="analysis/crossclade_domain_search"
FOLDSEEK_MODULE="foldseek/10-941cd33"
T="${SLURM_CPUS_PER_TASK:-16}"
FS_FMT="query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"
RS_COLS="query+target+qlo+qhi+ql+tlo+thi+tl+pctid+pvalue"

fail() { echo "FATAL: $*" >&2; exit 1; }

git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$PINNED" || fail "$PINNED changed since $CODE_COMMIT"
[ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$PINNED")" ] || fail "uncommitted changes in $PINNED"
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -x "$RESEEK" ] || fail "$RESEEK is not executable"
(cd "$IN2" && sha256sum -c --quiet SHA256SUMS) || fail "$IN2 does not match its SHA256SUMS"
mkdir -p "$TMP"
module load "$FOLDSEEK_MODULE"
command -v foldseek >/dev/null || fail "foldseek not on PATH after module load $FOLDSEEK_MODULE"

# 1. Reseek.
#
# Job 32502912 (2026-11-27) died here in 36 s: -convert took the targets DIRECTORY happily (1,138
# chains) but -search on the queries directory failed with
#   "ReadStdioFile failed, attempted 4 bytes, read 0 bytes, errno=21"   (errno 21 = EISDIR)
# The Reseek docs say a directory is a valid STRUCTS argument for both -search and -db, so this
# looks like a v3.01 bug rather than a usage error. The binary needs CPU features the login node
# lacks ("Illegal instruction" there), so the query form cannot be probed outside a job.
#
# Therefore: try the documented query forms in order and use the first that produces hits, and
# record which one worked in the provenance. The failing mode costs ~30 s, so the fallback is cheap.
#   .files  - documented list convention: "a text file containing filenames, one per line"
#   .bcb    - allowed by the STRUCTS grammar; the form -db already requires
#   dir     - what the docs say should work; kept last so a fixed Reseek is still used as intended
"$RESEEK" -convert "$IN2/targets" -bcb "$TMP/targets.bcb" -log "$TMP/reseek_convert.log"
find "$IN2/queries" -name '*.pdb' | sort > "$TMP/queries.files"
[ -s "$TMP/queries.files" ] || fail "no .pdb files in $IN2/queries"
"$RESEEK" -convert "$IN2/queries" -bcb "$TMP/queries.bcb" -log "$TMP/reseek_convert_queries.log" \
  || echo "note: -convert of the queries failed; the .bcb query form will be skipped" >&2

RESEEK_QUERY_FORM=""
for form in files bcb dir; do
  case "$form" in
    files) q="$TMP/queries.files" ;;
    bcb)   q="$TMP/queries.bcb" ;;
    dir)   q="$IN2/queries" ;;
  esac
  [ -e "$q" ] || continue
  rm -f "$TMP/reseek.tsv"
  if "$RESEEK" -search "$q" -db "$TMP/targets.bcb" -verysensitive -stats superfamily -pvalue 1 \
       -columns "$RS_COLS" -output "$TMP/reseek.tsv" -threads "$T" \
       -log "$TMP/reseek_search_$form.log" && [ -s "$TMP/reseek.tsv" ]; then
    RESEEK_QUERY_FORM="$form"
    cp "$TMP/reseek_search_$form.log" "$TMP/reseek_search.log"
    echo "reseek query form that worked: $form"
    break
  fi
  echo "reseek query form '$form' produced no hits; trying the next" >&2
done
[ -n "$RESEEK_QUERY_FORM" ] || fail "every reseek query form failed; see $TMP/reseek_search_*.log"

# 2. Foldseek on the AFDB models.
foldseek easy-search "$IN2/queries" "$IN2/targets" "$TMP/foldseek_af.tsv" "$TMP/fs_tmp" \
  -e 10 --exhaustive-search 1 --mask 0 --threads "$T" --format-output "$FS_FMT" -v 1

{
  echo "job: xclade_round2 (cross-clade domain search, round 2)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "node: $(hostname)"
  echo "reseek: $RESEEK sha256 $(sha256sum < "$RESEEK" | cut -d' ' -f1); $("$RESEEK" -version 2>&1 | head -1 || true)"
  echo "reseek: -verysensitive -stats superfamily -pvalue 1 -columns $RS_COLS"
  echo "reseek query form: $RESEEK_QUERY_FORM (queries given as a .files list, a .bcb, or a directory)"
  echo "foldseek_version: $(foldseek version) ($FOLDSEEK_MODULE)"
  echo "foldseek: easy-search -e 10 --exhaustive-search 1 --mask 0; format $FS_FMT"
  echo "input SHA256SUMS sha256 $(sha256sum < "$IN2/SHA256SUMS" | cut -d' ' -f1)"
  echo "queries: $(ls "$IN2/queries" | wc -l); targets: $(ls "$IN2/targets" | wc -l)"
  for f in reseek foldseek_af; do echo "rows $f: $(wc -l < "$TMP/$f.tsv"); queries with rows: $(cut -f1 "$TMP/$f.tsv" | sort -u | wc -l)"; done
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP/reseek.tsv" "$TMP/foldseek_af.tsv" "$TMP/provenance.txt" "$TMP"/reseek_*.log "$OUT/"
echo "wrote $OUT"
