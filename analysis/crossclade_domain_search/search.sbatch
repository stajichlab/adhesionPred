#!/bin/bash -l
# Cross-clade domain search, round 1 (strategy A of docs/reports/2026-11-25-plm-survey-rhodotorula-adhesion.md).
# Needs the owner's go-ahead for EACH submission. Heavy work is a SLURM job, never on the login node.
#
# Two methods on the same queries, decoys and targets:
#   1. Foldseek 3Di+AA search, with 3Di predicted from sequence by ProstT5 on the GPU (no structures needed);
#   2. MMseqs2 sequence search (-s 7.5), the sequence baseline the structure alphabet must beat.
# The reversed decoys (decoys.faa) are searched exactly like the queries; their hits give the empirical false-hit rate.
# Searches are exhaustive (small target set) and kept to E <= 10 so the summary can apply the declared cut-offs afterwards.
#
# Inputs:  $IN from build_inputs.py (queries.faa, decoys.faa, targets.faa, SHA256SUMS);
#          $WEIGHTS = ProstT5 weights from `foldseek databases ProstT5 <dir> <tmp>` (downloaded once, on the login node).
# Output:  $OUT/{foldseek,mmseqs}_{queries,decoys}.tsv, targets_3di.faa (reusable), provenance.txt. $OUT must NOT exist.
#
# Submit (after the go-ahead; CODE_COMMIT = the commit that holds this script; the log directory must exist):
#   sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=/path/in,WEIGHTS=/path/prostt5,OUT=/path/new_out \
#     -o /path/logs/%x.%j.log -e /path/logs/%x.%j.log analysis/crossclade_domain_search/search.sbatch
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH -J xclade_search
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone that holds this script) before sbatch}"
: "${CODE_COMMIT:?export CODE_COMMIT (the reviewed commit) before sbatch}"
: "${IN:?export IN (the build_inputs.py output directory) before sbatch}"
: "${WEIGHTS:?export WEIGHTS (the ProstT5 weights directory) before sbatch}"
: "${OUT:?export OUT (a new output directory) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/xclade_search.$$"
trap 'rm -rf "$TMP"' EXIT
PINNED="analysis/crossclade_domain_search"
FOLDSEEK_MODULE="foldseek/10-941cd33-gpu"
MMSEQS_MODULE="mmseqs2/17-b804f"
T="${SLURM_CPUS_PER_TASK:-16}"
FS_FMT="query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"
MM_FMT="query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"

fail() { echo "FATAL: $*" >&2; exit 1; }

check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- "$PINNED" \
    || fail "$PINNED changed since $CODE_COMMIT; the run would not be the reviewed one"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- "$PINNED")" ] || fail "uncommitted changes in $PINNED"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT"
}

check_pin
[ ! -e "$OUT" ] || fail "$OUT exists; nothing is overwritten"
[ -d "$WEIGHTS" ] || fail "ProstT5 weights directory $WEIGHTS is missing"
(cd "$IN" && sha256sum -c --quiet SHA256SUMS) || fail "$IN does not match its SHA256SUMS"
mkdir -p "$TMP"

module load "$FOLDSEEK_MODULE" "$MMSEQS_MODULE"
command -v foldseek >/dev/null || fail "foldseek not on PATH after module load $FOLDSEEK_MODULE"
command -v mmseqs >/dev/null || fail "mmseqs not on PATH after module load $MMSEQS_MODULE"
nvidia-smi -L || fail "no GPU visible"

# 1. ProstT5 3Di prediction (GPU), then Foldseek 3Di+AA search (CPU).
for s in queries decoys targets; do
  foldseek createdb "$IN/$s.faa" "$TMP/$s" --prostt5-model "$WEIGHTS" --gpu 1 --threads "$T" -v 1
done
for s in queries decoys; do
  foldseek search "$TMP/$s" "$TMP/targets" "$TMP/aln_$s" "$TMP/fs_tmp_$s" -e 10 --exhaustive-search 1 --threads "$T" -v 1
  foldseek convertalis "$TMP/$s" "$TMP/targets" "$TMP/aln_$s" "$TMP/foldseek_$s.tsv" --format-output "$FS_FMT" -v 1
done
foldseek lndb "$TMP/targets_h" "$TMP/targets_ss_h"
foldseek convert2fasta "$TMP/targets_ss" "$TMP/targets_3di.faa" -v 1
[ "$(grep -c '^>' "$TMP/targets_3di.faa")" -eq "$(grep -c '^>' "$IN/targets.faa")" ] || fail "3Di FASTA does not cover every target"

# 2. MMseqs2 sequence baseline.
for s in queries decoys; do
  mmseqs easy-search "$IN/$s.faa" "$IN/targets.faa" "$TMP/mmseqs_$s.tsv" "$TMP/mm_tmp_$s" \
    -s 7.5 -e 10 --exhaustive-search 1 --threads "$T" --format-output "$MM_FMT" -v 1
done

{
  echo "job: xclade_search (cross-clade domain search, round 1)"
  echo "commit: $(git -C "$PROJ_ROOT" rev-parse HEAD) (pinned: $CODE_COMMIT)"
  echo "node: $(hostname); gpu: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
  echo "module_list: $(module -t list 2>&1 | tr '\n' ' ')"
  echo "foldseek_version: $(foldseek version)"
  echo "mmseqs_version: $(mmseqs version)"
  echo "prostt5_weights: $WEIGHTS ($(du -sh "$WEIGHTS" | cut -f1); files sha256 $(cd "$WEIGHTS" && find . -type f | sort | xargs sha256sum | sha256sum | cut -d' ' -f1))"
  echo "foldseek: createdb --prostt5-model --gpu 1; search -e 10 --exhaustive-search 1; format $FS_FMT"
  echo "mmseqs: easy-search -s 7.5 -e 10 --exhaustive-search 1; format $MM_FMT"
  echo "input SHA256SUMS sha256 $(sha256sum < "$IN/SHA256SUMS" | cut -d' ' -f1)"
  for f in foldseek_queries foldseek_decoys mmseqs_queries mmseqs_decoys; do echo "rows $f: $(wc -l < "$TMP/$f.tsv")"; done
} > "$TMP/provenance.txt"

[ ! -e "$OUT" ] || fail "$OUT appeared while the job ran; nothing is overwritten"
mkdir -p "$OUT"
cp "$TMP"/foldseek_{queries,decoys}.tsv "$TMP"/mmseqs_{queries,decoys}.tsv "$TMP/targets_3di.faa" "$TMP/provenance.txt" "$OUT/"
echo "wrote $OUT"
