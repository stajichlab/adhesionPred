#!/bin/bash -l
# Task 3.10, job 0 (DRAFT, not submitted): make the new work directories, one per proteome of inputs.tsv.
# Copies (does not move) from the old work directories under _workdir/sorting_hat:
#   modules/ except the four Pfam modules (they are rebuilt by job 2), raw/pfam/{domtbl.txt,provenance.json},
#   and status/pfam_hydrophobin.json for the 8 proteomes that have a committed status file.
# Every copied module file is checked against its .sha256 file. Nothing in the old directories is changed.
# A new directory that already exists stops the job. Output: $NEW/<proteome>/{modules,raw/pfam,status}.
# Submit: sbatch --export=ALL,PROJ_ROOT=/path/to/clone -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log 00_stage_workdirs.sbatch
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=0:30:00
#SBATCH -J csh_r310_stage
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
. "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/lib.sh"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_r310_stage.$$"
trap 'rm -rf "$TMP" "$NEW"/.*.tmp.$$' EXIT
mkdir -p "$TMP" "$NEW"
check_pin
make_venv "$TMP"
for p in $(tail -n +2 "$INPUTS" | cut -f1); do
  echo "== $p"
  [ ! -e "$NEW/$p" ] || fail "$NEW/$p exists; nothing is overwritten"
  [ -d "$OLD/$p/modules" ] || fail "$OLD/$p/modules not found"
  check_fasta "$p"
  W="$TMP/$p"
  mkdir -p "$W/modules" "$W/raw/pfam" "$W/status"
  for f in "$OLD/$p"/modules/*; do
    b=$(basename "$f")
    case "$b" in *.sha256) continue ;; pfam_adhesion.*|pfam_allergen.*|pfam_hydrophobin.*|pfam_hsba.*) continue ;; esac
    sidecar_ok "$f" || fail "$f does not match its .sha256 file"
    cp "$f" "$f.sha256" "$W/modules/"
  done
  for b in domtbl.txt provenance.json; do cp "$OLD/$p/raw/pfam/$b" "$W/raw/pfam/$b"; done
  grep -q "\"sha256\": \"$PFAM_HMM_SHA256\"" "$W/raw/pfam/provenance.json" || fail "$p: Pfam-A.hmm sha256 in provenance.json is not $PFAM_HMM_SHA256"
  grep -q "\"release\": \"$PFAM_RELEASE\"" "$W/raw/pfam/provenance.json" || fail "$p: Pfam release in provenance.json is not $PFAM_RELEASE"
  if [ "$(field "$p" hyd)" = yes ]; then
    # the committed status file and the one on HPCC must hold the same data, or the committed file is not the one the engine reads
    s="$OLD/$p/status/pfam_hydrophobin.json"
    sidecar_ok "$s" || fail "$s does not match its .sha256 file"
    "$PY" "$R/refresh_checks.py" json-equal "$s" "$CALL_STATUS_DIR/hydrophobin/pfam_hydrophobin.$p.json"
    cp "$s" "$s.sha256" "$W/status/"
  fi
  ls -l "$W/modules" | tail -n +2
  cp -a "$W" "$NEW/.$p.tmp.$$"
  mv "$NEW/.$p.tmp.$$" "$NEW/$p"
done
echo "staged $(tail -n +2 "$INPUTS" | wc -l) work directories under $NEW"
