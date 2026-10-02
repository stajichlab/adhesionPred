#!/bin/bash -l
# Cut six Pfam models and search them against the Coccidioides proteomes (--cut_ga).
# Models: PF05730 (CFEM), PF04681 (Bys1), PF01185, PF06766 and PF28987 (hydrophobins, PF28987 =
# DewD), PF28404 (PRA3 family, ARB_05178).
# Small job (about one minute for 7 proteomes). Not submitted by this script.
# Environment (no path is derived from the script location):
#   OUT_DIR    output directory (required)
#   MANIFEST   TSV with columns name, fasta (extra columns are ignored); default
#              $OUT_DIR/proteomes.tsv
#   PFAM_HMM   default: the central database link
#              /bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm
#              The link is resolved and the release directory name is recorded.
#              The shared funannotate_db copy was Pfam 38.0 when this was first run and has
#              no PF28404. Do not use it unless you checked its release.
#   CPU        default 4
# Output: <name>.domtbl.gz per proteome, known_families.hmm, pfam_provenance.json.
# All outputs are written to a temporary directory and moved only when every step worked.
# Optional SLURM header, only if you run this as a job:
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --time=1:00:00
set -euo pipefail

stop() { echo "STOP: $*" >&2; exit 2; }

[ -n "${OUT_DIR:-}" ] || stop "OUT_DIR is not set"
PFAM_HMM="${PFAM_HMM:-/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm}"
MANIFEST="${MANIFEST:-$OUT_DIR/proteomes.tsv}"
CPU="${CPU:-4}"
ACCS="PF05730 PF04681 PF01185 PF06766 PF28987 PF28404"

[ -s "$PFAM_HMM" ] || stop "Pfam-A.hmm not found: $PFAM_HMM"
[ -s "$MANIFEST" ] || stop "manifest not found: $MANIFEST"
module load hmmer/3.4 || stop "cannot load hmmer/3.4"
mkdir -p "$OUT_DIR" || stop "cannot create OUT_DIR: $OUT_DIR"
TMP=$(mktemp -d "$OUT_DIR/.tmp_01.XXXXXX") || stop "OUT_DIR is not writable: $OUT_DIR"
trap 'rm -rf "${TMP:?}"' EXIT

REAL=$(readlink -f "$PFAM_HMM")
RELEASE=$(basename "$(dirname "$REAL")")

# One pass over the file to read the exact accession (with version) of each model.
KEYS="$TMP/keys.txt"
: > "$KEYS"
grep -E "^ACC +(PF05730|PF04681|PF01185|PF06766|PF28987|PF28404)\." "$PFAM_HMM" \
  | awk '{print $2}' > "$KEYS" || true
for acc in $ACCS; do
  grep -q "^${acc}\." "$KEYS" || stop "$acc is not in $PFAM_HMM (release $RELEASE)"
done
hmmfetch -f "$PFAM_HMM" "$KEYS" > "$TMP/known_families.hmm" 2> "$TMP/hmmfetch.err" \
  || stop "hmmfetch failed on $PFAM_HMM (the database may be changing): $(tail -2 "$TMP/hmmfetch.err" | tr '\n' ' ')"
N=$(grep -c '^NAME' "$TMP/known_families.hmm")
[ "$N" -eq 6 ] || stop "expected 6 models, fetched $N"
echo "release: $RELEASE; models: $N"

while IFS=$'\t' read -r name fasta _rest; do
  case "$name" in ''|'#'*) continue ;; esac
  [ -s "$fasta" ] || stop "FASTA not found for $name: $fasta"
  hmmsearch --cut_ga --cpu "$CPU" --noali -o /dev/null --domtblout "$TMP/$name.domtbl" \
      "$TMP/known_families.hmm" "$fasta" || stop "hmmsearch failed for $name"
  gzip -n -c "$TMP/$name.domtbl" > "$TMP/$name.domtbl.gz"
  rm -f "$TMP/$name.domtbl"
  echo "$name: $(zcat "$TMP/$name.domtbl.gz" | grep -vc '^#') domain hits"
done < "$MANIFEST"

SHA=$(sha256sum "$TMP/known_families.hmm" | awk '{print $1}')
{
  printf '{\n  "pfam_hmm_path": "%s",\n  "pfam_hmm_resolved": "%s",\n' "$PFAM_HMM" "$REAL"
  printf '  "pfam_release_dir": "%s",\n' "$RELEASE"
  printf '  "known_families_sha256": "%s",\n  "models": [\n' "$SHA"
  first=1
  while read -r acc name; do
    [ $first -eq 1 ] || printf ',\n'
    first=0
    printf '    {"acc": "%s", "name": "%s"}' "$acc" "$name"
  done < <(awk '/^NAME/{n=$2} /^ACC/{print $2, n}' "$TMP/known_families.hmm")
  printf '\n  ]\n}\n'
} > "$TMP/pfam_provenance.json"

rm -f "$TMP/keys.txt" "$TMP/hmmfetch.err"
for f in "$TMP"/*; do mv -f "$f" "$OUT_DIR/$(basename "$f")"; done
echo "done: $OUT_DIR"
