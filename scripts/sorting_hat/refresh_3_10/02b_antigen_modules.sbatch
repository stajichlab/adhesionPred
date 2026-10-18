#!/bin/bash -l
# Task 3.10, job 2b (run on HPCC on 2026-10-17; results: docs/reports/2026-10-18-task-3-10-refresh-results.md): write the `antigen_lookup` module for the 8 run proteomes that have none
# (every run proteome except Cimm_RS). Outside Coccidioides the module writes state `not_applicable` for every
# protein (modules/lookups.py antigen_rows), so the engine no longer sees the call cocci_specificity_rank_top15 as
# not_assessable (module absent) and D5 stops turning the other_* calls and primary_label into not_assessable.
# Owner decision 2026-10-17, option (a): no code or config change. The structural fix is issue #93 (not done here).
# It runs exactly the command that made the Cimm_RS record (same ranking, same protein map, applicable taxon 246410),
# with the proteome FASTA and --taxon of each proteome. Nothing else differs. Checks:
#   - the sha256 of the two inputs are those that reproduce the artefact_hash of the Cimm_RS record;
#   - the rebuilt params_hash and artefact_hash of each new record equal those of the Cimm_RS record;
#   - every row is not_applicable (na_invalid rows, for a protein that is not a valid sequence, are counted and shown);
#     no protein is ok, not_in_reference or error, so no protein is called;
#   - the other module files of the work directory are unchanged (sha256 before and after).
# Nothing is written to $NEW unless every proteome passes. An existing antigen_lookup file stops the job.
# Output: $NEW/<proteome>/modules/antigen_lookup.{json,tsv.gz} (+ .sha256). Needs job 0 (and job 2 for the work directories to be final).
# Submit: sbatch --export=ALL,PROJ_ROOT=/path/to/clone,WHEEL_DIR=/path/to/wheels -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log 02b_antigen_modules.sbatch
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=0:20:00
#SBATCH -J csh_r310_antigen
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT before sbatch}"
. "$PROJ_ROOT/scripts/sorting_hat/refresh_3_10/lib.sh"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_r310_antigen.$$"
trap 'rm -rf "$TMP" "$NEW"/*/modules/.*.tmp.$$' EXIT
mkdir -p "$TMP"
# Inputs of the Cimm_RS record (found 2026-10-17; their two sha256 reproduce its artefact_hash 040ba813...)
RANKING="$PROJ_ROOT/analysis/cocci_antigens/cocci_antigen_ranking.tsv"
RANKING_SHA256="cf04b9d4d26f25b9c78cd073276aebf4bd19966693232a1e3c5903b5369318bc"
PROTEIN_MAP="$DATA_ROOT/_workdir/cocci_spherule/ref/protein_map.tsv"
PROTEIN_MAP_SHA256="ca0615211e865559207f98f00e3c676769552331871625ef3dbcea7483c83806"
APPLICABLE=246410
REF="$NEW/Cimm_RS/modules/antigen_lookup.json"
check_pin
[ "$(sha256sum "$RANKING" | cut -d' ' -f1)" = "$RANKING_SHA256" ] || fail "$RANKING sha256 differs from $RANKING_SHA256"
[ "$(sha256sum "$PROTEIN_MAP" | cut -d' ' -f1)" = "$PROTEIN_MAP_SHA256" ] || fail "$PROTEIN_MAP sha256 differs from $PROTEIN_MAP_SHA256"
[ -f "$REF" ] || fail "$REF not found (the Cimm_RS record is the reference)"
make_venv "$TMP"
PROTS="Afum_A1163 Afum_Af293_Fungi5k Afum_Af293_UniProt Afum_W72310 Bder_ER3 Calb_SC5314 Scer_S288C Scer_S288C_nodubious"
# the identity of the reference, from its own json
ref_id() { "$PY" -I -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["params_hash"], d["artefact_hash"], d["params"]["applicable_taxa"], d["params"]["genes_with_several_ranking_rows"])' "$1"; }
REF_ID="$(ref_id "$REF")"; echo "reference (Cimm_RS): $REF_ID"
# 1. build in scratch
for p in $PROTS; do
  echo "== $p"
  W="$NEW/$p"
  [ -d "$W/modules" ] || fail "$W/modules not found; run job 0"
  [ "$(field "$p" run)" = yes ] && [ "$p" != Cimm_RS ] || fail "$p is not a run proteome of inputs.tsv"
  for ext in json json.sha256 tsv.gz tsv.gz.sha256; do [ ! -e "$W/modules/antigen_lookup.$ext" ] || fail "$W/modules/antigen_lookup.$ext exists; nothing is overwritten"; done
  check_fasta "$p"
  ( cd "$W/modules" && sha256sum -- *.json *.tsv.gz | sort -k2 ) > "$TMP/$p.before"
  mkdir -p "$TMP/$p"
  "$PY" -m cellsurface_sorting_hat.modules.cli antigen --fasta "$OLD/$p.faa" --workdir "$TMP/$p" --taxon "$(field "$p" taxon)" \
    --applicable-taxa "$APPLICABLE" --ranking "$RANKING" --protein-map "$PROTEIN_MAP"
  M="$TMP/$p/modules"
  [ -f "$M/antigen_lookup.json" ] && [ -f "$M/antigen_lookup.tsv.gz" ] || fail "$p: the module was not written"
  [ "$(ls "$M" | wc -l)" = 4 ] || fail "$p: the command wrote files other than the four of antigen_lookup: $(ls "$M" | tr '\n' ' ')"
  [ "$(ref_id "$M/antigen_lookup.json")" = "$REF_ID" ] || fail "$p: identity differs from the Cimm_RS record: $(ref_id "$M/antigen_lookup.json") (reference $REF_ID)"
  sidecar_ok "$M/antigen_lookup.json" && sidecar_ok "$M/antigen_lookup.tsv.gz" || fail "$p: a .sha256 does not match its file"
  "$PY" -I - "$M/antigen_lookup.tsv.gz" "$M/antigen_lookup.json" "$(field "$p" n_proteins)" <<'PYEOF' || fail "$p: the antigen table is not all not_applicable"
import collections, gzip, json, sys
table, record, want = sys.argv[1], sys.argv[2], sys.argv[3]
rows = [line.rstrip("\n").split("\t") for line in gzip.open(table, "rt", encoding="utf-8")]
head, body = rows[0], rows[1:]
state = head.index("state")
counts = collections.Counter(r[state] for r in body)
meta = json.load(open(record))
print("states:", dict(counts), "n_rows:", len(body), "run_state:", meta["run_state"], "note:", meta["note"])
bad = {k: v for k, v in counts.items() if k not in ("not_applicable", "na_invalid")}
ok = (not bad and counts["not_applicable"] > 0 and len(body) == meta["n_rows"] and meta["run_state"] == "ok"
      and (want == "-" or len(body) == int(want)))
sys.exit(0 if ok else 1)
PYEOF
done
# 2. copy (temporary name, then rename; never replaces a file)
for p in $PROTS; do
  for ext in json.sha256 tsv.gz.sha256 json tsv.gz; do
    [ ! -e "$NEW/$p/modules/antigen_lookup.$ext" ] || fail "$NEW/$p/modules/antigen_lookup.$ext exists; nothing is overwritten"
    t="$NEW/$p/modules/.antigen_lookup.$ext.tmp.$$"
    cp "$TMP/$p/modules/antigen_lookup.$ext" "$t" && mv "$t" "$NEW/$p/modules/antigen_lookup.$ext"
  done
  ( cd "$NEW/$p/modules" && sha256sum -- $(ls *.json *.tsv.gz | grep -v '^antigen_lookup\.') | sort -k2 ) > "$TMP/$p.after"
  diff -q <(grep -v ' antigen_lookup\.' "$TMP/$p.before") "$TMP/$p.after" > /dev/null || fail "$p: another module file changed"
  echo "wrote antigen_lookup in $NEW/$p/modules (other modules unchanged)"
done
