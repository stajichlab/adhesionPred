# Shared settings and checks for the Task 3.10 job scripts. It is sourced by them, never run alone.
# DRAFT (2026-10-17): nothing here has been run on HPCC. Read README.md first.
# Needs PROJ_ROOT: a clone of branch sorting-hat-engine-changes on HPCC. Paths are absolute.
: "${PROJ_ROOT:?export PROJ_ROOT (an HPCC clone of branch sorting-hat-engine-changes)}"
DATA_ROOT="${DATA_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
OLD="$DATA_ROOT/_workdir/sorting_hat"                     # old work directories: read only
NEW="${NEW:-$DATA_ROOT/_workdir/sorting_hat_refresh_3_10}" # new work directories: every output goes here
TAXDUMP="$OLD/taxdump"
PHASEC_DIR="$DATA_ROOT/_workdir/step1_compare/phasec"
R="$PROJ_ROOT/scripts/sorting_hat/refresh_3_10"
INPUTS="$R/inputs.tsv"
CONFIG="$PROJ_ROOT/src/cellsurface_sorting_hat/categories.yaml"
FAMILY_TABLE="$PROJ_ROOT/data/sorting_hat/family_table.tsv"
TRUTH_DIR="$PROJ_ROOT/analysis/calibration_truth/repeat_call_truth"
CALL_STATUS_DIR="$PROJ_ROOT/docs/reports/data/sorting_hat/call_status"
# The runs are tied to these values. Override only on purpose (and say so in the report).
CODE_COMMIT="${CODE_COMMIT:-51252bbdd11ec2fe5d81e2e6ec205b7b034f0752}"
FAMILY_SHA256="${FAMILY_SHA256:-f15c6f72dc79a79a2f4a128db311551353260efb80f8d63bd843e56cb75031af}"
CATEGORIES_SHA256="${CATEGORIES_SHA256:-68a35bad9dedcb8c34e0368153817edeb1ef489da0e12ed9d9243d4e52cc3721}"
PFAM_RELEASE="38.2"
PFAM_HMM_SHA256="4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e"  # of Pfam-A.hmm in the old provenance.json files

fail() { echo "FATAL: $*" >&2; exit 1; }

# prots COLUMN: the proteomes with "yes" in that column of inputs.tsv
prots() { awk -F'\t' -v col="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==col) c=i; next} $c=="yes"{print $1}' "$INPUTS"; }
# field PROT COLUMN: one value of inputs.tsv
field() { awk -F'\t' -v p="$1" -v col="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==col) c=i; next} $1==p{print $c}' "$INPUTS"; }

# The code and the data must be the ones the owner reviewed.
check_pin() {
  git -C "$PROJ_ROOT" cat-file -e "${CODE_COMMIT}^{commit}" || fail "commit $CODE_COMMIT is not in $PROJ_ROOT"
  git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD -- src data/sorting_hat analysis/cocci_repeats analysis/calibration_truth/repeat_call_truth \
    || fail "src, data/sorting_hat or a detector or truth file changed since $CODE_COMMIT; the runs would not be the reviewed ones"
  [ -z "$(git -C "$PROJ_ROOT" status --porcelain -- src data analysis scripts)" ] || fail "uncommitted changes in $PROJ_ROOT"
  [ "$(sha256sum "$FAMILY_TABLE" | cut -d' ' -f1)" = "$FAMILY_SHA256" ] || fail "family_table.tsv sha256 differs from $FAMILY_SHA256"
  [ "$(sha256sum "$CONFIG" | cut -d' ' -f1)" = "$CATEGORIES_SHA256" ] || fail "categories.yaml sha256 differs from $CATEGORIES_SHA256"
  echo "pin ok: HEAD $(git -C "$PROJ_ROOT" rev-parse HEAD); code pinned to $CODE_COMMIT; family_table $FAMILY_SHA256; categories $CATEGORIES_SHA256"
}

# make_venv TMP: a Python 3.12 virtual environment on node-local scratch (the package needs pyyaml and numpy;
# /usr/bin/python3.12 has neither). Sets PY and PYTHONPATH. The package itself is not installed: it runs from src.
make_venv() {
  /usr/bin/python3.12 -m venv "$1/venv"
  if [ -n "${WHEEL_DIR:-}" ]; then
    "$1/venv/bin/pip" install --quiet --no-index --find-links "$WHEEL_DIR" pyyaml numpy
  else
    "$1/venv/bin/pip" install --quiet pyyaml numpy
  fi
  PY="$1/venv/bin/python"
  export PYTHONPATH="$PROJ_ROOT/src"
  "$PY" -c 'import yaml, numpy; print("pyyaml", yaml.__version__, "numpy", numpy.__version__)'
}

# sidecar_ok FILE: the .sha256 next to FILE holds the sha256 of FILE
sidecar_ok() { [ "$(sha256sum < "$1" | cut -d' ' -f1)" = "$(tr -d '[:space:]' < "$1.sha256")" ]; }

# put SRC DST: copy under a temporary name, then rename; never replaces an existing file
put() {
  [ ! -e "$2" ] || fail "$2 exists; nothing is overwritten"
  cp "$1" "$2.tmp.$$" && mv "$2.tmp.$$" "$2"
}

# the number of proteins in the FASTA of a proteome must equal n_rows of its tm table (the same FASTA made both)
check_fasta() {
  local n rows
  n=$(grep -c '^>' "$OLD/$1.faa")
  rows=$(grep -o '"n_rows": [0-9]*' "$OLD/$1/modules/tm.json" | grep -o '[0-9]*$')
  [ "$n" = "$rows" ] || fail "$1: $OLD/$1.faa has $n proteins but modules/tm.json has n_rows $rows"
  local want; want=$(field "$1" n_proteins)
  if [ "$want" != "-" ] && [ "$n" != "$want" ]; then fail "$1: $n proteins; the committed provenance says $want"; fi
}
