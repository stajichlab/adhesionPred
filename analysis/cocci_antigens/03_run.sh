#!/bin/bash
set -uo pipefail
PROJ_ROOT="${PROJ_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
source "$PROJ_ROOT/analysis/_common/paths.sh"

cd "$WORKDIR"
W=cocci_antigens
PAN=$COCCI_PANGENOME
module load MMseqs2 2>/dev/null
[ -s iedb_antigens.fa ] || python3 "$PROJ_ROOT/analysis/cocci_antigens/00_download_iedb_antigens.py" \
  --out iedb_antigens.fa --meta iedb_antigens_meta.tsv
[ -s $W/iedb.m8 ] || mmseqs easy-search $PAN/input_run2/CimmitisRS_FungiDB.fasta iedb_antigens.fa $W/iedb.m8 $W/tmp_iedb \
  --min-seq-id 0.3 -c 0.5 --cov-mode 0 -v 1 --max-seqs 3 --format-output query,target,pident,qcov >/dev/null 2>&1
echo "IEDB homology hits: $(wc -l < $W/iedb.m8)"
python3 "$PROJ_ROOT/analysis/cocci_antigens/02_score_antigens.py" --work $W --antigens /dev/null --iedb-hits $W/iedb.m8 --out $W/cocci_antigen_ranking.tsv
