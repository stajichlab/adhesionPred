#!/bin/bash
set -uo pipefail
cd /bigdata/stajichlab/jstajich/projects/adhesionPred_review
W=cocci_antigens
PAN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome
module load MMseqs2 2>/dev/null
[ -s $W/iedb.m8 ] || mmseqs easy-search $PAN/input_run2/CimmitisRS_FungiDB.fasta iedb_antigens.fa $W/iedb.m8 $W/tmp_iedb \
  --min-seq-id 0.3 -c 0.5 --cov-mode 0 -v 1 --max-seqs 3 --format-output query,target,pident,qcov >/dev/null 2>&1
echo "IEDB homology hits: $(wc -l < $W/iedb.m8)"
python3 02_score_antigens.py --work $W --antigens /dev/null --iedb-hits $W/iedb.m8 --out $W/cocci_antigen_ranking.tsv
