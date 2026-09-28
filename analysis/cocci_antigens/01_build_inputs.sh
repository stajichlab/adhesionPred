#!/bin/bash
# Build the inputs the antigen ranking needs, on HPCC.
#   - ID map from the pangenome's FungiDB reference to the Fungi_5k proteome (for SignalP/Pfam)
#   - cross-reactivity searches against the fungi that actually confound cocci serology
#   - anchor -> reference gene mapping
# Run:  srun -p epyc -c 8 --mem 16G -t 60 ./01_build_inputs.sh
set -euo pipefail
PROJ_ROOT="${PROJ_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
source "$PROJ_ROOT/analysis/_common/paths.sh"

OUT=${OUT:-$WORKDIR/cocci_antigens}
PAN=$COCCI_PANGENOME
F5K=$FUNGI5K_INPUT
mkdir -p "$OUT" && cd "$OUT"
module load MMseqs2 2>/dev/null || true
T=${SLURM_CPUS_PER_TASK:-8}
REF=$PAN/input_run2/CimmitisRS_FungiDB.fasta

# 1. map pangenome reference ids -> Fungi_5k ids (same genome, different annotation ids)
F5K_RS=$(ls $F5K/Coccidioides_immitis_RS.proteins.fa 2>/dev/null || ls $F5K/*immitis_RS*.proteins.fa | head -1)
echo "Fungi_5k RS proteome: $F5K_RS"
[ -s idmap.m8 ] || mmseqs easy-search "$REF" "$F5K_RS" idmap.m8 tmp_id \
   --min-seq-id 0.95 -c 0.9 --cov-mode 0 -v 1 --max-seqs 5 \
   --format-output query,target,pident,qcov >/dev/null

# 2. cross-reactivity: the fungi that confound coccidioidomycosis serology
: > crossreact.m8
for PAT in Histoplasma Blastomyces Paracoccidioides Aspergillus_fumigatus; do
  FA=$(ls $F5K/${PAT}*.proteins.fa 2>/dev/null | head -1) || true
  [ -z "${FA:-}" ] && { echo "  no proteome for $PAT"; continue; }
  echo "  cross-reactivity vs $(basename "$FA")"
  mmseqs easy-search "$REF" "$FA" "cr_$PAT.m8" "tmp_cr_$PAT" \
     --min-seq-id 0.3 -c 0.5 --cov-mode 0 -v 1 --max-seqs 3 \
     --format-output query,target,pident,qcov,evalue >/dev/null
  awk -v g="$PAT" 'BEGIN{OFS="\t"}{print g,$0}' "cr_$PAT.m8" >> crossreact.m8
done

# 3. human cross-reactivity (kept from the original score)
[ -s human.faa ] || curl -sL "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome:UP000005640" -o human.faa
[ -s human.m8 ] || mmseqs easy-search "$REF" human.faa human.m8 tmp_hs \
   --min-seq-id 0.4 -c 0.5 --cov-mode 0 -v 1 --max-seqs 3 --format-output query,target,pident,qcov >/dev/null

# 4. anchors -> reference genes
[ -s anchors.m8 ] || mmseqs easy-search ../cocci_anchors.fa "$REF" anchors.m8 tmp_anc \
   --min-seq-id 0.4 -c 0.4 --cov-mode 0 -v 1 --format-output query,target,pident,qcov >/dev/null

wc -l idmap.m8 crossreact.m8 human.m8 anchors.m8
