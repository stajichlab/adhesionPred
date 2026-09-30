#!/usr/bin/bash -l
# Ask the SOWgp locus question of the three long-read assemblies where no SOWgp protein was
# called (report section 4.2: VFC140, Cpos1038, Cpos3700): is there sequence similarity to the
# SOWgp seed anywhere in the assembly?
# If the protein call is missing because:
#   a) the locus is absent        -> no tblastn hit
#   b) the gene model failed      -> tblastn hit, but no protein in the proteome at that position
#   c) the repeat array collapsed -> tblastn hit spanning the array with truncation
# then the tblastn hit is the discriminator between "real absence" and "annotation artifact".
# Uses the full-length seed proteins (>=300 aa: RS, CiB10637, CiB10992, Silveira + published).
# Output: sowgp_tblastn_{VFC140,Cpos1038,Cpos3700}.out
# Run: sbatch -N 1 -c 8 --mem 8G --time 1:0:0 -p short 07_sowgp_locus_tblastn.sh

module load ncbi-blast/2.14.1
set -euo pipefail
cd "$(dirname "$0")"
LONGREAD=/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc

# build a small db of the assembly contigs for each strain, then tblastn with the SOWgp seeds
for strain in VFC140 Cpos1038 Cpos3700; do
    fa=$(ls "$LONGREAD/$strain"/*.scaffolds.fa 2>/dev/null | head -1)
    [ -n "$fa" ] || fa=$(ls "$LONGREAD/$strain"/*.contigs.fsa 2>/dev/null | head -1)
    echo "== $strain ($fa)"
    makeblastdb -in "$fa" -dbtype nucl -out "/tmp/sowgp_${strain}" -title "$strain" -nowarn
    tblastn -query sowgp_seed.fa -db "/tmp/sowgp_${strain}" \
        -out sowgp_tblastn_${strain}.out -outfmt 6 \
        -evalue 1e-5 -seg yes -max_target_seqs 50 -num_threads 8
    echo "   hits: $(wc -l < sowgp_tblastn_${strain}.out)"
done
