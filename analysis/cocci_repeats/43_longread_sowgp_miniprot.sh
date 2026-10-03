#!/usr/bin/bash
# Splice-aware protein-to-genome alignment of the SOWgp seed proteins to the 5 UArizona long-read
# assemblies (miniprot). Independent of the funannotate gene calls.
# Output: longread_miniprot/<strain>.sowgp.gff  (GFF3 with ##STA translated protein lines)
set -euo pipefail
module load miniprot
LONGREAD=/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc
OUT=longread_miniprot
mkdir -p $OUT
for strain in CiB10637 CiB10992 VFC140 Cpos1038 Cpos3700; do
    fa=$(ls $LONGREAD/$strain/*.scaffolds.fa | head -1)
    miniprot -t 8 --gff -G 1000 $fa sowgp_seed.fa > $OUT/$strain.sowgp.gff 2> $OUT/$strain.sowgp.log
    echo "$strain $(grep -c -P '\tmRNA\t' $OUT/$strain.sowgp.gff) mRNA records"
done
