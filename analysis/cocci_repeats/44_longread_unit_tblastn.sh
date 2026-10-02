#!/usr/bin/bash
# Count SOWgp repeat units directly in the 5 UArizona long-read assemblies, without gene calls.
# Query: one internal RS unit (CIMG_04613 aa 130-176, 47 aa). tblastn, SEG off, E <= 1e-5.
# Output: longread_unit_tblastn/<strain>.tsv (outfmt 6) and <strain>.scaffolds.fa.fai is not written.
set -euo pipefail
module load ncbi-blast/2.13.0+
LONGREAD=/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc
OUT=longread_unit_tblastn
mkdir -p $OUT
: "${SCRATCH:?SCRATCH not set}"
cat > $SCRATCH/unit.faa <<'Q'
>RS_unit2_aa130-176
PTDCYGDCEDGYDYSPPPPPKKYGDCDYDDGYCDGPSKTSMKPEPPK
Q
for strain in CiB10637 CiB10992 VFC140 Cpos1038 Cpos3700; do
    fa=$(ls $LONGREAD/$strain/*.scaffolds.fa | head -1)
    makeblastdb -in $fa -dbtype nucl -out $SCRATCH/db_$strain -logfile /dev/null
    tblastn -query $SCRATCH/unit.faa -db $SCRATCH/db_$strain -seg no -evalue 1e-5 \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
        > $OUT/$strain.tsv
    echo "$strain $(wc -l < $OUT/$strain.tsv) hits"
done
