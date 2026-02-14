#!/usr/bin/bash -l
module load fasta

QUERY=positive_to_skip.pep
for a in $(ls *AnnotatedProteins.fasta)
do
	fasta36 -m 8c -E 1e-5 $QUERY $a > $(basename $a .fasta).positivesearch.FASTA.tab
done
