#!/usr/bin/bash -l
module load samtools

for file in $(ls *.sampled_ids.txt)
do
	in=$(basename $file .sampled_ids.txt)
	samtools faidx -r $file $in.fasta > $in.negative.pep
done
