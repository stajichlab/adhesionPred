module load samtools
module load fasta
module load phmmer
phmmer --domtblout FLO11_S288C.allScer.domtbl FLO11_S288C.pep scer_proteins.faa > FLO11_S288C.allScer.phmmer
fasta36 -m 8c -E 1e-10 FLO11_S288C.pep scer_proteins.faa > FLO11_S288C.FASTA.tab
python filter_blast_results.py FLO11_search.FASTA.tab | cut -f2 | sort | uniq | samtools faidx -r - scer_proteins.faa > FLO11_scer_hits.faa
