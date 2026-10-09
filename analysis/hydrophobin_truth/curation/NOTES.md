# Evidence for the T4 curation of the 61 unlabelled hydrophobin calls (2026-10-08)

This directory holds evidence only. No row is called hydrophobin or not hydrophobin here. The rubric is `../curation_rubric.md`. Nothing was committed to git.

## Files

| File | Content |
|---|---|
| `curation_evidence.tsv` | One row per protein (61). Lines A to E, with details. UniProt and FungiDB mapping. |
| `literature.tsv` | Papers per species (21 rows). PMID, DOI, what was read, gene IDs as written, protein-level evidence. |
| `literature_gene_map.tsv` | 98 rows. Each paper gene ID, its UniProt sequence, and the best BLASTP hit in the species proteome. |
| `scripts/` | `extract61.py`, `line_a.py`, `assemble.py`, `litmap.py`, `lit_genes.tsv` (the paper gene list used by `litmap.py`). |

Downloads (UniProt FASTA and JSON, PubMed full-text dumps) are in `/scratch/jstajich/29395251/hydro_curation/dl/`. That directory is node-local scratch of job 29395251. It is not kept.

## Method per line

- **Line A.** `line_a.py`. All Cys positions are found. For each window of k consecutive Cys (k = 8, 9, 10), every choice of 8 Cys that keeps the first and last Cys of the window is tested. A choice qualifies if the C2-C3 gap is 0 and the C6-C7 gap is 0. Line A = yes if a choice qualifies, R0 is `called` (from `evidence_sheets.tsv`) and length <= 300. `A_detail` gives the first qualifying window (C1-C8 positions and gaps), the number of qualifying choices, or all Cys gaps when none qualifies.
  - Test: P28346 (RodA, E. nidulans), Q04571 (RodL, N. crassa), P41746 (RodA, A. fumigatus), P52754 (HFB1, T. reesei), P16933 (SC3) all give yes, with C2C3 = 0 and C6C7 = 0 at the expected Cys.
  - Test on all 174 known hydrophobins (`positives_t2t3.faa`, R0 assumed called): motif found in 166. 161 are yes. 5 have the motif but are longer than 300 aa. 8 have no qualifying window (for example Q4WE22 RodD, P52755, J5JTB1 Hyd2B, W8NQ71).
  - Caution: with k = 9 or 10, a protein with many Cys can qualify by combination. `A_detail` says which window size qualified. Example: F00FD2C2_001945-T1 (19 Cys) qualifies only with k = 10.
- **Line B.** `strict_families` non-empty in `evidence_sheets.tsv`.
- **Line C.** Yes only if a curated record or a paper names this gene a hydrophobin and the gene maps by sequence. Sources checked: (1) the mapped UniProt entry: counts if reviewed, or if the name has PubMed evidence, or if a comment has ECO:0000269; (2) the FungiDB release-31 `gene_product` of the mapped gene (only 4 species are in FungiDB-31); (3) the paper genes in `literature_gene_map.tsv`. The mapping rule is >= 95% identity and >= 90% coverage of **both** sequences (coverage merged over HSPs). A match with >= 95% identity and >= 90% coverage of only one sequence is written `PARTIAL` and does not count.
- **Line D.** BLASTP (E <= 1e-3) of each protein against the 174 known hydrophobins, best hit by bitscore. Reverse: BLASTP of the 174 against the protein's own proteome, best hit. Yes if the pair is reciprocal and identity >= 40% (top HSP) and coverage >= 70% of both (merged HSPs). The clause "one-to-one ortholog in a species where hydrophobins are annotated" was **not checked**.
- **Line E.** (1) Swiss-Prot 2023_03, BLASTP E <= 1e-3, best hit by bitscore. Strong = identity >= 40% and query coverage >= 70%. The hit counts as non-hydrophobin if its accession is not one of the 174 and its title has no "hydrophob" or "rodlet" and is not "Uncharacterized". (2) `hmmscan --cut_ga` against Pfam-A 38.2 (`/bigdata/operations/pkgadmin/srv/projects/db/pfam/current`). All domains and envelope coverage are listed. Hydrophobin-class models are the 7 in `data/sorting_hat/family_table.tsv` (PF01185, PF06766, PF28987, PF22354, PF29785, PF29802, PF29465). "Accounts for most of the sequence" is taken as non-hydrophobin domain envelopes covering > 50% of the full length (signal peptide not removed). This threshold is my reading of the rubric word "most". E = yes if (1) or (2) holds.
- **UniProt mapping.** Afum_Af293 and Afum_A1163 proteome IDs are UniProt accessions, used as is. Others: BLASTP against UniProtKB 2026_03 entries downloaded by `organism_id` (Calb 237561, Cimm 246410, Bder 559297, Bbas 655819, Fful 5499, Fgra 229533, Pexp 27334), or by `taxonomy_id` (P. ostreatus 5322, because PC9 1137139 has 1 entry; A. fumigatus 746128 for W72310, which has only 577 own entries). So W72310 and PC9 proteins can map to an entry of another strain. The `OX=` of the hit is written in `uniprot_mapping`.

## Counts (61 proteins)

| Item | Count |
|---|---|
| Line A yes | 24 |
| Line B yes | 21 |
| Line C yes | 0 |
| Line C PARTIAL (not counted) | 2: Afum_W72310 KAK9636619.1 (54 aa; 100% identity to RodD AFUA_5G01490 / AFUB_050030 of Valsecchi 2018, but only 27% of RodD is covered); Pexp_MD-8 FA144B67_000394-T1 (99.3% identity to HfbA A0A0A2K7M3 of Luciano-Rosario 2022, 90% of FA144B67_000394-T1 but 65% of HfbA covered) |
| Line D yes | 6 |
| Line E yes | 8 (6 by a strong Swiss-Prot non-hydrophobin hit, 2 by Pfam coverage only: F0349401_004011-T1 Thaumatin, F0349401_005634-T1 ARB_05566_N) |
| Line A yes and line E yes | 1: Bbas F1BB8A46_003331-T1 (Swiss-Prot AMP1_MIRJA antimicrobial peptide, 52% identity over 87%) |
| Named hydrophobin in a paper or curated record, by the mapping rule | 0 |
| Mapped to a UniProt entry (rule met) | 49 |
| UniProt PARTIAL | 5 (KAK9636619.1, F00FD2C2_006466-T1, FEBB419A_012600-T1, F0349401_006476-T1, FA144B67_000394-T1) |
| No UniProt entry at >= 95% identity | 7 |
| Mapped UniProt entry named "Hydrophobin..." (all unreviewed, name by RuleBase, ProtNLM or EMBL submitter; do not count for C) | 12 |
| Mapped UniProt entry that is Swiss-Prot | 1 (C6_00820W_A = SUN41 Q59NP5) |
| FungiDB-31 product mapped | 19; none says hydrophobin |
| Another (non-hydrophobin) Pfam-A domain at GA | 16 |
| Hit to the 174 known hydrophobins at E <= 1e-3 | 24 |

## Observations to check (facts, not decisions)

1. `evidence_sheets.tsv` column `blast_known_hydrophobins` says "no hit at 1e-3" for all 61 rows. My BLASTP of the same sequences against `positives_t2t3.faa` (makeblastdb, E <= 1e-3) gives hits for 24 rows. I did not find the cause in `evidence_sheets.py`. The column should be re-checked.
2. Several rows look like gene-model differences, not new genes: FA144B67_000394-T1 vs HfbA (UniProt HfbA is 227 aa, the proteome model 164 aa); KAK9636619.1 is a 54 aa fragment of RodD; FEBB419A_012600-T1 (358 aa) contains UniProt A0A9Q8PL52 "Hydrophobin" (99 aa) at 97% over 28% of the proteome protein; F0349401_006476-T1 (96 aa) is 100% identical to 35% of UniProt A0A1C3YND2 "Cutinase" (FGRAMPH1_01G15975).
3. Near misses to named paper genes (below the 95% rule): W72310 KAK9641117.1 vs RodE 84.6% (100% of KAK9641117.1, 81% of RodE); Bbas F1BB8A46_000504-T1 vs Hyd2C J4WMI6 75.4% over 100%/100%; Post F4DD442B_010567-T1 vs Hydph8 (strain PC15) 92.8% over 100%/100%; Post F4DD442B_008436-T1 vs Hydph18 (PC15) 98.9% over 80%/64%. The P. ostreatus UniProt entries that cite Xu 2021 are mostly strain PC15. The proteome is PC9. Strain differences can push true orthologs below 95%.
4. Paper genes with no hit in the proteome: RodG in Afum_A1163, HfbD in Pexp_MD-8, Hyd1F in Bbas, Hydph11 in Post_PC9. HCf-3 and HCf-5 have no >= 95% hit in Fful_Race5. These genes may be missing from the proteome annotation. I did not check the genomes.
5. No paper was found that names a hydrophobin gene in C. immitis, B. dermatitidis, C. albicans or S. cerevisiae (searches listed in `literature.tsv`).

## What failed or was not done

- NCBI E-utilities did not respond from this node (empty reply). FGSG, PEXP (Pe21) and other paper IDs were resolved through UniProt instead.
- FungiDB release-31 has only A. fumigatus Af293, C. albicans SC5314, C. immitis RS and F. graminearum PH-1 of these species. FungiDB-31 F. graminearum uses FGRAMPH1_ IDs, not FGSG_.
- Full text was read only where PMC had it: Quarantin 2019, Valsecchi 2018, Luciano-Rosario 2022, Izumi 2026. All others are abstract only (stated per row in `literature.tsv`). Supplementary tables were not opened in this task.
- In the PMC text of Valsecchi 2018 and Luciano-Rosario 2022 the italic gene names are lost. Names were matched to IDs by the UniProt entries that cite the paper (Valsecchi) or by the stated order (Luciano-Rosario, HfbA to HfbG). The Luciano-Rosario Pe21 IDs have no sequence in the text; the Swiss-Prot sequences used are from strain MD-8.
- The one-to-one ortholog clause of line D was not checked.
- Line C for "orthologous gene product" in another species was not checked beyond line D.

## Commands

```bash
source /etc/profile.d/modules.sh; module load hmmer/3.4 ncbi-blast/2.14.0+
W=/scratch/jstajich/29395251/hydro_curation; T=analysis/hydrophobin_truth
/usr/bin/python3.12 -I scripts/extract61.py $T $W/work                       # q61.faa from the proteome FASTAs in run_list.tsv
/usr/bin/python3.12 -I scripts/line_a.py $W/work/test_t2.faa -               # test on P28346, Q04571, P41746, P52754, P16933
/usr/bin/python3.12 -I scripts/line_a.py $T/positives_t2t3.faa -             # all 174
/usr/bin/python3.12 -I scripts/line_a.py $W/work/q61.faa $T/evidence_sheets.tsv > $W/work/lineA.tsv
hmmscan --cut_ga --cpu 1 --domtblout pfam61.domtbl -o /dev/null /bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm q61.faa
blastp -query q61.faa -db /bigdata/operations/pkgadmin/srv/projects/db/Swissprot/2023_03/uniprot_sprot.fasta -evalue 1e-3 -max_target_seqs 25 \
  -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle' > sp61.tsv
makeblastdb -in $T/positives_t2t3.faa -dbtype prot -out pos174
blastp -query q61.faa -db pos174 -evalue 1e-3 -max_target_seqs 174 -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' > q61_vs_pos.tsv
# per proteome P: makeblastdb of the proteome; blastp positives_t2t3.faa vs P (-evalue 1e-3 -max_target_seqs 5) > pos_vs_P.tsv
# UniProt: curl "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:<taxon>&format=fasta&compressed=true" (taxonomy_id:5322 and taxonomy_id:746128 as above)
# per proteome: blastp q61 subset vs UniProt set (-evalue 1e-5 -max_target_seqs 10) > up/hits_P.tsv ; vs FungiDB-31 *_AnnotatedProteins.fasta > fdb/hits_P.tsv
# UniProt JSON per hit: curl https://rest.uniprot.org/uniprotkb/<acc>.json
# Papers to UniProt: curl "https://rest.uniprot.org/uniprotkb/stream?query=lit_pubmed:<pmid>&fields=accession,id,reviewed,gene_names,organism_id,length,protein_name&format=tsv"
/usr/bin/python3.12 -I scripts/litmap.py scripts/lit_genes.tsv $W/dl/litup/lit_seqs.faa $W/work/db $T/evidence_sheets.tsv curation/literature_gene_map.tsv
/usr/bin/python3.12 -I scripts/assemble.py $T $W/work $W/dl/upjson curation/literature_gene_map.tsv curation/curation_evidence.tsv
```

PubMed searches (PubMed MCP tool, 2026-10-08): `hydrophobin[Title/Abstract] AND (Fusarium graminearum OR Gibberella zeae)` (8 hits); `... AND Beauveria bassiana` (28); `... AND (Cladosporium fulvum OR Fulvia fulva OR Passalora fulva)` (7); `... AND Pleurotus ostreatus` (38); `... AND Penicillium expansum` (2); `... AND Aspergillus fumigatus` (42); `... AND (Coccidioides OR Blastomyces OR Candida albicans OR Saccharomyces cerevisiae)` (15); `hydrophobin AND (Coccidioides OR Blastomyces OR Ajellomyces)` (3); `hydrophobin[Title] AND Aspergillus fumigatus[Title/Abstract]` (10). Not every hit was read. Papers were chosen by title for gene lists or gene characterisation.
