# Family-matched negatives: counts

Date: 2026-10-10. UniProt release 2026_03. Per-family cap 10. Sampling seed 20261010.

## Totals

| Item | Count |
|---|---|
| Positives in positives.tsv | 200 |
| Positives with UniProt accession (queried) | 174 |
| UniProt positives with at least one Pfam domain | 152 |
| Positives skipped (NCBI protein IDs, no UniProt mapping) | 26 |
| Distinct Pfam families used by positives | 76 |
| Families with at least one new negative | 67 |
| New negatives (family_negatives.faa) | 516 |
| Distinct new negative accessions | 516 |
| New negatives in candidates.tsv | 0 (excluded from pools) |
| Pool rows (all families, with duplicates) that are candidates | 542 |
| Pool rows that are positives | 19 |
| Pool rows (all families, with duplicates) | 2543 |

## MMseqs2 clusters (easy-cluster, --min-seq-id 0.3 -c 0.5 --cov-mode 0)

| Item | Count |
|---|---|
| Input sequences (200 positives + 516 new negatives) | 716 |
| Clusters containing positives | 82 |
| Clusters containing new negatives | 157 |
| Clusters containing both | 49 |
| Positive-only clusters | 33 |
| Negative-only clusters | 108 |

## Per family

Pool = reviewed fungal Swiss-Prot entries with the Pfam, after removing the KW-0020 and cc_allergen:allergic hits in the query; The pool column subtracts candidate and positive hits. Negatives already picked for an earlier family (families in sorted order) are removed before sampling, so New negatives can be below the cap.

| Pfam | Positives with this Pfam | Query hits | Pool (not positive/candidate) | Hits in candidates.tsv | Hits that are positives | New negatives |
|---|---|---|---|---|---|---|
| PF00010 | 1 | 21 | 12 | 9 | 0 | 10 |
| PF00012 | 4 | 70 | 46 | 24 | 0 | 10 |
| PF00026 | 3 | 111 | 66 | 43 | 2 | 10 |
| PF00034 | 1 | 23 | 17 | 5 | 1 | 10 |
| PF00056 | 1 | 8 | 5 | 3 | 0 | 5 |
| PF00076 | 1 | 280 | 211 | 69 | 0 | 10 |
| PF00081 | 3 | 11 | 8 | 3 | 0 | 8 |
| PF00082 | 24 | 151 | 132 | 17 | 2 | 10 |
| PF00085 | 8 | 32 | 20 | 12 | 0 | 10 |
| PF00107 | 2 | 131 | 109 | 22 | 0 | 10 |
| PF00113 | 8 | 25 | 19 | 4 | 2 | 10 |
| PF00128 | 2 | 41 | 25 | 14 | 2 | 10 |
| PF00160 | 5 | 114 | 100 | 14 | 0 | 10 |
| PF00171 | 2 | 36 | 27 | 9 | 0 | 10 |
| PF00183 | 1 | 9 | 5 | 4 | 0 | 5 |
| PF00188 | 1 | 15 | 11 | 4 | 0 | 10 |
| PF00199 | 1 | 30 | 20 | 10 | 0 | 10 |
| PF00258 | 1 | 28 | 22 | 6 | 0 | 10 |
| PF00262 | 2 | 6 | 4 | 2 | 0 | 4 |
| PF00297 | 1 | 13 | 9 | 4 | 0 | 9 |
| PF00326 | 3 | 42 | 38 | 4 | 0 | 10 |
| PF00328 | 2 | 25 | 8 | 16 | 1 | 8 |
| PF00428 | 8 | 22 | 14 | 8 | 0 | 10 |
| PF00544 | 1 | 46 | 43 | 3 | 0 | 10 |
| PF00545 | 2 | 20 | 16 | 4 | 0 | 10 |
| PF00686 | 2 | 11 | 4 | 7 | 0 | 4 |
| PF00722 | 2 | 17 | 8 | 9 | 0 | 8 |
| PF00723 | 2 | 17 | 14 | 3 | 0 | 10 |
| PF00728 | 2 | 5 | 3 | 2 | 0 | 3 |
| PF00732 | 2 | 46 | 25 | 21 | 0 | 10 |
| PF00736 | 1 | 3 | 2 | 1 | 0 | 2 |
| PF00838 | 2 | 15 | 13 | 2 | 0 | 10 |
| PF00923 | 3 | 5 | 3 | 2 | 0 | 3 |
| PF00933 | 2 | 73 | 64 | 8 | 1 | 10 |
| PF01116 | 1 | 8 | 6 | 2 | 0 | 6 |
| PF01915 | 2 | 74 | 65 | 8 | 1 | 10 |
| PF01965 | 1 | 10 | 5 | 5 | 0 | 5 |
| PF02128 | 1 | 68 | 67 | 1 | 0 | 10 |
| PF02136 | 2 | 13 | 12 | 1 | 0 | 10 |
| PF02777 | 3 | 16 | 12 | 4 | 0 | 4 |
| PF02826 | 1 | 25 | 18 | 7 | 0 | 10 |
| PF02866 | 1 | 8 | 5 | 3 | 0 | 0 |
| PF03358 | 1 | 12 | 9 | 3 | 0 | 9 |
| PF03952 | 8 | 25 | 19 | 4 | 2 | 9 |
| PF04681 | 1 | 4 | 4 | 0 | 0 | 4 |
| PF05199 | 2 | 52 | 26 | 26 | 0 | 10 |
| PF05922 | 20 | 130 | 118 | 10 | 2 | 10 |
| PF06628 | 1 | 19 | 10 | 9 | 0 | 5 |
| PF07249 | 2 | 16 | 13 | 3 | 0 | 10 |
| PF07335 | 2 | 6 | 6 | 0 | 0 | 6 |
| PF07504 | 1 | 64 | 63 | 1 | 0 | 10 |
| PF08240 | 2 | 168 | 144 | 24 | 0 | 10 |
| PF08534 | 7 | 5 | 2 | 3 | 0 | 2 |
| PF09260 | 2 | 16 | 7 | 7 | 2 | 5 |
| PF10287 | 1 | 5 | 2 | 3 | 0 | 2 |
| PF10290 | 1 | 5 | 2 | 3 | 0 | 0 |
| PF10587 | 1 | 3 | 2 | 1 | 0 | 0 |
| PF12296 | 1 | 3 | 3 | 0 | 0 | 3 |
| PF12697 | 1 | 34 | 32 | 2 | 0 | 10 |
| PF13409 | 1 | 14 | 11 | 3 | 0 | 10 |
| PF13561 | 2 | 88 | 73 | 15 | 0 | 10 |
| PF13589 | 1 | 11 | 5 | 6 | 0 | 3 |
| PF13848 | 2 | 7 | 5 | 2 | 0 | 3 |
| PF13933 | 1 | 4 | 3 | 1 | 0 | 3 |
| PF14310 | 1 | 65 | 58 | 6 | 1 | 10 |
| PF14845 | 1 | 5 | 3 | 2 | 0 | 0 |
| PF16541 | 13 | 2 | 1 | 1 | 0 | 1 |
| PF18011 | 1 | 8 | 5 | 3 | 0 | 0 |
| PF22354 | 1 | 11 | 10 | 1 | 0 | 10 |
| PF22701 | 1 | 2 | 2 | 0 | 0 | 2 |
| PF25312 | 1 | 0 | 0 | 0 | 0 | 0 |
| PF25484 | 1 | 0 | 0 | 0 | 0 | 0 |
| PF27508 | 1 | 0 | 0 | 0 | 0 | 0 |
| PF29664 | 1 | 30 | 26 | 4 | 0 | 10 |
| PF30672 | 3 | 5 | 5 | 0 | 0 | 5 |
| PF30996 | 1 | 0 | 0 | 0 | 0 | 0 |

Pfam names: PF00010 Helix-loop-helix DNA-binding domain; PF00012 Hsp70 protein; PF00026 Eukaryotic aspartyl protease; PF00034 Cytochrome c; PF00056 lactate/malate dehydrogenase, NAD binding domain; PF00076 RNA recognition motif; PF00081 Iron/manganese superoxide dismutases, alpha-hairpin domain; PF00082 Subtilase family; PF00085 Thioredoxin; PF00107 Zinc-binding dehydrogenase; PF00113 Enolase, C-terminal TIM barrel domain; PF00128 Alpha amylase, catalytic domain; PF00160 Cyclophilin type peptidyl-prolyl cis-trans isomerase/CLD; PF00171 Aldehyde dehydrogenase family; PF00183 Hsp90 protein; PF00188 Cysteine-rich secretory protein family; PF00199 Catalase; PF00258 Flavodoxin; PF00262 Calreticulin family; PF00297 Ribosomal protein L3; PF00326 Prolyl oligopeptidase family; PF00328 Histidine phosphatase superfamily (branch 2); PF00428 60s Acidic ribosomal protein; PF00544 Pectate lyase; PF00545 ribonuclease; PF00686 Starch binding domain; PF00722 Glycosyl hydrolases family 16; PF00723 Glycosyl hydrolases family 15; PF00728 Glycosyl hydrolase family 20, catalytic domain; PF00732 GMC oxidoreductase; PF00736 EF-1 guanine nucleotide exchange domain; PF00838 Translationally controlled tumour protein; PF00923 Transaldolase/Fructose-6-phosphate aldolase; PF00933 Glycosyl hydrolase family 3 N terminal domain; PF01116 Fructose-bisphosphate aldolase class-II; PF01915 Glycosyl hydrolase family 3 C-terminal domain; PF01965 DJ-1/PfpI family; PF02128 Fungalysin metallopeptidase (M36); PF02136 Nuclear transport factor 2 (NTF2) domain; PF02777 Iron/manganese superoxide dismutases, C-terminal domain; PF02826 D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain; PF03358 NADPH-dependent FMN reductase; PF03952 Enolase, N-terminal domain; PF04681 Blastomyces yeast-phase-specific protein; PF05199 GMC oxidoreductase; PF05922 Peptidase inhibitor I9; PF06628 Catalase-related immune-responsive; PF07249 Cerato-platanin; PF07335 Fungal chitosanase of glycosyl hydrolase group 75; PF07504 Fungalysin/Thermolysin Propeptide Motif; PF08240 Alcohol dehydrogenase GroES-like domain; PF08534 Redoxin; PF09260 Alpha-amylase, domain C; PF10287 Cell wall protein YJL171C/Tos1, C-terminal; PF12296 Hydrophobic surface binding protein A; PF12697 Alpha/beta hydrolase family; PF13409 Glutathione S-transferase, N-terminal domain; PF13561 Enoyl-(Acyl carrier protein) reductase; PF13589 Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase; PF13848 Thioredoxin-like domain; PF13933 Putative peptidase family; PF14310 Fibronectin type III-like domain; PF16541 Alternaria alternata allergen 1; PF22354 Hydrophobin; PF22701 Mal s 1 allergenic protein-like; PF29664 BtgE/BtgC glycoside hydrolase domain; PF30672 PhiA-like secreted protein.

## Skipped positives

26 positives have NCBI protein IDs (source ncbi_protein or compare, IDs like AAM77471.1 or XP_003030591.1). They were not assigned Pfam domains, so they add no family. The ID mapping step (scripts/step1b_idmap.py) used UniProt idmapping from EMBL-GenBank-DDBJ_CDS and RefSeq_Protein to UniProtKB. It mapped 0 of 26. The EMBL-to-UniParc step then returned HTTP 400, so no other route was tried. The 26 are kept in the clustering input.

