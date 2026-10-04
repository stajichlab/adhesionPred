# Spherule follow-up: specificity, gene models, Cys domains and signal peptides

*2026-10-03. Branch `cocci-extreme-genes`. Follows `2026-10-03-cocci-spherule-surface-table.md` (PR #56).
Scripts: `analysis/cocci_spherule/10_` to `13_`. Counts are computed from the output files.*

## 1. What I checked and what came out

The first report left a set of 45 genes to check: 26 specific
extreme genes without a signal peptide (set A) and 25 Cys-rich spherule-up genes without a signal peptide that pass
the specificity flag (set B). 6 genes are in both. All 45 pass the ranking's specificity flag
(45 of 45).

1. **The ranking's specificity flag means "no hit at 30% identity", not "no homolog".** The flag
   (`cocci_antigens/01_build_inputs.sh`) searched one proteome per genus (*Histoplasma*, *Blastomyces*,
   *Paracoccidioides*, *A. fumigatus*) and human, and kept hits with at least 30% identity and 50% coverage
   of query and target. I re-tested the 45 genes with diamond (very sensitive) against the same proteomes.
   44 of 45 still have no hit that meets those criteria. CIMG_13156
   meets them, so the flag missed it. Only 8 of 45 genes have any diamond hit (E-value at most 1e-3)
   in those five proteomes, and these are distant or partial (section 2b).
   The phmmer search (section 2) finds hits in more proteomes: 26 of 45 genes have
   a hit outside *Coccidioides* at E-value at most 1e-5, 16 of them in a dimorphic relative or an
   outgroup (class C1). I read these as family-level or domain-level matches. They are not orthologs.
   Example: CIMG_04662 is a P-type cation ATPase. Its best hits in human and in *Histoplasma* are other
   P-type ATPases at 26% to 27% identity over about 90% of the length, below the flag's cut-off. phmmer finds it
   in 16 dimorphic-relative proteomes
   for the same reason. Whether it has an ortholog elsewhere needs a reciprocal-best-hit or tree analysis.
2. **19 genes have no phmmer hit outside *Coccidioides*** (class C3). This is the strongest evidence of restriction here: no match of any kind at E-value 1e-5 in 87 proteomes. Only 3 of them are
   supported by a second annotation of the same genome and are not rare in the pangenome
   (CIMG_13082, CIMG_13230, CIMG_13657).
3. **Gene models.** Only 5 of 45 proteins have a close match (at least 95% identity and 90% coverage on
   both sides) in the pangenome RS annotation. 11 match partly or with lower identity. 29 have
   no match at all (diamond ultra-sensitive, E-value at most 1e-3). 21 are absent from the Fungi5k
   annotation of RS. 13 are in at most 1% of the 488 proteomes, so they are essentially RS-only.
   All 45 genes are up in spherules in the RNA-seq, so reads map to the RefSeq transcript sequences. That does
   not show that the RefSeq gene models are right or that the genes encode proteins.
4. **Phobius calls more signal peptides than SignalP.** Among the 564 spherule-up genes, Phobius calls
   42 and SignalP calls 23. All 23 SignalP calls are also Phobius calls.
   19 genes are Phobius-only, 15 of them with no TM helix. Phobius is usually
   less specific than SignalP, so read the Phobius-only genes as candidates.
5. **No CFEM domain.** None of the 564 spherule-up proteins has a CFEM hit (PF05730, Pfam gathering
   threshold). 6 of the 45 follow-up genes have any Pfam domain.

## 2. Methods

| Check | Tool | Against |
|---|---|---|
| Orthologs | phmmer (HMMER 3.4), E-value at most 1e-3 searched, 1e-5 used for "strict" | 69 Onygenales proteomes of Fungi5k, the two *C. posadasii* references (Silveira, C735-TKO), and one proteome for each of 16 outgroup genera (`search_proteomes.tsv`). The Fungi5k annotation of RS itself is excluded from the hit counts and used only for the model check |
| Gene model | diamond blastp ultra-sensitive, E-value at most 1e-3 | pangenome RS proteome (`C3CED3A_*`) |
| Domains | hmmscan, Pfam-A, gathering thresholds | all 564 spherule-up proteins |
| Signal peptide | Phobius 1.01 | all 564 spherule-up proteins |

Class rules (E-value only, so family-level and domain-level matches count): C1 = strict hit in a dimorphic relative (*Histoplasma*, *Blastomyces*, *Paracoccidioides*,
*Emergomyces*, *Emmonsia*, *Ajellomyces*) or an outgroup. C2 = strict hit only in other Onygenales. C3 = no
strict hit outside *Coccidioides*. Model agrees = 95% identity and 90% coverage of query and subject.

## 2b. The specificity flag, re-tested with its own criteria

Same proteomes as the flag, diamond very-sensitive, hits kept at E-value 1e-3. The flag's criteria:
identity at least 30%, coverage at least 50% of query and of target. The table lists every follow-up gene
with a hit, with its best three hits.

| Gene | Best hits: group identity% / query cov / target cov | Meets the flag's criteria |
|---|---|---|
| CIMG_04662 | human 27% / 0.90 / 0.91; human 27% / 0.90 / 0.87; human 26% / 0.91 / 0.90 | no |
| CIMG_06250 | human 26% / 0.51 / 0.44 | no |
| CIMG_11763 | Aspergillus_fumigatus 30% / 0.84 / 0.77; Paracoccidioides 28% / 0.33 / 0.23 | no |
| CIMG_11958 | Blastomyces 45% / 0.29 / 0.31; Blastomyces 38% / 0.35 / 0.16 | no |
| CIMG_12971 | Histoplasma 28% / 0.55 / 0.88; Histoplasma 35% / 0.41 / 0.62; Histoplasma 38% / 0.36 / 0.96 | no |
| CIMG_13156 | Histoplasma 44% / 0.38 / 0.97; Blastomyces 42% / 0.43 / 0.62; Histoplasma 40% / 0.48 / 0.69 | yes |
| CIMG_13214 | Histoplasma 43% / 0.78 / 0.12; Blastomyces 32% / 0.64 / 0.39 | no |
| CIMG_13220 | Paracoccidioides 26% / 0.80 / 0.24; Histoplasma 27% / 0.80 / 0.24 | no |

`flag_check_hits.tsv` has all 66 hits. CIMG_11763 has a hit at about 30% identity and 77% to 84%
coverage in *A. fumigatus*. It is at the boundary, and I did not check the exact identity.

## 3. The 45 follow-up genes

Hits column: number of proteomes with a strict hit in *Coccidioides* (other than RS), dimorphic relatives,
other Onygenales and outgroups.

| Gene | Class | Set | aa | log2FC | Cys | Prevalence | Model vs pangenome RS | In Fungi5k RS | Hits cocc/dimorph/other ony/outgroup | Pfam | Phobius SP |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CIMG_04662 | C1 | A | 929 | 8.85 | 0.018 | 0.99 | agree | yes | 2/16/52/16 | Cation_ATPase(PF13246);Cation_ATPase_C(P | no |
| CIMG_12615 | C1 | A | 301 | 7.19 | 0.013 | 0.0615 | no hit | yes | 0/1/8/0 | - | no |
| CIMG_08364 | C1 | A | 179 | 7.12 | 0.022 | 0.352 | partial (100.0%, cov 1.00) | yes | 0/2/3/1 | - | no |
| CIMG_06250 | C1 | A | 320 | 6.53 | 0.019 | 1 | agree | yes | 2/0/13/3 | Methyltransf_21(PF05050);Methyltransf_22 | yes |
| CIMG_13125 | C1 | A | 139 | 6.28 | 0.029 | 0.002 | partial (74.4%, cov 0.56) | yes | 2/1/6/0 | - | no |
| CIMG_13214 | C1 | AB | 133 | 6.26 | 0.045 | 0.844 | partial (49.5%, cov 0.80) | yes | 0/7/8/1 | - | no |
| CIMG_07389 | C1 | B | 150 | 5.07 | 0.053 | 0.199 | agree | no | 0/1/5/0 | - | no |
| CIMG_13156 | C1 | B | 259 | 4.76 | 0.058 | 0.266 | partial (37.6%, cov 0.45) | yes | 2/16/51/16 | Histone(PF00125) | no |
| CIMG_11958 | C1 | B | 173 | 4.55 | 0.052 | 0.314 | partial (96.3%, cov 0.47) | no | 0/4/1/0 | - | no |
| CIMG_13539 | C1 | B | 116 | 4.16 | 0.06 | 0.002 | no hit | yes | 0/2/4/0 | - | no |
| CIMG_13220 | C1 | B | 198 | 3.79 | 0.056 | 0.295 | partial (100.0%, cov 0.54) | yes | 2/12/22/0 | F-box(PF00646) | no |
| CIMG_12971 | C1 | B | 272 | 3.57 | 0.04 | 0.266 | partial (37.9%, cov 0.32) | yes | 2/16/51/15 | Histone(PF00125) | no |
| CIMG_09338 | C1 | B | 143 | 3.45 | 0.063 | 0.859 | partial (100.0%, cov 0.75) | no | 2/0/3/2 | - | no |
| CIMG_13151 | C1 | B | 280 | 3.41 | 0.057 | 0.002 | no hit | yes | 0/2/5/0 | - | no |
| CIMG_11763 | C1 | B | 357 | 2.45 | 0.045 | 0.709 | partial (25.4%, cov 0.40) | yes | 2/4/16/4 | F-box(PF00646) | no |
| CIMG_13457 | C1 | B | 95 | 2.25 | 0.042 | 0.221 | no hit | no | 0/4/0/0 | - | no |
| CIMG_10877 | C2 | A | 163 | 7.22 | 0.0061 | 0.002 | no hit | no | 0/0/1/0 | - | no |
| CIMG_13332 | C2 | A | 263 | 7.1 | 0.0076 | 0.0615 | no hit | yes | 0/0/5/0 | - | no |
| CIMG_12573 | C2 | A | 238 | 6.86 | 0.0084 | 0.0615 | no hit | yes | 0/0/5/0 | - | no |
| CIMG_12595 | C2 | AB | 100 | 6.55 | 0.05 | 0.002 | no hit | no | 0/0/2/0 | - | yes |
| CIMG_13004 | C2 | A | 289 | 6.48 | 0.01 | 0.0615 | no hit | yes | 0/0/2/0 | - | no |
| CIMG_12994 | C2 | A | 130 | 6.47 | 0 | 0.002 | no hit | yes | 0/0/1/0 | - | no |
| CIMG_13499 | C2 | A | 205 | 6.37 | 0.019 | 0.0061 | no hit | yes | 0/0/1/0 | - | no |
| CIMG_12854 | C2 | A | 253 | 6.02 | 0.028 | 0.998 | partial (86.6%, cov 0.89) | yes | 2/0/4/0 | - | no |
| CIMG_12625 | C2 | A | 121 | 6 | 0.017 | 0.0615 | no hit | yes | 0/0/6/0 | - | no |
| CIMG_12995 | C2 | B | 116 | 5.31 | 0.043 | 0.266 | no hit | yes | 1/0/1/0 | - | no |
| CIMG_13657 | C3 | A | 162 | 8.26 | 0.012 | 0.0615 | no hit | yes | 0/0/0/0 | - | no |
| CIMG_12484 | C3 | A | 124 | 7.16 | 0.016 | 0.0779 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12808 | C3 | AB | 143 | 6.63 | 0.042 | 0.0143 | no hit | no | 0/0/0/0 | - | yes |
| CIMG_11522 | C3 | AB | 85 | 6.62 | 0.047 | 0.0492 | no hit | no | 0/0/0/0 | - | no |
| CIMG_13230 | C3 | A | 104 | 6.62 | 0.019 | 0.201 | agree | yes | 0/0/0/0 | - | no |
| CIMG_10007 | C3 | A | 125 | 6.59 | 0.008 | 0.17 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12953 | C3 | A | 170 | 6.4 | 0.012 | 0.002 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12716 | C3 | A | 77 | 6.38 | 0 | 0.002 | no hit | no | 0/0/0/0 | - | yes |
| CIMG_03730 | C3 | AB | 41 | 6.34 | 0.049 | 0.002 | no hit | no | 0/0/0/0 | - | no |
| CIMG_10013 | C3 | AB | 114 | 6.18 | 0.044 | 0.0369 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12654 | C3 | A | 92 | 6.05 | 0.033 | 0.002 | no hit | yes | 0/0/0/0 | - | no |
| CIMG_12913 | C3 | B | 75 | 4.54 | 0.04 | 0.0061 | no hit | yes | 0/0/0/0 | - | no |
| CIMG_13084 | C3 | B | 200 | 3.6 | 0.04 | 0.65 | partial (100.0%, cov 0.82) | no | 1/0/0/0 | - | no |
| CIMG_12945 | C3 | B | 130 | 2.91 | 0.061 | 0.0164 | no hit | no | 0/0/0/0 | - | no |
| CIMG_03569 | C3 | B | 100 | 2.48 | 0.04 | 0.111 | no hit | no | 0/0/0/0 | - | yes |
| CIMG_03444 | C3 | B | 133 | 2.27 | 0.045 | 0.15 | no hit | no | 0/0/0/0 | - | no |
| CIMG_11847 | C3 | B | 178 | 2.12 | 0.045 | 0.189 | no hit | no | 0/0/0/0 | - | no |
| CIMG_04176 | C3 | B | 159 | 2.07 | 0.05 | 0.0041 | no hit | no | 0/0/0/0 | - | no |
| CIMG_13082 | C3 | B | 111 | 2.02 | 0.045 | 0.0963 | agree | no | 0/0/0/0 | - | yes |

## 4. Class C3: no hit outside *Coccidioides* (19 genes)

These are the candidates for *Coccidioides*-restricted genes. Most lack support in a second annotation.

| Gene | Class | Set | aa | log2FC | Cys | Prevalence | Model vs pangenome RS | In Fungi5k RS | Hits cocc/dimorph/other ony/outgroup | Pfam | Phobius SP |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CIMG_13657 | C3 | A | 162 | 8.26 | 0.012 | 0.0615 | no hit | yes | 0/0/0/0 | - | no |
| CIMG_12484 | C3 | A | 124 | 7.16 | 0.016 | 0.0779 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12808 | C3 | AB | 143 | 6.63 | 0.042 | 0.0143 | no hit | no | 0/0/0/0 | - | yes |
| CIMG_11522 | C3 | AB | 85 | 6.62 | 0.047 | 0.0492 | no hit | no | 0/0/0/0 | - | no |
| CIMG_13230 | C3 | A | 104 | 6.62 | 0.019 | 0.201 | agree | yes | 0/0/0/0 | - | no |
| CIMG_10007 | C3 | A | 125 | 6.59 | 0.008 | 0.17 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12953 | C3 | A | 170 | 6.4 | 0.012 | 0.002 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12716 | C3 | A | 77 | 6.38 | 0 | 0.002 | no hit | no | 0/0/0/0 | - | yes |
| CIMG_03730 | C3 | AB | 41 | 6.34 | 0.049 | 0.002 | no hit | no | 0/0/0/0 | - | no |
| CIMG_10013 | C3 | AB | 114 | 6.18 | 0.044 | 0.0369 | no hit | no | 0/0/0/0 | - | no |
| CIMG_12654 | C3 | A | 92 | 6.05 | 0.033 | 0.002 | no hit | yes | 0/0/0/0 | - | no |
| CIMG_12913 | C3 | B | 75 | 4.54 | 0.04 | 0.0061 | no hit | yes | 0/0/0/0 | - | no |
| CIMG_13084 | C3 | B | 200 | 3.6 | 0.04 | 0.65 | partial (100.0%, cov 0.82) | no | 1/0/0/0 | - | no |
| CIMG_12945 | C3 | B | 130 | 2.91 | 0.061 | 0.0164 | no hit | no | 0/0/0/0 | - | no |
| CIMG_03569 | C3 | B | 100 | 2.48 | 0.04 | 0.111 | no hit | no | 0/0/0/0 | - | yes |
| CIMG_03444 | C3 | B | 133 | 2.27 | 0.045 | 0.15 | no hit | no | 0/0/0/0 | - | no |
| CIMG_11847 | C3 | B | 178 | 2.12 | 0.045 | 0.189 | no hit | no | 0/0/0/0 | - | no |
| CIMG_04176 | C3 | B | 159 | 2.07 | 0.05 | 0.0041 | no hit | no | 0/0/0/0 | - | no |
| CIMG_13082 | C3 | B | 111 | 2.02 | 0.045 | 0.0963 | agree | no | 0/0/0/0 | - | yes |

Genes with a second-annotation match and prevalence above 1% are the best supported: CIMG_13082, CIMG_13230, CIMG_13657.
A phmmer miss at E-value 1e-5 does not prove absence. Short proteins (under 150 aa) have little power.

## 5. Phobius-only signal peptides among spherule-up genes (19)

| Gene | Product | aa | log2FC 48 h | SignalP SP prob | Phobius TM | Cys | Specific flag |
|---|---|---|---|---|---|---|---|
| CIMG_12808 | uncharacterized protein CIMG_12808 | 143 | 6.63 | 0.00011 | 0 | 0.042 | yes |
| CIMG_12595 | uncharacterized protein CIMG_12595 | 100 | 6.55 | 0.32 | 0 | 0.05 | yes |
| CIMG_06250 | uncharacterized protein CIMG_06250 | 320 | 6.53 | 0.27 | 0 | 0.019 | yes |
| CIMG_12716 | uncharacterized protein CIMG_12716 | 77 | 6.38 | 0.0056 | 0 | 0 | yes |
| CIMG_13336 | uncharacterized protein CIMG_13336 | 58 | 5.87 | 0.017 | 0 | 0.017 | yes |
| CIMG_05158 | uncharacterized protein CIMG_05158 | 168 | 5.67 | 0.031 | 0 | 0.012 | yes |
| CIMG_13542 | uncharacterized protein CIMG_13542 | 140 | 5.3 | 7.2e-05 | 0 | 0 | yes |
| CIMG_13129 | uncharacterized protein CIMG_13129 | 216 | 5.01 | 0 | 1 | 0.023 | yes |
| CIMG_11498 | uncharacterized protein CIMG_11498 | 243 | 3.63 | 1e-06 | 3 | 0.012 | yes |
| CIMG_08344 | uncharacterized protein CIMG_08344 | 230 | 3.5 | 0.0044 | 0 | 0.022 | yes |
| CIMG_09765 | oxidoreductase | 243 | 3.45 | 0 | 0 | 0 | no |
| CIMG_11666 | carboxypeptidase S1 | 651 | 3.39 | 0.29 | 1 | 0.015 | no |
| CIMG_12951 | uncharacterized protein CIMG_12951 | 164 | 3.26 | 0.021 | 1 | 0.012 | yes |
| CIMG_09292 | uncharacterized protein CIMG_09292 | 372 | 2.61 | 2e-06 | 0 | 0.011 | no |
| CIMG_03569 | uncharacterized protein CIMG_03569 | 100 | 2.48 | 0 | 0 | 0.04 | yes |
| CIMG_07302 | uncharacterized protein CIMG_07302 | 370 | 2.33 | 0 | 0 | 0.027 | yes |
| CIMG_12424 | uncharacterized protein CIMG_12424 | 171 | 2.15 | 0.002 | 0 | 0.023 | yes |
| CIMG_13082 | uncharacterized protein CIMG_13082 | 111 | 2.02 | 0.00077 | 0 | 0.045 | yes |
| CIMG_02660 | flavin-containing amine oxidasedeh | 544 | 2.01 | 4.6e-05 | 0 | 0.013 | no |

## 6. Limits

- The orthology search is sensitive but not complete. phmmer used single-sequence queries. A profile search
  (jackhmmer, or an HMM built from the hits) would find more.
- The outgroup set is 16 genera, one proteome each. It does not span the fungal tree.
- Absence from the pangenome RS annotation or from Fungi5k RS means the other pipelines did not predict the
  gene. It does not show that the RefSeq model is wrong, and the reverse is also open.
- The flag uses one proteome per confounder genus and a 30% identity cut-off. That is a cross-reactivity
  proxy for serology. It is not a test for homology or orthology, and it should not be read as "genus-specific".
  The phmmer classes C1 and C2 are not orthology calls either.
- Phobius-only calls can be false positives. TM-containing proteins can also be signal anchors.
- Two replicates, one strain (first report, section 7).
- No experiment tests any surface localisation here.

## 7. Next

- Map the 19 C3 genes to the RS genome and check exon structure against the RNA-seq reads
  (the kallisto run does not show gene structure).
- Decide whether the table's "specific" column should be renamed (for example "no 30% identity hit in four
  confounder fungi or human") so it is not read as "no homolog". That change belongs in the first report.
- Decide whether a hand-checked gene list is worth sending for a proteomics or localisation test.
