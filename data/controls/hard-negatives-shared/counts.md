# Counts: shared hard-negative set (task 07)

*Written 2026-10-08 by the task 07 agent (claude-opus-5-5). Source table: `controls.tsv`. Clusters:
MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5`, run on the hard negatives together with
104 adhesin positives (`adhesins.tsv` E1/E2 rows with an accession, plus the task 08 adhesin rows).
A cluster id counts once per group. Every row is a negative (`label=negative`).*

## Summary

| Item | Value |
|---|---|
| Rows (proteins) | 161 |
| Sequences in `sequences.faa` | 161 |
| Independent clusters (30% identity, 50% coverage) | 91 |
| N1 rows (characterised non-adhesive function) | 137 |
| N2 rows (adhesin-like, no adhesion evidence only) | 24 |
| Rows from non-yeast species (not Saccharomycotina) | 134 |
| Rows carried over from `adhesins.tsv` (`hard_negative`) | 40 of 42 (2 rejected) |
| Rows carried over from task 08 | 1 of 10 (9 rejected) |
| New rows | 120 |
| Rows that share a 30% cluster with an adhesin positive | 6 (flagged, kept): EXG1, exg1 (with gp43); bglH, bglI, bglJ, bglK (with two PA14 beta-glucosidases labelled E2 adhesin by domain) |
| Candidates rejected | 36 lines in `rejected.tsv` (some lines hold several proteins) |

## Size targets

Task target: at least 20 clusters per stratum and per clade that we report. COMMON-RULES rule 7:
the `estimated` status needs at least 20 negative clusters **in one species** for the call that
is measured.

| Target | Met? |
|---|---|
| 20 clusters per stratum | Met for `laccase_etc` (23) only. Not met for the other 8 strata (3 to 16 clusters). |
| 20 clusters per clade | Met for Eurotiales (27), Pezizomycotina other than Eurotiales and Onygenales (32), Saccharomycotina (23). Not met for Basidiomycota (14), Onygenales (15), Taphrinomycotina (5), other fungi (2). |
| 20 clusters per stratum and clade cell | Not met in any cell. The largest cell is `laccase_etc` x Onygenales (13). |
| 20 negative clusters in one species (rule 7, `estimated`) | Met for *S. cerevisiae* only (23 clusters, all strata pooled). Next: *A. fumigatus* 16, *V. dahliae* 12, *A. nidulans* 10. Pooling strata for one call is a choice for the owner. |

I did not pad any stratum. Families with many paralogs (GH16 Crh, GH72 Gel/Gas, GH18 chitinases,
AA1 laccases, CBM1 cellulases) add rows but few clusters.

### Per stratum

| Group | Rows | Sequences | Clusters | N1 rows | N2 rows | Clusters to 20 |
|---|---|---|---|---|---|---|
| domain_non_member | 23 | 23 | 16 | 12 | 11 | 4 |
| gpi_wall_enzyme | 21 | 21 | 8 | 21 | 0 | 12 |
| laccase_etc | 38 | 38 | 23 | 37 | 1 | 0 |
| mucin_sensor | 4 | 4 | 4 | 3 | 1 | 16 |
| repeat_non_adhesin | 12 | 12 | 7 | 11 | 1 | 13 |
| secretory_non_surface | 5 | 5 | 3 | 5 | 0 | 17 |
| st_linker_enzyme | 21 | 21 | 8 | 18 | 3 | 12 |
| wall_hydrolase | 28 | 28 | 15 | 24 | 4 | 5 |
| wall_structural_other | 9 | 9 | 7 | 6 | 3 | 13 |
| **all** | 161 | 161 | 91 | 137 | 24 | |

### Clusters per stratum and clade (rows in brackets)

| Stratum | Basidiomycota | Eurotiales | Onygenales | Pezizomycotina_other | Saccharomycotina | Taphrinomycotina | other_fungi | All |
|---|---|---|---|---|---|---|---|---|
| domain_non_member | 0 | 4 (7) | 0 | 13 (16) | 0 | 0 | 0 | 16 (23) |
| gpi_wall_enzyme | 0 | 6 (11) | 0 | 2 (4) | 5 (5) | 1 (1) | 0 | 8 (21) |
| laccase_etc | 7 (11) | 5 (7) | 13 (14) | 3 (4) | 1 (1) | 1 (1) | 0 | 23 (38) |
| mucin_sensor | 0 | 0 | 0 | 0 | 4 (4) | 0 | 0 | 4 (4) |
| repeat_non_adhesin | 1 (1) | 1 (2) | 0 | 1 (1) | 3 (5) | 0 | 1 (3) | 7 (12) |
| secretory_non_surface | 2 (3) | 0 | 0 | 2 (2) | 0 | 0 | 0 | 3 (5) |
| st_linker_enzyme | 3 (5) | 5 (5) | 0 | 5 (11) | 0 | 0 | 0 | 8 (21) |
| wall_hydrolase | 1 (1) | 6 (6) | 2 (3) | 6 (11) | 3 (3) | 3 (3) | 1 (1) | 15 (28) |
| wall_structural_other | 0 | 0 | 0 | 0 | 7 (9) | 0 | 0 | 7 (9) |
| **All** | 14 (21) | 27 (38) | 15 (17) | 32 (49) | 23 (27) | 5 (5) | 2 (4) | 91 (161) |

### Per clade

| Group | Rows | Sequences | Clusters | N1 rows | N2 rows | Clusters to 20 |
|---|---|---|---|---|---|---|
| Basidiomycota | 21 | 21 | 14 | 21 | 0 | 6 |
| Eurotiales | 38 | 38 | 27 | 34 | 4 | 0 |
| Onygenales | 17 | 17 | 15 | 15 | 2 | 5 |
| Pezizomycotina_other | 49 | 49 | 32 | 39 | 10 | 0 |
| Saccharomycotina | 27 | 27 | 23 | 20 | 7 | 0 |
| Taphrinomycotina | 5 | 5 | 5 | 4 | 1 | 15 |
| other_fungi | 4 | 4 | 2 | 4 | 0 | 18 |
| **all** | 161 | 161 | 91 | 137 | 24 | |

### Per species

| Group | Rows | Sequences | Clusters | N1 rows | N2 rows | Clusters to 20 |
|---|---|---|---|---|---|---|
| Agaricus bisporus (5341) | 1 | 1 | 1 | 1 | 0 | 19 |
| Aspergillus fumigatus (746128) | 24 | 24 | 16 | 20 | 4 | 4 |
| Aspergillus nidulans (162425) | 10 | 10 | 10 | 10 | 0 | 10 |
| Aspergillus niger (5061) | 1 | 1 | 1 | 1 | 0 | 19 |
| Blastomyces dermatitidis (5039) | 1 | 1 | 1 | 0 | 1 | 19 |
| Botrytis cinerea (40559) | 4 | 4 | 2 | 4 | 0 | 18 |
| Coccidioides immitis (5501) | 1 | 1 | 1 | 0 | 1 | 19 |
| Coccidioides posadasii (199306) | 3 | 3 | 3 | 3 | 0 | 17 |
| Colletotrichum lindemuthianum (290576) | 1 | 1 | 1 | 1 | 0 | 19 |
| Coprinopsis cinerea (5346) | 2 | 2 | 2 | 2 | 0 | 18 |
| Cryptococcus deneoformans (40410) | 1 | 1 | 1 | 1 | 0 | 19 |
| Cryptococcus neoformans (5207) | 5 | 5 | 5 | 5 | 0 | 15 |
| Drepanopeziza brunnea (698440) | 4 | 4 | 4 | 0 | 4 | 16 |
| Fomitopsis schrenkii (2126942) | 1 | 1 | 1 | 1 | 0 | 19 |
| Fusarium graminearum (5518) | 5 | 5 | 5 | 5 | 0 | 15 |
| Fusarium oxysporum (5507) | 1 | 1 | 1 | 1 | 0 | 19 |
| Histoplasma capsulatum (5037) | 1 | 1 | 1 | 1 | 0 | 19 |
| Histoplasma ohiense (2902605) | 1 | 1 | 1 | 1 | 0 | 19 |
| Malassezia globosa (76773) | 1 | 1 | 1 | 1 | 0 | 19 |
| Melanocarpus albomyces (204285) | 1 | 1 | 1 | 1 | 0 | 19 |
| Metarhizium anisopliae (5530) | 3 | 3 | 3 | 3 | 0 | 17 |
| Metarhizium robertsii (568076) | 1 | 1 | 1 | 1 | 0 | 19 |
| Microsporum canis (63405) | 1 | 1 | 1 | 1 | 0 | 19 |
| Moniliophthora perniciosa (153609) | 1 | 1 | 1 | 1 | 0 | 19 |
| Mortierella alpina (64518) | 2 | 2 | 1 | 2 | 0 | 19 |
| Mycetinis scorodonius (182058) | 1 | 1 | 1 | 1 | 0 | 19 |
| Mycosarcoma maydis (5270) | 3 | 3 | 3 | 3 | 0 | 17 |
| Mycothermus thermophilus (85995) | 1 | 1 | 1 | 1 | 0 | 19 |
| Neurospora crassa (5141) | 1 | 1 | 1 | 1 | 0 | 19 |
| Penicillium rubens (1108849) | 1 | 1 | 1 | 1 | 0 | 19 |
| Phanerodontia chrysosporium (2822231) | 2 | 2 | 1 | 2 | 0 | 19 |
| Pleurotus eryngii (5323) | 1 | 1 | 1 | 1 | 0 | 19 |
| Podila clonocystis (979688) | 1 | 1 | 1 | 1 | 0 | 19 |
| Podospora anserina (2587412) | 1 | 1 | 1 | 1 | 0 | 19 |
| Pyricularia oryzae (318829) | 2 | 2 | 2 | 2 | 0 | 18 |
| Rhizomucor miehei (4839) | 1 | 1 | 1 | 1 | 0 | 19 |
| Saccharomyces cerevisiae (4932) | 27 | 27 | 23 | 20 | 7 | 0 |
| Schizosaccharomyces pombe (4896) | 5 | 5 | 5 | 4 | 1 | 15 |
| Talaromyces marneffei (37727) | 1 | 1 | 1 | 1 | 0 | 19 |
| Talaromyces purpureogenus (1266744) | 1 | 1 | 1 | 1 | 0 | 19 |
| Thermothelomyces thermophilus (78579) | 1 | 1 | 1 | 1 | 0 | 19 |
| Trametes cinnabarina (5643) | 1 | 1 | 1 | 1 | 0 | 19 |
| Trametes hirsuta (5327) | 1 | 1 | 1 | 1 | 0 | 19 |
| Trichoderma harzianum (5544) | 5 | 5 | 4 | 5 | 0 | 16 |
| Trichoderma reesei (51453) | 6 | 6 | 5 | 3 | 3 | 15 |
| Trichophyton rubrum (5551) | 9 | 9 | 8 | 9 | 0 | 12 |
| Verticillium dahliae (27337) | 12 | 12 | 12 | 9 | 3 | 8 |
| **all** | 161 | 161 | 91 | 137 | 24 | |

### Rows and clusters per module served

| serves | Rows | Clusters | Leakage none | Leakage other |
|---|---|---|---|---|
| step1 | 9 | 7 | 9 | 0 |
| repeat | 82 | 47 | 82 | 0 |
| family_domain | 29 | 20 | 6 | partial=20, in_reference=3 |
| adhesion_level | 156 | 88 | 156 | 0 |

### Per stratum, N1 rows only

| Group | Rows | Sequences | Clusters | N1 rows | N2 rows | Clusters to 20 |
|---|---|---|---|---|---|---|
| domain_non_member | 12 | 12 | 9 | 12 | 0 | 11 |
| gpi_wall_enzyme | 21 | 21 | 8 | 21 | 0 | 12 |
| laccase_etc | 37 | 37 | 23 | 37 | 0 | 0 |
| mucin_sensor | 3 | 3 | 3 | 3 | 0 | 17 |
| repeat_non_adhesin | 11 | 11 | 6 | 11 | 0 | 14 |
| secretory_non_surface | 5 | 5 | 3 | 5 | 0 | 17 |
| st_linker_enzyme | 18 | 18 | 7 | 18 | 0 | 13 |
| wall_hydrolase | 24 | 24 | 13 | 24 | 0 | 7 |
| wall_structural_other | 6 | 6 | 4 | 6 | 0 | 16 |
| **all** | 137 | 137 | 76 | 137 | 0 | |

### Per stratum, without rows that share a cluster with a positive

| Group | Rows | Sequences | Clusters | N1 rows | N2 rows | Clusters to 20 |
|---|---|---|---|---|---|---|
| domain_non_member | 19 | 19 | 15 | 12 | 7 | 5 |
| gpi_wall_enzyme | 21 | 21 | 8 | 21 | 0 | 12 |
| laccase_etc | 38 | 38 | 23 | 37 | 1 | 0 |
| mucin_sensor | 4 | 4 | 4 | 3 | 1 | 16 |
| repeat_non_adhesin | 12 | 12 | 7 | 11 | 1 | 13 |
| secretory_non_surface | 5 | 5 | 3 | 5 | 0 | 17 |
| st_linker_enzyme | 21 | 21 | 8 | 18 | 3 | 12 |
| wall_hydrolase | 26 | 26 | 14 | 23 | 3 | 6 |
| wall_structural_other | 9 | 9 | 7 | 6 | 3 | 13 |
| **all** | 155 | 155 | 89 | 136 | 19 | |

Note on the `serves` table: `adhesion_level` means a call that asks "adhesin or not" after step 1.
`step1` rows are signal-peptide or membrane proteins that are **not** at the surface (ER, vacuole,
plasma membrane sensors). A secreted or wall enzyme is a **positive** for step 1 (surface or
secreted), so it does not serve as a step 1 negative. `repeat` is set when UniProt lists 2 or more
Repeat features, or when a 30-residue window is at least 50% Ser+Thr, or by stratum. `family_domain`
is set when the protein has a Pfam domain of a `pfam_adhesion` family in
`data/sorting_hat/family_table.tsv`.

## What exists outside this set (task step 1 inventory)

| Source | Content | Counts |
|---|---|---|
| `data/curated/adhesins/hard_negative_seeds.tsv` | Seed list | 29 data rows (28 *S. cerevisiae* gene queries, 1 *A. fumigatus* abr1) plus 2 comment lines and a header; the task file says 31. All 29 resolve to rows of `adhesins.tsv`. |
| `data/curated/adhesins/adhesins.tsv`, `cls=hard_negative` | Draft hard negatives | 42 rows: 32 N1, 10 N2; 28 *S. cerevisiae*, 10 *A. fumigatus*, 1 each *H. capsulatum*, *C. immitis*, *T. marneffei*, and 1 *A. fumigatus* (taxon 746128, aspf2). 40 are in `controls.tsv`; MSB2 and aspf2 are in `rejected.tsv`. |
| `data/controls/repeat-adhesins-other-species/hard_negatives_other_species.tsv` (task 08) | Look-alike negatives | 10 rows. 1 kept (Blastomyces ENG2, N2); 9 rejected (hydrophobins x4, Hsp60, Cpl1, Pyricularia MSB2, Blumeria LIP1, Zymoseptoria GT2). |
| `_workdir/step1_compare/phasec/eval_table.tsv.gz` (Phase C, step 1) | GO-derived step 1 negatives | N-sec rows: *S. cerevisiae* 1545, *S. pombe* 1203 (+2 N-int,N-sec), *A. nidulans* 1064, *A. fumigatus* 1060 (+2 shared with *A. nidulans*), *C. albicans* 961, *U. maydis* 827, *C. neoformans* 25. PM-TM rows (labelled P-ext, stratum PM-TM): *C. albicans* 22, *S. pombe* 4, *S. cerevisiae* 3, *A. nidulans* 2, *A. fumigatus* 1, *U. maydis* 1. These are step 1 negatives by GO rule, not curated look-alikes. They are not copied here. |
