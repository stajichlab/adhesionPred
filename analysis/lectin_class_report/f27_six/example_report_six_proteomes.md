# Lectin-like family report (no call)

Date: 2026-10-10. Spec: docs/superpowers/specs/2026-10-10-lectin-class-design.md (option c-i; sections 3.1 and 6.2).

This report makes no call. It writes no status file. It does not claim a function for any protein. Every model here is reported only (spec decision D1, option c).

## Method

- Search: `hmmsearch --cut_ga` (HMMER 3.4 (Aug 2023); http://hmmer.org/) with the models of the table below against Pfam 38.2 (`/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`).
- The cutoff is the Pfam GA cutoff of each model. A hit is a protein with a sequence score at or above the sequence GA (the `--tblout` table). Scores in the tables are sequence scores.
- Column `domain_ge_ga` (hits table) and `n_hit_domain_ge_ga` (counts table) show whether the best domain also reaches the domain GA. The tool's own Pfam table comes from the domain table, so the tool can count fewer proteins for a model than this report (for example RicinB_lectin_2 in C. neoformans H99).
- Proteomes: 6. No signal peptide or other condition is applied to the hits.
- Signal peptide data: supplied (`analysis/lectin_class_report/f27_six/signalp6_prediction_results.txt`).

## PF07938 hits per proteome

Label of every PF07938 hit: **PF07938 (AAL-like / FleA-family) domain hit**.

Fixed text for every PF07938 hit: "This is a hit to the PF07938 (AAL-like / FleA-family) domain model. It does not show that the protein is on the cell surface or secreted. It does not show that the protein binds a host glycan. Host-glycan binding is shown in papers for FleA and SapL1 only; the role of FleA is disputed."

| proteome | proteins | PF07938 hit proteins | best score | proteins with score >= 100 |
|---|---|---|---|---|
| Arthroderma_uncinatum_CBS_119779 | 7939 | 2 | 291.1 | 2 |
| Microsporum_canis_CBS_113480 | 8161 | 2 | 296.6 | 1 |
| Myotisia_sp._PD_48_PD_48 | 7564 | 2 | 330.6 | 1 |
| Nannizzia_gypsea_CBS_118893 | 7970 | 1 | 290.4 | 1 |
| Onygena_corvina_CBS_281.48 | 7101 | 1 | 291.0 | 1 |
| Trichophyton_benhamiae_CBS_112371 | 7420 | 1 | 288.3 | 1 |

Total PF07938 hit proteins: 9 in 6 of 6 proteomes.

## Lectin-like families (not called)

Counts per model over all input proteomes. H is the host-interaction evidence level of spec 3.1. U is the whole-proteome ubiquity of spec 3.1 (653 genera of the discovery pass). H and U are copied from the spec. This script does not measure them. The column "proteomes with hit" is for the proteomes of this run.

| set | model | name | H | U (spec 3.1) | hit proteins | proteomes with hit | best score |
|---|---|---|---|---|---|---|---|
| M1 | PF07938 | Fungal_lectin | H1a | 392 genera (60.0%); 104 (15.9%) at score 100 or more | 9 | 6 | 330.6 |
| M2 | PF01476 | LysM | H1b | not scanned; at least 511 (78%), pool S value | 97 | 6 | 99.5 |
| M2 | PF08881 | CVNH | H2 | not scanned; at least 233 (36%), pool S value | 3 | 3 | 94.2 |
| M2 | PF07367 | FB_lectin | H2 | 68 genera (10.4%) | 0 | 0 | NA |
| M2 | PF21691 | LDL | H0 | 36 genera (5.5%) | 3 | 3 | 34.9 |
| M2 | PF00337 | Gal-bind_lectin | H0 | 31 genera (4.7%) | 0 | 0 | NA |
| M2 | PF21595 | CCL2-like | H0 | 2 genera (0.3%) | 0 | 0 | NA |
| M2 | PF18647 | Fungal_lectin_2 | H0 | not scanned; at least 170 (26%), pool S value | 6 | 4 | 175.6 |
| M2 | PF22803 | GBD_Y3 | H0 | not scanned; at least 51 (8%), pool S value | 0 | 0 | NA |
| M2 | PF01419 | Jacalin | H0 | not scanned; at least 256 (39%), pool S value | 2 | 1 | 31.0 |
| M2 | PF00652 | Ricin_B_lectin | H0 | not scanned; at least 218 (33%), pool S value | 0 | 0 | NA |
| M2 | PF14200 | RicinB_lectin_2 | H0 | not scanned; at least 122 (19%), pool S value | 0 | 0 | NA |
| M2 | PF13385 | Laminin_G_3 | H0 | not scanned; at least 54 (8%), pool S value | 6 | 6 | 31.9 |
| M2 | PF01607 | CBM_14 | H1b | 31 genera (4.7%) | 2 | 2 | 37.3 |
| M2 | PF00187 | Chitin_bind_1 | H1b-like | 572 genera (87.6%) | 29 | 6 | 61.1 |

M1 here means the models with H1a evidence (PF07938 only, spec 3.1). M1 is reported, not called. M2 models are reported, not called.

## Per proteome and per model

Hit proteins (best score) per proteome. A dash means no hit. The full table is `lectin_report_counts.tsv`.

| proteome | Fungal_lectin | LysM | CVNH | FB_lectin | LDL | Gal-bind_lectin | CCL2-like | Fungal_lectin_2 | GBD_Y3 | Jacalin | Ricin_B_lectin | RicinB_lectin_2 | Laminin_G_3 | CBM_14 | Chitin_bind_1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Arthroderma_uncinatum_CBS_119779 | 2 (291.1) | 18 (62.0) | - | - | 1 (27.6) | - | - | - | - | - | - | - | 1 (31.9) | - | 4 (40.8) |
| Microsporum_canis_CBS_113480 | 2 (296.6) | 25 (95.6) | 1 (91.3) | - | - | - | - | 2 (171.3) | - | - | - | - | 1 (29.0) | 1 (32.7) | 7 (53.6) |
| Myotisia_sp._PD_48_PD_48 | 2 (330.6) | 14 (89.0) | - | - | 1 (34.9) | - | - | 1 (158.1) | - | 2 (31.0) | - | - | 1 (30.8) | - | 5 (45.6) |
| Nannizzia_gypsea_CBS_118893 | 1 (290.4) | 19 (99.5) | 1 (94.2) | - | - | - | - | 2 (173.2) | - | - | - | - | 1 (28.3) | 1 (37.3) | 8 (61.1) |
| Onygena_corvina_CBS_281.48 | 1 (291.0) | 7 (52.1) | - | - | - | - | - | - | - | - | - | - | 1 (31.3) | - | 2 (48.0) |
| Trichophyton_benhamiae_CBS_112371 | 1 (288.3) | 14 (84.3) | 1 (93.7) | - | 1 (27.7) | - | - | 1 (175.6) | - | - | - | - | 1 (27.9) | - | 3 (42.1) |

## (a2) group: PF07938 score >= 100 and a signal peptide

This group is named for the record only. It is the group of spec option (a2). No call is made. No function is claimed. Orthology to FleA is not checked. The spec states that the floor of 100 was fitted to the same data that found the group.

PF07938 hits at score 100 or more: 6 with a signal peptide, 0 without, 1 with no signal peptide data.

| proteome | protein | score | label |
|---|---|---|---|
| Arthroderma_uncinatum_CBS_119779 | FA49FE15_006535-T1 | 291.1 | PF07938 (AAL-like / FleA-family) domain hit |
| Microsporum_canis_CBS_113480 | FCE2ACBB_004168-T1 | 296.6 | PF07938 (AAL-like / FleA-family) domain hit |
| Myotisia_sp._PD_48_PD_48 | FC452C90_000306-T1 | 330.6 | PF07938 (AAL-like / FleA-family) domain hit |
| Nannizzia_gypsea_CBS_118893 | F94BF6A1_004981-T1 | 290.4 | PF07938 (AAL-like / FleA-family) domain hit |
| Onygena_corvina_CBS_281.48 | F1215351_005656-T1 | 291.0 | PF07938 (AAL-like / FleA-family) domain hit |
| Trichophyton_benhamiae_CBS_112371 | F3F37487_006503-T1 | 288.3 | PF07938 (AAL-like / FleA-family) domain hit |

Fixed text for every row: "This is a hit to the PF07938 (AAL-like / FleA-family) domain model. It does not show that the protein is on the cell surface or secreted. It does not show that the protein binds a host glycan. Host-glycan binding is shown in papers for FleA and SapL1 only; the role of FleA is disputed."

Spec F27 names six full-length PF07938 proteins (score 288 to 331) in six Onygenales genera. No paper on the function of these proteins was found: a PubMed search for a dermatophyte AAL-like or fucose-binding lectin gave 0 hits (spec F27, review 3). "Not found" is not "does not exist". Kerr 2016 (PMID 27058347) reports putative FleA orthologs in dermatophytes, found by BLAST of AAL against NCBI sequences (S1 Fig legend). Orthology of these proteins to FleA is not checked by this report.

### Named proteins

These proteins were named by the person who ran the report (--name-list). Each row shows its PF07938 result in the input proteomes. No call is made. No function is claimed.

| protein | note | proteome | PF07938 score | signal peptide | label |
|---|---|---|---|---|---|
| FC452C90_000306-T1 | Myotisia sp. PD_48 (F27) | Myotisia_sp._PD_48_PD_48 | 330.6 | yes | PF07938 (AAL-like / FleA-family) domain hit |
| FCE2ACBB_004168-T1 | Microsporum canis (F27) | Microsporum_canis_CBS_113480 | 296.6 | yes | PF07938 (AAL-like / FleA-family) domain hit |
| FA49FE15_006535-T1 | Arthroderma uncinatum (F27) | Arthroderma_uncinatum_CBS_119779 | 291.1 | yes | PF07938 (AAL-like / FleA-family) domain hit |
| F1215351_005656-T1 | Onygena corvina (F27) | Onygena_corvina_CBS_281.48 | 291.0 | yes | PF07938 (AAL-like / FleA-family) domain hit |
| F94BF6A1_004981-T1 | Nannizzia gypsea (F27) | Nannizzia_gypsea_CBS_118893 | 290.4 | yes | PF07938 (AAL-like / FleA-family) domain hit |
| F3F37487_006503-T1 | Trichophyton benhamiae (F27) | Trichophyton_benhamiae_CBS_112371 | 288.3 | yes | PF07938 (AAL-like / FleA-family) domain hit |

## Notes

- M3 house-keeping models are not scanned. They are listed so the reader sees what is left out:
  - PF03388 Lectin_leg-like: ER-to-Golgi cargo receptors (Emp46p, Emp47p).
  - PF00262 Calreticulin: calnexin and calreticulin.
  - PF07915 and PF02157 PRKCSH and Man-6-P_recep: P-type lectin (MRH) models.
  - PF11721 Malectin: ER glycoprotein quality control.
- Models owned by the narrow tier (PA14 PF07691, GLEYA PF10528, Candida_ALS_N PF11766) are not part of this report.
- Known gap (spec decision D9, no custom model): paracoccin of Paracoccidioides (UniProt C1G842) is a published GlcNAc-binding, laminin-binding protein (Coltri 2006, PMID 16476564). Its only Pfam hit is a partial GH18 chitinase domain (Glyco_hydro_18 PF00704). Its predicted lectin sites lie inside that domain. No lectin model of this report hits it, so this report does not show it.
- Coccidioides, Histoplasma and Blastomyces: the literature searches of spec F17 found no fungal lectin literature. The hits in this report are domain hits only.
- The list of models comes from a name search for lectin-like words (spec F4). A carbohydrate-binding protein whose model has no such word can be missed.
- The report items 4 and 5 of spec 6.2 (UniLectin3D sequence table, T4 outcomes) are not made by this script.
- Signal peptide data and the name list are matched on the protein ID alone, not on the proteome name. Equal IDs in two proteomes would share one value.
