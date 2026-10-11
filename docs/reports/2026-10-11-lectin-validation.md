# Lectin-like Pfam models: recovery of known fungal lectins and rejection of look-alikes

*2026-10-11. Author: Claude Code (claude-sonnet-5-5). Branch `lectin-validation` (from `lectin-report-script`, merged with `lectin-class-spec-r5`). Task: measure only. The owner decision of the lectin spec (D1, option (c)) stands: there is no lectin call in v1 and no status file. This work creates no call, no status file and no threshold. Code and tables: `analysis/lectin_validation/` (log: `LOG.md`; unchecked claims: `unverified_notes.md`). Tests: `tests/lectin_validation/` (23 pass). Every number below comes from a file in `analysis/lectin_validation/results/` or `controls/`.*

## 1. Summary

1. **Truth set.** 56 primary positives (56 sequences, 30 clusters at 30% identity, 37 species) and 166 reasoned negatives (166 sequences, 84 clusters, 54 species). A separate weak set of 497 `assumed` negatives (430 clusters). 10 PA14 positives are kept apart (the PA14 models are not among the 15).
2. **Recall.** The union of the 15 M1/M2 models hits 39 of 56 positives (0.70; Wilson 0.57 to 0.80) and 19 of 30 clusters (0.63; 0.46 to 0.78). PF07938 alone hits 5 of 56 (3 of 30 clusters; 0.10; 0.03 to 0.26).
3. **The recall is mostly coverage of known structures.** 34 of the 56 positives (18 of 30 clusters) have 30% identity or more to a Pfam seed sequence (`in_reference`). All 34 are hit. Of the 22 positives without such a seed hit, 5 are hit (3 of 14 clusters; 0.08 to 0.48).
4. **False hits.** The union of the 15 models hits 0 of 136 N1 negatives (0 of 74 clusters; 0.00 to 0.05), 0 of 18 ER house-keeping proteins (0 of 8 clusters; 0.00 to 0.32), 0 of 497 assumed negatives (0 of 430 clusters; 0.00 to 0.01). It hits 8 of 12 enzymes with a lectin-like module (alpha-galactosidases with a ricin-B-like model hit; 1 of 2 clusters). PF07938 alone hits none of the 166 reasoned negatives (0 of 84 clusters; 0.00 to 0.04). **No negative was chosen as a close relative of PF07938**, so this is not a test of specificity against relatives.
5. **Misses.** 17 positives (11 clusters) are missed by all 15 models: paracoccin (1), and 16 structure-backed mushroom lectins. Three PVL-like lectins hit only FG-GAP_3 (not a lectin model). Tectonin 2 and FIP-Fve hit Pfam models that are not among the 15. Six beta-trefoil lectins are below the GA; five of them have sub-GA hits to RicinB_lectin_2 that lie 0.7 to 4.4 bits below it. Five have no lectin-family hit (CBM67-like, two PhoSL, CML-like, monocot-like). I found no missed positive that is a truncated gene model.
6. **PF07938 score floors** (50, 100, 200) change no false hit (there is none). At 50 and 100 they lose no positive. At 200 they lose two of the five PF07938-family positives (AAL, 158.6; *Arthrobotrys* AOL, 132.6).
7. **The report script** (`lectin_report.py`) gives the same hit sets and scores as my own `hmmsearch` for 10 class FASTA files and 15 models (150 cells, 0 mismatches).
8. **Status.** Not written, by decision. If a call existed, the best status that any species entry could reach is `smoke`. No species has more than 3 primary positive clusters, and no species has 20 reasoned negative clusters (maximum 19). `estimated` is not reachable (section 7).

## 2. What this is and what it is not

The class is defined by Pfam models. The truth is partly the structures from which the same families were described. So the recall below is mostly the coverage of known lectins by their own families. It is not a measure of how well the models find unknown lectins. T4 (review judgement) and any label that comes from a model hit are not used as truth. The 15 models are the M1 and M2 models of spec 3.1. They are reported, not called. The five M3 models (PF03388, PF00262, PF07915, PF02157, PF11721) are scanned for the negatives and are not in the union.

## 3. Data

### 3.1 Positives (primary: 56 rows)

Source rules: every row has a source and a quote (`controls/controls.tsv`). Tier T1: I read in PubMed an abstract that reports glycan binding for the protein. Tier T2: UniLectin3D curated structure, or Swiss-Prot entry with an ECO:0000269 sugar-binding statement. All rows are fungal proteins with a UniProt accession (one UniParc ID, UPI000150AC21, has a sequence from the earlier spec work).

| Stratum | Rows | Sequences | Clusters | Species | Tier | Sources |
|---|---|---|---|---|---|---|
| `structure_backed` | 39 | 39 | 23 | 27 | T2 | UniLectin3D fungal structures (193 structures fetched 2026-10-11; 55 of 56 accessions re-found), RCSB ligand names |
| `host_binding_h1a` | 4 | 4 | 3 | 3 | T1 | FleA Q4WW81 (PMIDs 23695231, 24340081, 27058347), SapL1 A0A084FYP2 (34373510), *T. rubrum* LysM1 F2ST24 and LysM2 F2SIK4 (30947244) |
| `chitin_masking_h1b` | 10 | 10 | 5 | 6 | T1 | Ecp6 (20724636), Mg1LysM and Mg3LysM (21467214), Slp1 (22267486), ELP1 and ELP2 (27174033), Vd2LysM (27911046), Blys2 and Blys5 (28873459), Avr4 (17153926) |
| `literature_binding` | 2 | 2 | 2 | 2 | T1, T2 | CNL B2ZRS9 (19100814, glycan array); Cgl-I Q06100 (Swiss-Prot ECO:0000269, 8999822) |
| `known_gap` | 1 | 1 | 1 | 1 | T1 | Paracoccin C1G842 (16476564, 24743161, 25474158) |
| **All primary** | **56** | **56** | **30** | **37** | 16 T1, 40 T2 | |
| `pa14_owned_elsewhere` (apart) | 10 | 10 | 6 | 4 | T2 | UniLectin3D PA14 class |

Notes.
- The 43 UniLectin3D rows (39 + FleA + SapL1 + Ecp6 + Mg1LysM) are the 43 sequences of spec 4.4 (PA14, Emp L-type and *Ralstonia* removed). One *Ralstonia* row (bacterial) was rejected. The two Emp rows are negatives (below).
- Beyond the 43 UniLectin3D rows I added 13 rows, all with a PMID or a Swiss-Prot record that I verified: F9XMT4, G2X4U8, G4N906, H1UVH3, H1VSQ9, J4URT3, J5J7Y8, Q00363, F2ST24, F2SIK4, B2ZRS9, Q06100 and C1G842. See `controls/controls.tsv`, column `source`.
- **New in this work, not in the spec:** the two *T. rubrum* LysM proteins bind N-linked oligosaccharides of human skin glycoproteins (abstract of PMID 30947244). That meets the H1a definition of spec 3.1. They are LysM (an M2 model) in an Onygenales species. Apart from paracoccin (a known gap), they are the only host-glycan-binding proteins in this set outside the FleA family.
- `host_binding_h1a` is small. FleA binding is disputed (Sakai 2019, inherited, not checked in the PDF).
- Rejected candidates (76 lines, `controls/rejected.tsv`): BbLec1 (no verifiable accession), *Histoplasma* and *Talaromyces* activities without a gene, CFP28 (no binding shown), G1X5Y2 and D4AL26 (tier T3: a deletion phenotype only, and annotation by similarity), duplicates, fragments, and the *Ralstonia* lectin.

### 3.2 Negatives

| Class | Reason (source) | Rows | Sequences | Clusters | Species |
|---|---|---|---|---|---|
| `er_housekeeping` | ER or Golgi house-keeping lectin-fold protein: calnexin, OS-9 (MRH), Emp46p, Emp47p (UniProt FUNCTION text; spec 3.1 row M3; PMID 12134087). 6 rows with ECO:0000269, 12 by similarity | 18 | 18 | 8 | 13 |
| `enzyme_lectin_like_domain` | Enzyme function: 8 alpha-galactosidases (PF00652 module; Swiss-Prot, by similarity or no code) and 4 PA14 beta-glucosidases (shared set, N2) | 12 | 12 | 2 | 7 |
| `n1_shared` | Characterised non-adhesive function with a source (shared set N1; 137 minus CNE1, counted once as `er_housekeeping`) | 136 | 136 | 74 | 44 |
| **All reasoned** | | **166** | **166** | **84** | **54** |
| `assumed` (weak, apart) | **No reason.** Reviewed fungal Swiss-Prot, no lectin keyword, no UniProt Pfam cross-reference to the 20 models, no lectin word; random sample, at most 15 per species | 497 | 497 | 430 | 148 |

- `er_housekeeping` proteins bind glycans. The label means "not a host-interaction lectin" (spec decision). It does not mean "no glycan binding".
- `enzyme_lectin_like_domain` is ambiguous. UniProt gives the alpha-galactosidases the keyword Lectin. I report it apart and never inside a pooled false-hit rate without the split.
- `assumed` was chosen with UniProt Pfam cross-references. So a zero in this class is partly by construction. It tests only that nothing common in Swiss-Prot hits by accident.
- One cluster mixes classes: paracoccin C1G842 clusters with two N1 negatives (Q12713, Q4U4T0), which are endochitinases (Q12713 *Trichoderma harzianum* endochitinase 33; Q4U4T0 *Metarhizium anisopliae* endochitinase 2). So a GH18 model for paracoccin would also call those N1 negatives.

## 4. Method

- Truth table: `build_truth.py` (inputs: UniLectin3D API, RCSB GraphQL, UniProt REST release 2026_03, `curated_positives.tsv`, the shared hard-negative set). Species and taxon ids: species rank of the UniProt organism in `/srv/projects/db/taxonomy/nodes.dmp`; names from `names.dmp` (checked in the tests).
- Clusters: MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5`, all 729 sequences together. A cluster is a hit when any member in the group hits.
- Scan: `hmmsearch --cut_ga` (HMMER 3.4) with the 20 models extracted from `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`. A protein is a hit when it has a sequence-level hit at GA (tblout). Score floors use the tblout sequence score.
- Wilson 95% intervals on clusters (and on sequences, shown for comparison; the sequence interval ignores cluster structure and is too narrow).
- Leakage: seed alignments of the 20 models from InterPro (versions equal to the local Pfam 38.2 file). `in_reference` = 30% identity or more, E 1e-3 or less, target coverage 0.5 or more (mmseqs2 `-s 7.5`). The method reproduces spec F28 (FleA 81.7%, SapL1 43.7%, *Arthrobotrys* AOL 36.6%, AAL 32.8%, *A. oryzae* AOL 100%).
- Miss analysis: `hmmsearch -E 10` of the 20 models and `hmmscan --cut_ga` against all of Pfam-A on the 17 missed positives.
- Re-run: `run_all.sh` (login node; the network steps are separate scripts and their outputs are in `data/`).

## 5. Results

### 5.1 Recall at the Pfam GA cutoff (primary positives; cluster-level Wilson in brackets)

| Model set | Sequences | Clusters |
|---|---|---|
| Union of 15 (M1+M2) | 39/56 (0.57-0.80) | 19/30 (0.46-0.78) |
| Five revision-1 models | 17/56 (0.20-0.43) | 7/30 (0.12-0.41) |
| PF07938 only | 5/56 (0.04-0.19) | 3/30 (0.03-0.26) |
| LysM PF01476 | 12/56 (0.13-0.34) | 5/30 (0.07-0.34) |
| RicinB_lectin_2 PF14200 | 6/56 (0.05-0.21) | 5/30 (0.07-0.34) |
| CVNH PF08881 | 4/56 | 2/30 |
| FB_lectin PF07367 | 5/56 | 1/30 |
| Gal-bind_lectin PF00337 | 5/56 | 1/30 |
| CBM_14 PF01607, LDL PF21691, CCL2-like PF21595 | 1/56 each | 1/30 each |
| PF18647, PF22803, PF01419, PF00652, PF13385, PF00187 | 0/56 | 0/30 |
| Five M3 models (union) | 0/56 | 0/30 |

The five revision-1 models and PF07938 give the numbers of spec F10 (16 of 43 and 5 of 43 sequences on the 43 UniLectin3D rows; 7 and 3 clusters). My 17 includes Cgl-I (a literature row).

By stratum (union of 15; clusters): `structure_backed` 13/23 (0.37-0.74); `host_binding_h1a` 3/3 (0.44-1.00); `chitin_masking_h1b` 5/5 (0.57-1.00); `literature_binding` 2/2; `known_gap` 0/1. T1 rows only: 15/16 sequences, 8/9 clusters (0.57-0.98). By phylum: Ascomycota 22/23 sequences (11/12 clusters), Basidiomycota 17/33 (8/18 clusters; 0.25-0.66). The phylum split mostly shows how I chose the rows: the literature rows are Ascomycota effectors of the LysM family, all in the seed.

Family recall (the expected model of a row: spec F11 class mapping, and the family named in the paper for literature rows): every family that has a model among the 15 is hit in full (for example PF07938 5/5, LysM 11/11, galectin 5/5, ricin-B-like 6/6; `results/family_recall.tsv`). **This is partly by construction.** The F11 mapping was made with Pfam hits, and the seeds contain the same or close sequences.

Rows with a model among the 15 and no seed hit (`not_seed_similar`): LysM 1/1, galectin 3/3, RicinB_lectin_2 1/1, i.e. 5 sequences in 3 clusters, all hit. That is the only seed-independent test of recall here. It is too small to interpret.

Hits close to the cutoff (`results/ga_margins.tsv`): the five smallest margins over GA are 5.5 (galectin Cgl-III), 6.6 (Avr4 on CBM_14), 8.1 (galectin), 11.6 (Mg1LysM on LysM), 13.2 bits.

### 5.2 False hits (negatives)

| Negative class | Union of 15 (sequences) | Union of 15 (clusters) | PF07938 only | Five M3 models |
|---|---|---|---|---|
| `n1_shared` | 0/136 (0.00-0.03) | 0/74 (0.00-0.05) | 0/136 | 0/136 |
| `er_housekeeping` | 0/18 (0.00-0.18) | 0/8 (0.00-0.32) | 0/18 | 18/18 |
| `enzyme_lectin_like_domain` | 8/12 (0.39-0.86) | 1/2 (0.09-0.91) | 0/12 | 0/12 |
| all reasoned | 8/166 (0.02-0.09) | 1/84 (0.00-0.06) | 0/166 (0-0.02) | 18/166 |
| `assumed` (weak) | 0/497 (0.00-0.01) | 0/430 (0.00-0.01) | 0/497 | 0/497 |

- The 8 enzyme hits are the alpha-galactosidases hit by Ricin_B_lectin PF00652 (8; scores 47.8 to 64.6) and RicinB_lectin_2 PF14200 (2, at 28.8, below 30). The PA14 beta-glucosidases hit no model (PA14 is not in the 15).
- The five M3 models hit all 18 ER house-keeping proteins (6 by PF00262, 10 by PF07915, 2 by PF03388). They do what the spec wants of them.
- **Meaning.** The 15 models do not hit the house-keeping lectins, the N1 set or the assumed set. That is expected from how these sets were chosen. None of them was picked as a close relative of PF07938 or of the other 14 families. The N1 set is not a look-alike set (spec 4.3 says the same). The result for the enzymes shows a protein-level false hit class: a ricin-B-like model hits enzymes that have a ricin-B-like module.
- Swiss-Prot background (`results/swissprot_background.summary.tsv`; description, not a measure): in all 38,451 reviewed fungal proteins the 15 models give 97 protein-model hits. PF07938 hits 7 proteins. All 7 are named fucose-specific lectins (AOL and AFL under two accessions each, AAL, D4AL26 and G1X5Y2). LysM hits 53. Most are named LysM effectors or LysM-domain proteins. A few are chitinases, a deacetylase or a killer toxin. The unlabelled hits are listed in `swissprot_background.unlabelled.tsv`. An unlabelled hit is not a false hit.
- D4AL26 (*A. benhamiae*, Swiss-Prot "Fucose-specific lectin", evidence by similarity to AAL) is 92.3% identical to the *T. benhamiae* member (F3F37487_006503-T1) of the six Onygenales PF07938 proteins in spec F27 (`results/d4al26_vs_six.m8`). So one of the six has a Swiss-Prot name by similarity. It has no experimental evidence. I did not check whether they are the same gene.

### 5.3 PF07938 score floors (spec 3.1 option (b))

| Floor | PF07938-family positives (sequences / clusters) | All primary positives, union with floor | Reasoned negatives | Assumed |
|---|---|---|---|---|
| GA (25.3) | 5/5; 3/3 | 39/56; 19/30 | 8/166 (0 from PF07938) | 0/497 |
| 50 | 5/5; 3/3 | 39/56; 19/30 | same | same |
| 100 | 5/5; 3/3 | 39/56; 19/30 | same | same |
| 200 | 3/5; 2/3 | 37/56; 18/30 | same | same |

The five PF07938 truth scores are 525.1, 489.1, 300.7, 158.6 and 132.6. A floor of 200 loses AAL and *Arthrobotrys* AOL. No negative hits PF07938 at any score, so these data cannot say which floor is better. All five positives are `in_reference` (30% identity or more to the 3-sequence seed), so they cannot set a cutoff (spec F28).

### 5.4 Leakage

| Group | Hit / n (sequences) | Hit / n (clusters) |
|---|---|---|
| Primary positives, `in_reference` | 34/34 (0.90-1.00) | 18/18 (0.82-1.00) |
| Primary positives, `none` | 5/22 (0.10-0.43) | 3/14 (0.08-0.48) |

The `none` group is mostly the families without a model (17 of its 22 rows are misses), so the lower recall there is confounded with "no model exists". The five hits are CNL, the three galectins and *T. rubrum* LysM2. In per-model terms: PF07938 has no positive without a seed hit (0 of 5), so its recall on unseen sequences is not measured. The `leakage` values in the table are `in_reference` or `none`. None is `unknown`. `none` means no hit at 30% identity by mmseqs2. It does not exclude remote homology.

### 5.5 Which known lectins the models miss, and why

17 positives (11 clusters) have no hit among the 15 models at GA (`results/misses.tsv`).

| Positives | Why missed |
|---|---|
| Paracoccin C1G842 (1) | No lectin model by decision D9. Its only Pfam hit is Glyco_hydro_18 (score 36.0, partial) |
| PVL-like lectins Q309D1, A0A1U7Q1Z0, H6CS64 (3) | No model in the 15. They hit FG-GAP_3 (79.9 to 98.9), which is not a lectin model |
| Tectonin 2 B0CZL6 (1) | A Pfam model (Tectonin PF19193) hits at 201.0. It is not among the 15 |
| FIP-Fve P80412 (1) | A Pfam model (Fve PF09259) hits at 128.6. It is not among the 15 |
| *Laetiporus*-type beta-trefoil lectins R4GRU4, R4GRU5, R4GRU6, R4GRU9, B0D650, Q7Z8V1 (6) | Below GA. Five have sub-GA hits to RicinB_lectin_2 (24.9, 24.7, 24.6, 21.2, 22.0 against GA 25.6). Q7Z8V1 hits at 15.3 and hits the toxin model ETX_MTX2 at GA |
| CBM67-like E7E2M2 (1) | Sub-GA hit to Laminin_G_3 (22.1 against GA 27.6), probably unrelated; Bac_rhamnosid_N at GA |
| PhoSL prism lectins A0A384E107, A0A3B6UEU4 (2) | No Pfam hit. The sequences are 40 and 43 aa |
| CML-like B3VS76, monocot-like I7H471 (2) | No Pfam hit |

I did not find a missed positive that is a truncated gene model. I did not compare lengths with full-length family members (`unverified_notes.md`, item 18). All 17 are in classes that spec F11 lists without a lectin model among the 15 (F11 notes the Tectonin and Fve hits).

For the host-interaction question: the host-binding proteins in this set are not missed, except paracoccin. FleA, SapL1 and the two *T. rubrum* LysM proteins are hit. Histoplasma and Coccidioides items have no sequence or no binding data, so they could not be tested (`rejected.tsv`).

### 5.6 Check of the report script

`analysis/lectin_class_report/lectin_report.py` was run on one FASTA per class (10 classes as 10 proteomes). For each class and each of the 15 models, the set of proteins with a hit and the best sequence score equal my own `hmmsearch` (0 of 150 cells differ; `results/report_script_check.tsv`). This checks the script on the M1 and M2 models only. It does not check the signal peptide path or the (a2) section.

## 6. What the numbers say and do not say

They say:
- The 15 models, at the Pfam GA, find the fungal lectins that belong to families with a model and a seed (39 of 56; all of the 34 `in_reference` rows).
- They hit none of the 136 N1, 18 ER house-keeping and 497 assumed negatives.
- They do not find the families that have no model (17 of 56 positives), and the models do not tell a lectin from an enzyme that carries a ricin-B-like module.
- Host-glycan-binding proteins in this set are hit (FleA, SapL1, two *T. rubrum* LysM), except paracoccin, which binds host TLR N-glycans and is a known gap.

They do not say:
- That the models recover unknown lectins. The truth is mostly in the seeds. Where a family has a model, the positives that are not close to a seed are 5 sequences in 3 clusters.
- That PF07938 is specific. No negative is a relative of PF07938.
- That a hit shows host binding, surface location or adhesion (spec fixed text, 6.2 item 3). Only the 16 T1 rows have such a paper, and the FleA role is disputed.
- Anything about a species. The truth is pooled over 37 species.
- That the pooled recall would hold in a proteome. The positives are curated, well-studied proteins. The Swiss-Prot background shows that the models also hit many proteins that nobody labelled; those are not scored.

Circularity. (1) The class is defined by Pfam models and the positives of mapped families are in the seeds. (2) The class-to-model mapping of spec F11 used Pfam hits. (3) The `none` leakage group is mostly families without a model, so seed-independence and "no model" are confounded. (4) The `assumed` set was selected with a Pfam filter. T4 and review judgement are used nowhere as truth.

## 7. Size against the floor, and the best status

| Target (COMMON-RULES rule 7) | This set |
|---|---|
| 20 positive clusters | 30 pooled (12 Ascomycota, 18 Basidiomycota; 9 with a T1 row; 14 without a seed hit). Not met in any one species: the most is 3 (*Coprinopsis cinerea*); *A. fumigatus* has 1 |
| 20 negative clusters | 84 reasoned pooled (74 are N1). Not met in any one species: the most is 19 (*A. fumigatus*, *S. cerevisiae*) |
| Half-width 0.10 or less | Not met. Union recall by cluster: half-width 0.16 (0.46 to 0.78) |

- The spec (4.4, statement 2) says that leakage other than `none` caps a status at `smoke`. 34 of the 56 positives (18 of 30 clusters) are `in_reference`. I did not re-read the code.
- **Plainly: if a lectin call existed, `smoke` would be the best status for any species entry, and `unvalidated` for any species without a positive** (150 of the 187 species in this table have none). *A. fumigatus* would be `smoke` on FleA (1 cluster, `in_reference`, 81.7% to a seed), which carries no information on sensitivity (spec 4.4). **I write no status.** Under option (c) there is no call to give a status to.
- The pooled sets are report tables, not species entries (hydrophobin decision D1, as the spec says).

## 8. Limits

- Pooled over species. No species-level measurement exists.
- 56 positives, 30 clusters. The intervals are wide. Wilson on clusters assumes independence of clusters at 30% identity.
- Strata are small: `host_binding_h1a` 3 clusters, `chitin_masking_h1b` 5, `known_gap` 1, `literature_binding` 2.
- UniLectin3D ligands were not checked per structure (`unverified_notes.md`). The 39 structure rows are T2.
- A7UNK4 (LDL) PDB ids and ligands were not retrieved (API cap).
- The negatives are not look-alikes of the models. The only look-alike class (`enzyme_lectin_like_domain`) has 2 clusters and is ambiguous.
- The 12 `er_housekeeping` rows without an experimental code have function text by similarity.
- The assumed set has no reason and a built-in bias.
- Some positives are Ascomycota effectors in the seeds of LysM. The Basidiomycota rows are mostly mushroom lectin structures.
- The report script check covers hit sets only.
- Data sources are live web resources (UniProt 2026_03, InterPro, UniLectin3D, RCSB, PubMed) on 2026-10-11. The data files are saved in `data/`.

## 9. Decisions for the owner

1. **Add the two *T. rubrum* LysM proteins (F2ST24, F2SIK4) to the lectin spec as H1a evidence in an Onygenales species?** They bind N-linked oligosaccharides of human skin glycoproteins (PMID 30947244). They hit LysM, an M2 model that is common in fungi (spec 3.1). This is a possible input to a later LysM sub-call (spec D2). It does not change decision D1.
2. **Accept `enzyme_lectin_like_domain` as an apart class?** It holds protein-level negatives that carry a real ricin-B-like module. A reader may call that a hit, not a false hit.
3. **Models outside the 15.** Tectonin PF19193 and Fve PF09259 hit two known lectins that the 15 miss. Should the report list them as M2? (This changes the model list of the spec. I did not change it.)
4. **The near-GA beta-trefoil family.** Five *Laetiporus*-type lectins score 0.7 to 4.4 bits below the RicinB_lectin_2 GA. I set no threshold. Whether this is a family-model gap is for the owner.
5. **Use of these tables later.** The set is a report-grade table. It is not a calibration set until the owner has read `controls/counts.md` and `unverified_notes.md` (COMMON-RULES review gate). The weakest rows are the 39 `structure_backed` T2 rows (ligand not checked per structure), Vd2LysM, and the `er_housekeeping` rows by similarity.

## 10. Files

- Report data: `analysis/lectin_validation/results/` (`counts.tsv`, `recall_by_group.tsv`, `false_hits.tsv`, `floors.tsv`, `family_recall.tsv`, `leakage_split.tsv`, `misses.tsv`, `ga_margins.tsv`, `species_clusters.tsv`, `per_row_hits.tsv`, `mixed_clusters.tsv`, `report_script_check.tsv`, `swissprot_background.*`, `d4al26_vs_six.m8`).
- Truth set: `analysis/lectin_validation/controls/` (`controls.tsv`, `sequences.faa`, `clusters.tsv`, `rejected.tsv`, `counts.md`, `README.md`, `unverified_notes.md` copy).
- Code: `analysis/lectin_validation/*.py`, `run_all.sh`. Tests: `tests/lectin_validation/`.
- Sources: PubMed (abstracts of the PMIDs named in section 3.1); UniProtKB release 2026_03; UniLectin3D (CC BY 4.0); RCSB PDB; InterPro (Pfam seeds); Pfam 38.2.
