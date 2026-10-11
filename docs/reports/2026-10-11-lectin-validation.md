# Lectin-like Pfam models: recovery of known fungal lectins and rejection of look-alikes

*2026-10-11. Author: Claude Code (claude-sonnet-5-5). Branch `lectin-validation`. Revision 2, after independent review 1 (`2026-10-11-lectin-validation-review-1.md`, verdict pass with changes; six must-fix items M1 to M6 and the minor items are fixed here). Task: measure only. The owner decision of the lectin spec (D1, option (c)) stands: there is no lectin call in v1 and no status file. This work creates no call, no status file and no threshold. Code and tables: `analysis/lectin_validation/` (log: `LOG.md`; unchecked claims: `unverified_notes.md`). Tests: `tests/lectin_validation/` (34 pass). Every number below comes from a file in `analysis/lectin_validation/results/` or `controls/`.*

## 1. Summary

1. **Truth set.** 59 primary positives (59 sequences, 32 clusters at 30% identity, 38 species) and 166 reasoned negatives (166 sequences, 84 clusters, 54 species). 10 PA14 positives are kept apart (the PA14 models are not among the 15). A weak `assumed` set (497) is kept only as a construction check (appendix A).
2. **Recall.** The union of the 15 M1/M2 models hits 42 of 59 positives (0.71; Wilson 0.59 to 0.81) and 21 of 32 clusters (0.66; 0.48 to 0.80). PF07938 alone hits 6 of 59 (4 of 32 clusters; 0.05 to 0.28).
3. **Part of this is by construction.** 10 positives (5 clusters) entered the candidate list only through a UniProt Pfam cross-reference to a model under test (`selection_route = model_gated_pfam_xref`): LysM effectors and one PF07938 lectin. All 10 are hit. This breaks the spirit of COMMON-RULES rule 3. Without them, the 49 rows that were found by UniLectin3D, a keyword or a literature search give 32 of 49 (0.51 to 0.77) and 17 of 28 clusters (0.42 to 0.76). **LysM 12/12, `chitin_masking_h1b` 5/5 clusters and the T1 result 17/18 are partly by construction.**
4. **Recall is mostly coverage of known domains.** 36 of 59 positives have 30% identity or more to a Pfam seed sequence (`in_reference`). Cluster-level split: the 19 clusters that hold at least one such row are hit 19/19 (40 of 40 sequences). The 13 clusters with no such row are hit **2 of 13** (0.04 to 0.42): CNL (B2ZRS9) and G1X5Y2. (Before G1X5Y2 was added the figure was 1 of 12, 0.01 to 0.35. The first version of this report said 3 of 14; it counted two clusters twice.) Four positives are 98% or more identical to a seed. Most `in_reference` rows are so because of one short domain (for example a 44-residue LysM domain).
5. **False hits.** The union of the 15 models hits 0 of 136 N1 negatives (0 of 74 clusters; 0.00 to 0.05) and 0 of 18 ER house-keeping proteins (0 of 8 clusters; 0.00 to 0.32; **by construction**, section 5.2). It hits 8 of 12 enzymes with a lectin-like module; those 8 are one cluster and one of them is identical to a PF00652 seed. PF07938 alone hits none of the 166 reasoned negatives. **No negative was chosen as a close relative of PF07938. This report makes no PF07938 specificity claim** (section 6).
6. **Misses.** 17 positives (11 clusters). 11 sequences (8 clusters) belong to classes with no model among the 15. 6 sequences (3 clusters) are ricin-B-type beta-trefoil lectins that score below the GA of an existing model, RicinB_lectin_2 (24.6 to 24.9 for three, 22.0, 21.2 and 15.3 against 25.6). PF14200 family recall with them: 6 of 12 sequences, 5 of 8 clusters (0.31 to 0.86). I set no threshold.
7. **PF07938 score floors.** No negative hits PF07938, so a floor changes no false hit. The six PF07938-family positives score 525.1, 489.1, 300.7, 158.6, 132.6 and 60.7. Floors of 50: 6/6. 100: 5/6 (G1X5Y2 is lost). 200: 3/6.
8. **Report script.** Same hit sets and scores as my own hmmsearch (0 mismatches). Only 18 of the 150 class-by-model cells are non-empty, so this checks parsing and grouping, not the models.
9. **Status.** Not written, by decision. If a call existed, `smoke` is the best any entry could reach (section 7).

## 2. What this is and what it is not

The class is defined by Pfam models. The truth is partly the structures from which the same families were described, and part of the literature rows were found with the same models. So the recall below is mostly the coverage of known lectins by their own families. It is not a measure of how well the models find unknown lectins. T4 (review judgement) and any label that comes from a model hit are not used as truth. The 15 models are the M1 and M2 models of spec 3.1. They are reported, not called. The five M3 models (PF03388, PF00262, PF07915, PF02157, PF11721) are scanned for the negatives and are not in the union.

## 3. Data

### 3.1 Positives (primary: 59 rows)

**One written rule for tiers (LOG decision 12).** T1: I read in PubMed an abstract that reports glycan or chitin binding for the named protein. T2: a UniLectin3D curated structure, or a Swiss-Prot reviewed entry whose FUNCTION text has an ECO:0000269 sugar-binding statement, whether or not the cited abstract shows the assay. Near-identical copies of one gene in another strain are not added again. Fragments under 50 aa are not added. Q06100 and G1X5Y2 are T2 under this one rule. The abstract of PMID 8999822 does not say which isolectin was expressed, and the abstract of PMID 35142828 reports a deletion phenotype and no binding assay. If the owner prefers a stricter T2, both rows move out together.

**Selection route (rule 3).** Column `selection_route`: `unilectin3d` (43 rows), `literature_or_keyword` (6: Slp1, Avr4, paracoccin, CNL, Cgl-I, CMRBL; the first two and paracoccin came from the spec, CNL and Cgl-I through the Swiss-Prot Lectin keyword, CMRBL from a PubMed search), `model_gated_pfam_xref` (10: F9XMT4, G2X4U8, H1UVH3, H1VSQ9, J4URT3, J5J7Y8, F2ST24, F2SIK4, P9WEM2, G1X5Y2). The last group entered the candidate list only because UniProt cross-references the PF01476 or PF07938 model to them (`make_acc_list.py`). A LysM effector that PF01476 misses could not enter this way.

| Stratum | Rows | Clusters | Species | Tier | Sources |
|---|---|---|---|---|---|
| `structure_backed` | 39 | 23 | 27 | T2 | UniLectin3D fungal structures (193 structures fetched 2026-10-11; 55 of 56 accessions re-found), RCSB ligand names |
| `host_binding_h1a` | 4 | 3 | 3 | T1 | FleA Q4WW81 (PMIDs 23695231, 24340081, 27058347), SapL1 A0A084FYP2 (34373510), *T. rubrum* LysM1 F2ST24 and LysM2 F2SIK4 (30947244) |
| `chitin_masking_h1b` | 11 | 5 | 6 | T1 | Ecp6 (20724636), Mg1LysM and Mg3LysM (21467214), Mgx1LysM (33797163), Slp1 (22267486), ELP1 and ELP2 (27174033), Vd2LysM (27911046), Blys2 and Blys5 (28873459), Avr4 (17153926) |
| `literature_binding` | 4 | 4 | 4 | 2 T1, 2 T2 | CNL B2ZRS9 (19100814); CMRBL G3JGU5 (36632960, glycan array, new); Cgl-I Q06100 (T2); G1X5Y2 (T2) |
| `known_gap` | 1 | 1 | 1 | T1 | Paracoccin C1G842 (16476564, 24743161, 25474158) |
| **All primary** | **59** | **32** | **38** | 18 T1, 41 T2 | |
| `pa14_owned_elsewhere` (apart) | 10 | 6 | 4 | T2 | UniLectin3D PA14 class |

Notes.
- **Strata share clusters.** The stratum cluster counts add to 36, not 32. Four clusters span two strata (F2ST24: h1a and h1b; P18891 and Q2UNX8: h1a and structure; Q9P4R8: literature and structure; `results/strata_shared_clusters.tsv`). Each cluster is counted once in the 32.
- The 43 UniLectin3D rows are the 43 sequences of spec 4.4. One *Ralstonia* row was rejected (bacterial). The two Emp rows are negatives.
- The two *T. rubrum* LysM proteins bind N-linked oligosaccharides of human skin glycoproteins **and also** chitin and N-linked oligosaccharides of fungal glycoproteins (abstract of PMID 30947244). The binding is not host-specific. The H1a label is my reading of the abstract. The owner decides (section 9).
- `host_binding_h1a` is small. FleA binding is disputed (Sakai 2019, inherited, not checked in the PDF).
- New in this revision: CMRBL (the gene CCM_03832 named in the abstract maps to UniProt G3JGU5), Mgx1LysM, G1X5Y2.
- Candidates not added, with reasons (`controls/rejected.tsv`, 78 lines): BbLec1 (the paper has no open full text; the gene ID could not be checked; J4UPT0 is a candidate only), CmLec4, CmRlec, RsLysM, Cg2LysM, TAL6, MfAvr4 and the Cercospora Avr4 orthologs (no verifiable accession), D4AL26 (no binding evidence), *Histoplasma* and *Talaromyces* activities without a gene, CFP28 (no binding shown), duplicates, fragments, and the *Ralstonia* lectin. UPI000150AC21 (a UniParc ID) is rejected: no sequence in the committed file.

### 3.2 Negatives

| Class | Reason (source) | Rows | Clusters | Species | Selection |
|---|---|---|---|---|---|
| `er_housekeeping` | ER or Golgi house-keeping lectin-fold protein: calnexin, OS-9 (MRH), Emp46p, Emp47p (UniProt FUNCTION text; spec 3.1 row M3; PMID 12134087). 6 rows with ECO:0000269, 12 by similarity | 18 | 8 | 13 | by UniProt xref to an M3 model |
| `enzyme_lectin_like_domain` | Enzyme function: 8 alpha-galactosidases (PF00652 module; Swiss-Prot, by similarity or no code) and 4 PA14 beta-glucosidases (shared set, N2) | 12 | 2 | 7 | Swiss-Prot keyword or xref; shared set |
| `n1_shared` | Characterised non-adhesive function with a source (shared set N1; 137 minus CNE1, counted once as `er_housekeeping`) | 136 | 74 | 44 | shared set |
| **All reasoned** | | **166** | **84** | **54** | |

- `er_housekeeping` proteins bind glycans. The label means "not a host-interaction lectin" (spec decision). It does not mean "no glycan binding".
- `enzyme_lectin_like_domain` is ambiguous. UniProt gives the alpha-galactosidases the keyword Lectin. The 8 are one cluster. P28351 is 100% identical to the PF00652 seed AGALA_ASPNG; the others are 47 to 96%. So the 8 hits are domain hits on a seed-family protein. They are never pooled into a headline rate without the split.
- One cluster mixes classes: paracoccin C1G842 clusters with two N1 endochitinases (Q12713, Q4U4T0). A GH18 model for paracoccin would also call those N1 negatives.
- The `assumed` set (497 reviewed fungal Swiss-Prot proteins) is not a false-hit set. See appendix A.

## 4. Method

- Truth table: `build_truth.py` (inputs: UniLectin3D API, RCSB GraphQL, UniProt REST release 2026_03, `curated_positives.tsv`, the shared hard-negative set). Species and taxon ids: species rank of the UniProt organism in `nodes.dmp`; names from `names.dmp` (names and ranks are checked in the tests; one row, A0A1U7Q1Z0, is at genus rank).
- Clusters: MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5`, all rows together. A cluster is a hit when any member in the group hits. **Groups that split rows by a flag (leakage) are built at the cluster level** (`lectin_metrics.cluster_seed_split`): a cluster with any flagged row goes to the flagged group. Row-level splits are also written, with a label.
- Scan: `hmmsearch --cut_ga` (HMMER 3.4), 20 models from the Pfam 38.2 file. Score floors use the tblout sequence score.
- Wilson 95% intervals on clusters (and on sequences; the sequence interval ignores cluster structure and is too narrow).
- Leakage: seed alignments of the 20 models from InterPro (versions equal to Pfam 38.2). `in_reference` = 30% identity or more, E 1e-3 or less, target coverage 0.5 or more (mmseqs2 `-s 7.5`). It reproduces spec F28. Review 1 re-ran it with 0 label differences. The identity is over the aligned part. Short domains can reach 30% by chance-like similarity to a bacterial, phage or algal seed domain.
- Miss analysis: `hmmsearch -E 10` and `hmmscan --cut_ga` against all of Pfam-A on the 17 missed positives.
- Re-run: `run_all.sh` (login node; includes the Swiss-Prot background scan; the network steps are separate scripts).

## 5. Results

### 5.1 Recall at the Pfam GA cutoff (primary positives; Wilson in brackets)

| Model set | Sequences | Clusters |
|---|---|---|
| Union of 15 (M1+M2), all primary | 42/59 (0.59-0.81) | 21/32 (0.48-0.80) |
| ...independent route only (49 rows) | 32/49 (0.51-0.77) | 17/28 (0.42-0.76) |
| ...model-gated route only (10 rows) | 10/10 (0.72-1.00) | 5/5 (0.57-1.00) |
| ...T1 rows only | 17/18 (0.74-0.99) | 9/10 (0.60-0.98) |
| Five revision-1 models | 18/59 (0.20-0.43) | 8/32 (0.13-0.42) |
| PF07938 only | 6/59 (0.05-0.20) | 4/32 (0.05-0.28) |
| LysM PF01476 | 13/59 | 5/32 |
| RicinB_lectin_2 PF14200 | 7/59 | 6/32 |
| Ricin_B_lectin PF00652 | 1/59 | 1/32 |
| Five M3 models (union) | 0/59 | 0/32 |

On the 43 UniLectin3D rows the revision-1 models and PF07938 give the numbers of spec F10 (16 of 43 and 5 of 43 sequences; 7 and 3 clusters).

By stratum (union of 15; clusters): `structure_backed` 13/23; `host_binding_h1a` 3/3; `chitin_masking_h1b` 5/5; `literature_binding` 4/4; `known_gap` 0/1. **The strata share clusters** (section 3.1), so these do not add up to 21. The 5/5 of `chitin_masking_h1b` is partly by construction (7 of its 11 rows came through the PF01476 cross-reference; Ecp6 and Mg1LysM are UniLectin3D rows; Slp1 and Avr4 came from the spec). By phylum: Ascomycota 25/26 sequences (13/14 clusters), Basidiomycota 17/33 (8/18 clusters). The phylum split mostly shows how I chose the rows.

**Family recall** (expected model: spec F11 class mapping, and the family named in the paper for literature rows; `results/family_recall.tsv`). Every family with a model is hit in full except RicinB_lectin_2: PF07938 6/6, LysM 12/12, galectin 5/5, ricin-B (PF00652) 1/1, CVNH 4/4, FB_lectin 5/5. **This is partly by construction**: the F11 mapping used Pfam hits, the seeds hold the same or close sequences, and the model-gated rows were found with the models. PF14200: **6 of 12 sequences, 5 of 8 clusters** (0.31 to 0.86) when the six beta-trefoil rows (the UniLectin3D class "Laetiporus trefoil lectin-like") are counted as this family; 6 of 6 sequences and 5 of 5 clusters without them. I map them to PF14200 because they are ricin-B-type beta-trefoil lectins and five of them hit PF14200 below GA; the choice is mine and is stated here. They fall in 3 clusters (R4GRU5 with 3 other *Boletus* paralogs, B0D650, Q7Z8V1), and the five near-GA rows in 2 clusters.

**Seed-independent test (cluster level).** Within families that have a model: PF07938 has 1 cluster with no seed-similar row (G1X5Y2, hit); PF14200 has 4 such clusters, 1 hit (CNL; 1 of 7 sequences); LysM, galectin, CVNH, FB_lectin and PF00652 have none. All primary positives: 2 of 13 clusters (sections 1 and 5.4).

Hits close to the cutoff (`results/ga_margins.tsv`): the five smallest margins over GA are 5.5 (galectin Cgl-III), 6.6 (Avr4 on CBM_14), 8.1 (galectin), 11.6 (Mg1LysM) and 13.2 bits.

### 5.2 False hits (negatives)

| Negative class | Union of 15 (sequences) | Union of 15 (clusters) | PF07938 only | Five M3 models |
|---|---|---|---|---|
| `n1_shared` | 0/136 (0.00-0.03) | 0/74 (0.00-0.05) | 0/136 | 0/136 |
| `er_housekeeping` | 0/18 (0.00-0.18) | 0/8 (0.00-0.32) | 0/18 | 18/18 |
| `enzyme_lectin_like_domain` | 8/12 (0.39-0.86) | 1/2 (0.09-0.91) | 0/12 | 0/12 |
| all reasoned | 8/166 (0.02-0.09) | 1/84 (0.00-0.06) | 0/166 (0.00-0.02) | 18/166 |

- The 8 enzyme hits are the 8 alpha-galactosidases, one cluster. They are hit by Ricin_B_lectin PF00652 (scores 47.8 to 64.6) and two also by RicinB_lectin_2 (28.8, above its GA of 25.6). P28351 is identical to a PF00652 seed. The PA14 beta-glucosidases hit no model.
- **The M3 result is by construction.** The 18 ER rows were selected by a UniProt cross-reference to an M3 model, and the M3 models hit them at GA. The 0 of 18 for the 15 models is also expected: M3 and M1/M2 are different families.
- **Meaning.** The 15 models do not hit the house-keeping lectins or the N1 set. That is expected from how these sets were chosen. None was picked as a close relative of PF07938 or the other 14 families. The N1 set is not a look-alike set (spec 4.3). The enzyme result shows a domain-level hit on a protein that is an enzyme.
- Swiss-Prot background (`results/swissprot_background.summary.tsv`; description, not a measure): in all 38,451 reviewed fungal proteins the 15 models give 97 protein-model hits (0.25% of proteins). PF07938 hits 7. All 7 are named fucose-specific lectins (AOL and AFL under two accessions each, AAL, D4AL26, G1X5Y2). LysM hits 53; most are named LysM effectors or LysM-domain proteins. An unlabelled hit is not a false hit.
- D4AL26 (*A. benhamiae*, Swiss-Prot "Fucose-specific lectin") has protein-level evidence: it was found in the culture secretome by mass spectrometry (UniProt, PMID 21919205). Its lectin function is by similarity to AAL (ECO:0000250). There is no binding evidence. It is 92.3% identical to the *T. benhamiae* member (F3F37487_006503-T1) of the six Onygenales PF07938 proteins in spec F27 (`results/d4al26_vs_six.m8`). I did not check whether they are the same gene. So one of the six is secreted in culture, with no binding evidence.

### 5.3 PF07938 score floors (spec 3.1 option (b))

| Floor | PF07938-family positives (sequences / clusters) | All primary, union with floor | Reasoned negatives |
|---|---|---|---|
| GA (25.3) | 6/6; 4/4 | 42/59; 21/32 | 8/166 (0 from PF07938) |
| 50 | 6/6; 4/4 | 42/59; 21/32 | same |
| 100 | 5/6; 3/4 | 41/59; 20/32 | same |
| 200 | 3/6; 2/4 | 39/59; 19/32 | same |

The six truth scores are 525.1, 489.1, 300.7, 158.6, 132.6 and 60.7 (G1X5Y2). No negative hits PF07938 at any score, so these data cannot say which floor is better. Five of the six are `in_reference` (30% identity or more to the 3-sequence seed), so they cannot set a cutoff (spec F28). G1X5Y2 has no seed hit, but one positive does not set a cutoff either.

### 5.4 Leakage (cluster level) and identity to the seed

| Group | Sequences hit | Clusters hit |
|---|---|---|
| Clusters with at least one `in_reference` row | 40/40 (0.91-1.00) | 19/19 (0.83-1.00) |
| Clusters with no `in_reference` row | 2/19 (0.03-0.31) | **2/13 (0.04-0.42)** |

Row-level counts (36 `in_reference`, 23 `none`; `results/leakage_split.tsv`) are labelled as such in the table, because clusters can hold both kinds of row (F2ST24, Q9P4R8). The 13 seed-independent clusters are mostly families without a model (17 of their 19 sequences are misses). So seed-independence and "no model" are confounded. The two hits are CNL and G1X5Y2.

Best identity to a seed of any of the 20 models (rows; `results/identity_bins.tsv`; clusters are not exclusive across bins):

| Identity to nearest seed | Rows hit |
|---|---|
| 98% or more | 4/4 |
| 90 to 97.9% | 1/1 |
| 50 to 89.9% | 15/15 |
| 30 to 49.9% | 16/16 |
| no hit at 30% | 6/23 |

So only 5 positives are near-identical to a seed. The 30 to 49.9% bin is mostly one LysM domain (about 44 residues) against bacterial, phage or algal seed domains. For such rows `in_reference` means "similar to a seed domain", not "in the training set". I kept the COMMON-RULES 30% rule. Whether to add a near-identity tier is an owner question (section 9).

### 5.5 Which known lectins the models miss, and why

17 positives (11 clusters) have no hit among the 15 models at GA (`results/misses.tsv`, column `reason`).

| Positives | Why missed |
|---|---|
| **Below GA of an existing model (6 sequences, 3 clusters):** beta-trefoil lectins R4GRU4, R4GRU5, R4GRU6, R4GRU9 (one cluster), B0D650, Q7Z8V1 | RicinB_lectin_2 scores 24.9, 24.7, 24.6, 21.2, 22.0 and 15.3 against GA 25.6. Q7Z8V1 also hits the toxin model ETX_MTX2 at GA. No threshold is set |
| **No model among the 15 (11 sequences, 8 clusters):** | |
| Paracoccin C1G842 (1) | No lectin model by decision D9. Its only Pfam hit is Glyco_hydro_18 (36.0) |
| PVL-like lectins Q309D1, A0A1U7Q1Z0, H6CS64 (3) | Hit FG-GAP_3 (79.9 to 98.9), not a lectin model |
| Tectonin 2 B0CZL6 (1) | A Pfam model (Tectonin PF19193) hits at 201.0. It is not among the 15 |
| FIP-Fve P80412 (1) | A Pfam model (Fve PF09259) hits at 128.6. It is not among the 15 |
| CBM67-like E7E2M2 (1) | Sub-GA hit to Laminin_G_3 (22.1 against 27.6), probably unrelated |
| PhoSL prism lectins A0A384E107, A0A3B6UEU4 (2) | No Pfam hit; the sequences are 40 and 43 aa |
| CML-like B3VS76, monocot-like I7H471 (2) | No Pfam hit |

I found no missed positive that is a truncated gene model; I did not compare lengths with full-length family members. Host-binding proteins in this set are hit (FleA, SapL1, two *T. rubrum* LysM), except paracoccin, which binds host TLR N-glycans. Histoplasma and Coccidioides items have no sequence or no binding data and could not be tested.

### 5.6 Check of the report script

`lectin_report.py` was run on one FASTA per class (10 classes as 10 proteomes). For each class and each of the 15 models, the protein set and best score equal my own `hmmsearch` (0 of 150 cells differ). **Only 18 of the 150 cells are non-empty.** Both sides use HMMER 3.4 `--cut_ga` and the same Pfam file. So this checks that the script parses and groups hits correctly. It is not an independent test of the models, and it does not cover the signal peptide path or the (a2) section.

## 6. Specificity of PF07938: what can be said

- All 166 reasoned negatives and the 497 construction-check rows give 0 PF07938 hits. **None of them is a relative of PF07938**, so this is not a specificity measure against look-alikes.
- I looked for look-alikes with a reason. In Swiss-Prot (reviewed) every protein with a PF07938 cross-reference is a fucose-specific lectin (7 of 7). Among all UniProt fungal proteins with a PF07938 cross-reference (949; `data/raw/uniprot_fungi_pf07938_xref.tsv.gz`), 52 have a non-lectin, non-"uncharacterized" name: serine/threonine and other kinases, manganese lipoxygenases, NIDO-domain and metallopeptidase-domain proteins, nucleoside phosphorylase-domain proteins, ADP-ribosylhydrolase. All 52 are TrEMBL (protein existence "inferred from homology" for 50, transcript level for 2). Their names are automatic annotation. I found no paper that gives a reason that any of them is not a fucose-binding lectin. **So no reasoned look-alike negative could be built.**
- J4KR73 (*B. bassiana*, "Peptidase M64, IgA", protein existence "inferred from homology", PF07938 cross-reference) is a candidate. It has no paper. It is listed in `unverified_notes.md` and `rejected.tsv`, not used as a negative. Q4WMH4 (spec 4.3) is another candidate.
- The multi-domain proteins with a PF07938 cross-reference (69 of the 949 entries have more than one Pfam family) show that the model hits domains inside larger proteins. This is a reason to build such negatives before any specificity claim, if the owner wants one (section 9).

## 7. Size against the floor, and the best status

| Target (COMMON-RULES rule 7) | This set |
|---|---|
| 20 positive clusters | 32 pooled (14 Ascomycota, 18 Basidiomycota; 10 with a T1 row; 13 with no seed-similar row). Not met in any one species: the most is 3 (*Coprinopsis cinerea*); *A. fumigatus* has 1 |
| 20 negative clusters (N1 rows only count, spec F3) | 74 N1 clusters pooled. Not met in any one species: the most is 16 (*S. cerevisiae*), then 15 (*A. fumigatus*) |
| Half-width 0.10 or less | Not met. Union recall by cluster: half-width 0.16 (0.48 to 0.80) |

- I read `status_from_measure` and `call_status._entries` (review 1 did too). A measure with fewer than 20 positive and 20 negative clusters, or a half-width above 0.10, is at most `smoke`. A leakage value other than `none` caps an entry at `smoke`. A species without a positive is `unvalidated`.
- **Plainly: if a lectin call existed, `smoke` would be the best status for any entry, pooled or by species, from the size rule alone. The leakage cap is a second reason** (36 of 59 positives are `in_reference`). `unvalidated` would apply to the 150 species without a positive (38 of 188 species have one). *A. fumigatus* would be `smoke` on FleA (1 cluster, `in_reference`). **I write no status.** Under option (c) there is no call.
- The pooled sets are report tables, not species entries (hydrophobin decision D1).

## 8. What the numbers say and do not say

They say:
- The 15 models, at the Pfam GA, find the fungal lectins in families that have a model and a seed (42 of 59), including the ones that came by literature search (32 of 49 independent-route rows).
- They hit none of the 136 N1 and 18 ER house-keeping proteins.
- They miss families with no model (11 sequences) and ricin-B-type lectins below the RicinB_lectin_2 GA (6 sequences). They cannot tell a lectin from an enzyme that carries a ricin-B-like domain.
- Host-glycan-binding proteins in this set are hit, except paracoccin.

They do not say:
- That the models recover unknown lectins. Seed-independent clusters: 2 of 13.
- That PF07938 is specific (section 6).
- That a hit shows host binding, surface location or adhesion (spec fixed text, 6.2 item 3).
- Anything about a species. The truth is pooled over 38 species.
- That LysM, `chitin_masking_h1b` or T1 recall is an independent result. It is partly by construction.

Circularity. (1) The class is defined by Pfam models and positives of mapped families are similar to the seeds. (2) The F11 class-to-model mapping used Pfam hits. (3) Ten positives were found through the models (`model_gated_pfam_xref`). (4) The `none` leakage group is mostly families without a model. (5) The ER negatives were found through M3 models. T4 and review judgement are used nowhere as truth.

## 9. Limits

- Pooled over species. 59 positives, 32 clusters. The intervals are wide. Wilson on clusters assumes clusters at 30% identity are independent.
- Strata are small and share clusters.
- UniLectin3D ligands were not checked per structure (39 rows are T2). A7UNK4 PDB ids were not retrieved (API cap).
- The negatives are not look-alikes of the models. The only look-alike class (`enzyme_lectin_like_domain`) is one cluster of 8 plus one of 4, and is ambiguous.
- The 12 `er_housekeeping` rows without an experimental code have function text by similarity.
- The report script check covers parsing of hit sets only.
- Data sources are live web resources (UniProt 2026_03, InterPro, UniLectin3D, RCSB, PubMed) on 2026-10-11. The data files are saved in `data/`, except the two raw Swiss-Prot lists (fetch command in `controls/README.md`).
- The control set is in `analysis/lectin_validation/controls/`, not in `data/controls/` as COMMON-RULES asks. This is an exception recorded here: the set is a measurement table, not a calibration set.

## 10. Decisions for the owner

1. **Candidate pool gated by the models (review 1, M2).** Accept the 10 model-gated rows with the disclosure above, or rebuild the literature positives from searches that do not use Pfam cross-references? Default here: keep, report both.
2. **Drop `assumed` from the false-hit tables?** Done here (appendix A only). Reverse if you disagree.
3. **T2 rule.** One rule applied to both: G1X5Y2 and Q06100 are in. If you want a stricter T2 (the abstract must show the assay), both leave and PF07938 has 5 positives again.
4. **PF07938 look-alikes.** Build reasoned look-alike negatives (kinase, lipoxygenase, peptidase M64 and similar multi-domain proteins) before any specificity claim? It needs a paper per protein or an experimental function. J4KR73 and Q4WMH4 are starting points.
5. **Leakage for short domains.** Keep the 30% rule for LysM-sized domains, or add a near-identity tier (for example 90%)? Only 5 positives are 90% or more identical to a seed.
6. **New candidates.** CMRBL is in. BbLec1 is out until its gene is checked in the full text (J4UPT0 is the candidate). Other candidates (CmLec4, CmRlec, RsLysM, Cg2LysM, TAL6, MfAvr4, Cercospora Avr4) need accessions.
7. **Add the two *T. rubrum* LysM proteins to the lectin spec as H1a evidence?** They bind host and fungal N-glycans and chitin. The binding is not host-specific. It could feed a later LysM sub-call (D2). It does not change D1.
8. **`enzyme_lectin_like_domain`.** Keep as an apart class? The question is whether the report labels proteins or domains.
9. **Models outside the 15.** Tectonin PF19193 hits Tectonin 2. Fve PF09259 hits FIP-Fve, whose lectin status is weaker. Decide separately.
10. **The beta-trefoil lectins.** Six rows (3 clusters; five rows in 2 clusters score 0.7 to 4.4 bits below the GA) are misses of an existing model at its GA. I set no threshold. Whether the PF14200 cutoff or the family needs attention is an owner question.
11. **Use of the tables later.** Report-grade only, not a calibration set, until the owner has read `controls/counts.md` and `unverified_notes.md`. The weakest rows: 39 `structure_backed` T2 rows, the 10 model-gated rows, Vd2LysM, and the `er_housekeeping` rows by similarity.

## Appendix A. `assumed` set: a construction check

497 reviewed fungal Swiss-Prot proteins (430 clusters, 148 species) with no lectin keyword, no UniProt Pfam cross-reference to any of the 20 models and no lectin word in name or function. They have no reason that they are not lectins. The filter removes every protein that UniProt already matches to these models, so a hit is possible only where the UniProt Pfam release and Pfam 38.2 differ. The result, 0 of 497 (0 of 430 clusters; 0.00 to 0.01) for the union and for PF07938, **is by construction and has no power**: the base rate in all of reviewed fungal Swiss-Prot is about 0.25% of proteins (97 of 38,451), so even an unfiltered sample of 497 would give about one hit. It shows only that the scan runs and that the filter works (`results/construction_check.tsv`). The Swiss-Prot background in section 5.2 is the better description.

## Files

- Results: `analysis/lectin_validation/results/` (`counts.tsv`, `recall_by_group.tsv`, `false_hits.tsv`, `construction_check.tsv`, `floors.tsv`, `family_recall.tsv`, `leakage_split.tsv`, `identity_bins.tsv`, `strata_shared_clusters.tsv`, `misses.tsv`, `ga_margins.tsv`, `species_clusters.tsv`, `per_row_hits.tsv`, `mixed_clusters.tsv`, `report_script_check.tsv`, `swissprot_background.*`, `d4al26_vs_six.m8`).
- Truth set: `analysis/lectin_validation/controls/` (`controls.tsv`, `sequences.faa`, `clusters.tsv`, `rejected.tsv`, `counts.md`, `README.md`, `unverified_notes.md`).
- Code: `analysis/lectin_validation/*.py`, `run_all.sh`. Tests: `tests/lectin_validation/`.
- Sources: PubMed (abstracts of the PMIDs named); UniProtKB 2026_03; UniLectin3D (CC BY 4.0); RCSB PDB; InterPro (Pfam seeds); Pfam 38.2.
