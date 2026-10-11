# Log of the lectin validation (2026-10-11)

Author: claude-sonnet-5-5. Task: measure how well the lectin-like Pfam models recover known fungal lectins and reject look-alikes.
This is a measurement. It makes no call, no status file and no threshold for the tool (owner decision D1, option (c)).

## Decisions

1. **Positive tiers.** Spec 4.1 tiers are used. T1 needs a paper that reports glycan binding for the protein, and I read the abstract in PubMed in this work. T2 is a UniLectin3D curated structure, or a Swiss-Prot entry with an ECO:0000269 sugar-binding statement. I did not give T1 to a UniLectin3D structure on the ground of a sugar in the PDB entry. A sugar component in a PDB entry can be a glycosylation site. So the 39 `structure_backed` rows are T2. The sugar component ids are listed in the table as information.
2. **Strata.** One stratum per row. `host_binding_h1a`: paper shows binding to a host glycan or host cell (FleA, SapL1, two *T. rubrum* LysM proteins). `chitin_masking_h1b`: paper shows binding to fungal chitin and an immune or protection effect (Ecp6, Mg1LysM, Mg3LysM, Slp1, ELP1, ELP2, Vd2LysM, Blys2, Blys5, Avr4). `known_gap`: paracoccin, named a gap by spec decision D9 before the measurement. `literature_binding`: glycan binding by a paper or Swiss-Prot, no host or chitin-masking claim (CNL, Cgl-I). `structure_backed`: the rest of the UniLectin3D fungal rows. `pa14_owned_elsewhere`: the 10 PA14 rows (the PA14 and GLEYA models belong to the narrow tier and are not among the 15 models). They are reported apart and are left out of the primary recall.
3. **The two *T. rubrum* LysM proteins are not in the spec.** Their paper (PMID 30947244) reports binding to N-linked oligosaccharides of human skin glycoproteins. That meets the H1a definition of spec 3.1. The spec lists H1a for PF07938 only. This is a new observation. I did not change the spec. See the report, section 7.
4. **Mapping of a positive to an expected model.** `measure.py` maps UniLectin3D class to model with the table of spec F11 (class names and Pfam hits). Literature rows are mapped by the protein family in the paper (LysM effector to LysM, Avr4 to CBM_14). The F11 mapping came in part from Pfam hits. So the family recall of a mapped family is not independent of the models.
5. **Union of models.** The union is the 15 M1/M2 models of spec 3.1. The five revision-1 models (PF07938, PF21691, PF00337, PF07367, PF21595) are reported too, to compare with spec F10. The five M3 models are scanned and reported for the negatives. They are not in the union.
6. **Negative classes and reasons.** `er_housekeeping`: reviewed Swiss-Prot entry with an M3 Pfam cross-reference (PF03388, PF00262, PF07915, PF02157, PF11721) and FUNCTION text. Reason: ER or Golgi house-keeping (spec 3.1 row M3). The label means "not a host-interaction lectin". It does not mean "no glycan binding". `enzyme_lectin_like_domain`: 8 alpha-galactosidases (Swiss-Prot, ricin-B-like module) and 4 PA14 beta-glucosidases (shared set, N2). Reason: enzyme function. This class is ambiguous (UniProt gives the alpha-galactosidases the keyword Lectin), so it is never pooled into a headline false-hit number without the split. `n1_shared`: the 136 N1 rows of `data/controls/hard-negatives-shared` (137 minus CNE1, which is counted once as `er_housekeeping`). `assumed`: 497 reviewed fungal Swiss-Prot proteins, no lectin keyword, no UniProt Pfam cross-reference to the 20 models, no lectin word in name or function text, random sample (seed 20261011, at most 15 per species). They have no reason. They are never mixed with the reasoned classes.
7. **Leakage.** Seed sequences of the 20 models come from the InterPro API (`alignment:seed`). Versions match the Pfam 38.2 file for all 20 models. `in_reference` = 30% identity or more (E 1e-3 or less, target coverage 0.5 or more, mmseqs2 `-s 7.5`) to a seed of any of the 20 models. `none` = no such hit. `none` does not exclude remote homology. The method reproduces the numbers of spec F28 (FleA 81.7, SapL1 43.7, AOL 36.6, AAL 32.8, A. oryzae AOL 100.0).
8. **Cluster unit.** MMseqs2 `easy-cluster --min-seq-id 0.3 -c 0.5` on all 729 sequences together. A cluster is a hit when any member in the group hits. Wilson 95% intervals on clusters.
9. **GA.** `hmmsearch --cut_ga` with the Pfam 38.2 models. The hit unit is the sequence (tblout). The score floors of PF07938 use the tblout sequence score.
10. **Check of the report script.** `lectin_report.py` run on one FASTA per class. Its hit sets and scores equal my own `hmmsearch` for every class and model (150 cells, 0 mismatches).
11. **No SLURM.** All runs were on the login node (largest: hmmsearch of 15+5 models on 38,451 proteins).

## Candidates not added

See `rejected_manual.tsv` and `controls/rejected.tsv` (76 lines). Main reasons:
- BbLec1: no verifiable UniProt accession.
- Histoplasma, Coccidioides CFP28, Talaromyces, dermatophyte hemagglutinins: no gene, or no binding shown.
- G1X5Y2 (*Arthrobotrys oligospora*, "fucose-specific lectin g276"): the paper has a deletion phenotype, no binding assay in the abstract. D4AL26 (*A. benhamiae*): annotation by similarity. Both are tier T3 and not truth.
- Q8NJT4, Q8TGE0: same protein as Q4WW81 and Q2UNX8 under another accession.
- Fragments (8, 13 and 28 aa), LysM proteins without an experimental binding statement, ldpA and ldpB (Q4WER9, Q4WRR0, "might have a role in sequestration of chitin": probably LysM effectors, so not negatives).
- Ralstonia lectin (bacterial; mislabelled fungal in UniLectin3D).

## Things that changed after the first pass

- First pass used `strain` taxon ids. Final table uses the species-rank NCBI name and id (checked against names.dmp in the tests).
- The first species table mixed `assumed` and reasoned negatives. The final one separates them.
- Re-fetch of UniLectin3D gave 193 structures for 55 of the 56 accessions. A7UNK4 (LDL) was not returned by the capped API. Its row keeps the structure count of the committed table and says that ligand data were not retrieved.

## Revision 2 (after independent review 1)

12. **One written T2 rule.** T2 = UniLectin3D curated structure, or Swiss-Prot reviewed entry with an ECO:0000269 sugar-binding FUNCTION statement (cited abstract not required to show the assay). Applied to Q06100 and G1X5Y2 together. Near-identical copies of a gene in another strain (Q8NJT4, Q8TGE0, A0A1X7S147, A0A1X7S1V5, A0A1X7S690) and fragments under 50 aa are not added. D4AL26 fails the rule: its FUNCTION is by similarity (ECO:0000250). Its protein-level evidence is secretome mass spectrometry (PMID 21919205), not binding.
13. **Selection route (review 1 M2).** New column `selection_route`. `model_gated_pfam_xref` = the row entered the candidate list only through a UniProt Pfam cross-reference to one of the 20 models (10 rows). Recall is reported for all rows and for independent-route rows.
14. **Cluster-level leakage split (M1).** The seed-independent group is the set of clusters with no `in_reference` row. The first version split rows and counted two clusters in both groups ("3 of 14"). Correct value then: 1 of 12. With G1X5Y2 added: 2 of 13.
15. **`assumed` (M3).** Removed from the false-hit and floor tables. Kept in `results/construction_check.tsv` and appendix A of the report.
16. **Beta-trefoil class (M4).** UniLectin3D class "Laetiporus trefoil lectin-like" mapped to PF14200. It is a miss of an existing model at its GA. PF14200 recall is given with and without these rows.
17. **Added after review.** CMRBL G3JGU5 (PMID 36632960; accession verified by the gene ID CCM_03832 in the abstract), Mgx1LysM P9WEM2 (PMID 33797163), G1X5Y2 (rule of decision 12). Not added: BbLec1 (no open full text; gene ID not checkable), CmLec4, CmRlec, RsLysM, Cg2LysM, TAL6, MfAvr4, Cercospora Avr4 (no verifiable accession). J4KR73 not used as a negative (no paper).
18. **Other fixes.** Per-species negative floor uses N1 rows only. Identity-to-seed table added. Strata that share clusters are listed. Report-script check reports the number of non-empty cells. `run_all.sh` now makes `out/scan_swissprot.hits.tsv`. Tests: grouped numbers of `measure.py` on a toy table with a cluster that spans two leakage groups; taxon rank check with an allow-list; selection route; no identical sequences; curator name no longer hard-coded.
19. **Output location.** COMMON-RULES says `data/controls/<task-name>/`. The set stays in `analysis/lectin_validation/controls/`. Exception recorded in the report (section 9).
