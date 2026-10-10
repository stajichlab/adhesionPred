# Lectin class design: review 2

*2026-10-10. Independent review of `docs/superpowers/specs/2026-10-10-lectin-class-design.md`, revision 2 (commit `76dc2f9`, branch `lectin-class-spec-r2`). Reviewer: Claude Code (claude-opus-5-5). Line numbers refer to the spec file at `76dc2f9`. The spec is not edited. Review files are in `analysis/lectin_class_spec/review2/`. No SLURM jobs were run.*

## Verdict

**Fail. The spec needs revision 3 before the owner sees it.**

Most of the 11 changes of review 1 are fixed in the files, not only in words (section 1). The numbers that I re-ran are correct (section 2).

One new finding breaks the main design decision. The ubiquity value U is measured on secretome pools (proteins with a signal peptide). The call `lectin_domain` is ungated: it reads the whole proteome. On the whole proteomes of the same 653 genera, Fungal_lectin PF07938 hits **392 of 653 genera (60.0%)**, not 28 (4%). Under the spec's own rule (common at 10% or more), PF07938 goes to M2 and M1 is empty. If PF07938 stays in an ungated call, the call fires in 60% of genera. That is the flooding that the rule was made to prevent. D1, owner question 1 and section 4.4 depend on this.

## 1. Review 1 changes: fixed in the files?

| # | Review 1 required change | Result | Evidence |
|---|---|---|---|
| 1 | FleA role | **Fixed.** F14 lists eight papers. It separates binding, disease outcome and inhibitor results. "The one protein with host-binding evidence" is gone. D7 option A does not decide the role. | spec:53, 282. I confirmed the eight PMIDs (section 2). Dussouy 2021 says "controversial role" in the abstract. |
| 2 | Paracoccin | **Fixed in part.** F20 added; hmmscan run; section 1, 4.4, 6.2 changed. But the founding paper (Coltri 2006) is missing, and the C1G842 mapping has a size conflict. See finding H2. | spec:59, 194, 259 |
| 3 | LysM/CVNH placement | **Fixed in words and rule.** One rule (H and U) is applied to every model in the table. LysM is excluded by ubiquity, with Ecp6 evidence. The D2 sub-call option is present. But U is measured on the wrong population (finding H1). | spec:68-100 |
| 4 | M1 coverage | **Fixed.** 16 of 43 sequences, 7 of 25 clusters; PF07938 5 of 43, 3 of 25; Ascomycota 5 of 12, 3 of 8. I reproduced all six numbers. The committed `cluster_check.py` prints only the five-model numbers. The claim at spec:12 that it reproduces "both counts" is not exact (finding L1). | `review2/cov.py` |
| 5 | `smoke` from one in-reference positive | **Fixed.** Statement 2 reads the code correctly (`status.py:105-133`, `call_status.py:96-100`). D4b offers no status file. | spec:199, 279 |
| 6 | No negative tests M1 | **Fixed.** N1 scan done (0 hits). Spec states that this says little. L3 asks for look-alike negatives. | spec:181, `m1_scan_hits_per_model.tsv` (N1 column: only Calreticulin 1) |
| 7 | Emp46/47 paper | **Fixed.** PMID 12134087 title and abstract match. | spec:61 |
| 8 | "FleA-orthology" | **Fixed.** Section 6.3 says family membership and reports 0 of 7 Aspergillaceae. | spec:265 |
| 9 | D5 precedent | **Fixed.** Both options are stated. | spec:127, 280 |
| 10 | Ungated call text | **Fixed** in 6.2 item 3. But D7 option A contradicts it (finding H3). | spec:256 vs 282 |
| 11 | L1-L7 | **Fixed.** Scripts, FASTA and clusters are committed. Malectin still says "from memory" (honest). | spec:92-93, 341-355 |

Recorded disagreement: the author applied the review's rule to every model instead of keeping five models. I agree with the method. Consistent application was what review 1 asked for. The result is wrong only because of how U was measured (H1).

## 2. Checks I ran

| Claim | Check | Result |
|---|---|---|
| F10 clusters | `mmseqs easy-cluster --min-seq-id 0.3 -c 0.5` (MMseqs2 17-b804f) on the committed FASTA; `cluster_check.py` | 56 sequences, 33 clusters, same partition as `cluster30`; after exclusions 43 / 25; five models 16 / 7. My cluster file is byte-identical to `unilectin_fungal_clusters.tsv`. **Match.** |
| F10 PF07938 and Ascomycota | `review2/cov.py` | PF07938 5 of 43 sequences, 3 of 25 clusters (clusters P18891 with SapL1; Q2UNX8 with FleA; G1XA82). Ascomycota 12 sequences, 8 clusters; PF07938 5 / 3. Most clusters in one species: *Coprinopsis* 3. **Match.** |
| F20 paracoccin Pfam | Fetched C1G842 from UniProt; `hmmscan --cut_ga` and `-E 1 --domE 1`, Pfam 38.2, HMMER 3.4 | One hit: Glyco_hydro_18 PF00704, domain score 35.8, alignment 31-176. No other hit even at E 1. **Match.** New: the hit covers **HMM positions 2-115 of 317**. C1G842 is 218 aa, TrEMBL, sequence version 2, PE=3. |
| F23 scan counts | Read `m1_scan_hits.tsv`, `m1_scan_hits_per_model.tsv` | N1: 0 M1-set and 0 M2 hits; 1 Calreticulin (CNE1). *A. fumigatus*: PF07938 Q4WW81 (489.1) and Q4WMH4 (30.2); FB_lectin Q4WAA0. Per-model counts as stated. **Match.** I did not re-run the six-proteome scan. |
| Q4WMH4, Q4WAA0 | UniProt REST | Q4WMH4: 108 aa, "Fucose-specific lectin", no signal peptide (a fragment-length AAL-like protein; FleA is 315 aa). Q4WAA0: 128 aa, uncharacterised, no signal peptide. |
| M3 PF07915 | F23 scan rows | PF07915 hits OS9_ASPFU (Q4WCG2) and OS9_EMENI (Q5BDB9), Swiss-Prot OS-9 homologs. Supports the M3 statement. |
| **U on whole proteomes (new)** | `review2/scan653.sh`: `hmmsearch --cut_ga` with PF07938, PF07367, PF00337, PF21595, PF21691 on all 653 discovery proteomes (`selection_discovery.tsv`, 6,847,098 proteins); `sum653.py`, `score653.py` | See table below. |

Whole-proteome scan (Pfam 38.2, GA cutoff, no signal-peptide condition), against pool S of F5:

| Model | Whole proteome: proteins, genera (share of 653) | Pool S (F5): proteins, genera |
|---|---|---|
| Fungal_lectin PF07938 | 997; **392 (60.0%)**; Onygenales 11 of 27, Aspergillaceae 4 of 7 | 32; 28 (4%) |
| FB_lectin PF07367 | 75; 68 (10.4%) | 0; 0 |
| LDL PF21691 | 39; 36 (5.5%) | 34; 31 |
| Gal-bind_lectin PF00337 | 41; 31 (4.7%) | 2; 2 |
| CCL2-like PF21595 | 2; 2 (0.3%) | 0; 0 |

PF07938 by sequence score: at 50 or more, 188 genera (28.8%); at 100 or more, 104 (15.9%); at 200 or more, 52 (8.0%); at 300 or more, 40 (6.1%). FleA scores 489. Output: `review2/whole_proteome_hits.tsv`, `review2/whole_proteome_genera.txt`.

PMIDs. I fetched PubMed metadata for 22 PMIDs that are new in revision 2 or that carry new claims: 24340081, 30356167, 36414265, 34373510, 33885578, 12781671, 25474158, 28846733, 34532342, 36599817, 20724636, 18452583, 21565701, 26455798, 34649678, 12134087, 30425696, 27406561, 30696749, 30688363, 29370403, 16476564 (the last is new, from me). All exist. Title, journal, year, volume, pages and DOI match the spec's citation for all 21 that the spec cites. I compared the spec's claims with the abstracts of 12 of them (30356167, 36414265, 34373510, 33885578, 25474158, 28846733, 26455798, 34649678, 30425696, 30688363, 18452583, 24340081). All claims match the abstracts. I read the full text of Pitangui 2021 (PMC8438136). Its statements in F20 match, with one exception (H2: the paper gives no protein length).

## 3. Findings, ranked by severity

### High

**H1. U is measured on a population that the call does not read. With U measured on whole proteomes, the stated rule gives an empty M1.**

- The rule (spec:71-73) defines U as "the share of the 653 genera with at least one **pool S** protein". Pool S contains only proteins with a predicted signal peptide.
- The call `lectin_domain` is ungated (spec:108, 113). It reads every protein of a proteome.
- The spec knows that the two populations differ for this family: FleA has no signal peptide (F12), and PF07938 is in 0 of 7 Aspergillaceae pool S sets (spec:265), although *A. fumigatus* has FleA.
- My scan: PF07938 is in 392 of 653 genera at GA (section 2). By the rule's own threshold (10%), PF07938 is common. M1 is then empty. FB_lectin is also above 10% (10.4%).
- The default of D1 (spec:275) and the answer to owner question 1 (spec:289) rest on "rare (28, 4%)" (spec:77). That value describes about 3% of the PF07938 proteins (32 of 997).
- Risk 10 (spec:307) says "at 1% PF07938 becomes common". The real issue is not the threshold value. Any threshold between 4.3% and 60% changes the result.

The rule is still principled if U and the gate match. Revision 3 must choose one of these, or another the author can defend:

- (a) **Gated call only** (`secreted_lectin_candidate`), with U from pool S. PF07938 is rare. But FleA, the protein that gives the H1 evidence, probably drops out (L2 must run first). H1 evidence then comes from SignalP-negative proteins and is applied to SignalP-positive proteins (6.3 already says these may be a different subgroup).
- (b) **Ungated call with a custom score cutoff** (for example full-length propeller hits). This breaks the v1 rule "GA only" (spec:73). A cutoff needs its own truth, which does not exist.
- (c) **No lectin call in v1.** All lectin-like models are report-only (M2), with H and whole-proteome U in the report. The call waits for L2 and a non-structure positive.

My view: (c) is the most defensible default with the present data, and (a) is second after L2. This is a judgement, not a measured result. The owner decides.

**H2. Paracoccin: the founding paper is missing, the C1G842 identity is not settled, and the "lectin sites" lie inside the GH18 domain.**

- Coltri 2006, Microbes Infect 8(3):704-13, PMID 16476564, [doi:10.1016/j.micinf.2005.09.008](https://doi.org/10.1016/j.micinf.2005.09.008), abstract: paracoccin was isolated as "an additional 70-kDa minor protein" that binds GlcNAc; it "was able to interact with laminin in a dose-dependent manner", and GlcNAc inhibited this. This is host-ECM glycoprotein binding by a fungal surface protein. The spec does not cite it. It is the closest published case of a fungal lectin with an adhesion-relevant host target in Onygenales.
- Size conflict. The native protein is 70 kDa (Coltri 2006). C1G842 is 218 aa (about 24 kDa unglycosylated). The spec says C1G842 "matches the paper's length" (spec:59). Pitangui 2021 (full text) gives no length. It gives a signal peptide at 22/23 and a GH18 domain at 29-218. Glycosylation or a different gene model could explain the gap. I did not check which. The spec must state the conflict.
- Partial GH18. The Pfam hit covers HMM positions 2-115 of 317. A 218-aa protein is short for a complete GH18 domain. The gene model may be truncated (TrEMBL, sequence version 2). Not checked.
- Lectin sites. Pitangui 2021 (full text) places the four predicted carbohydrate-binding regions at residues 30-37, 55-62, 98-107 and 120-126, all inside the GH18 region (29-218). The 3D model uses a GH18 chitinase (PDB 2XTK) as template. So the published evidence does not separate a lectin fold from GH18 substrate binding. The spec says "I did not check whether GlcNAc binding is a property of the GH18 fold" (spec:59, 312). The paper's own coordinates answer part of that question. Revision 3 should report this.
- Consequence for D9 (spec:284). A custom model for paracoccin would in practice be a GH18 chitinase model. GH18 chitinases are common in fungi. Such a model would flood the call and would not be a lectin model. D9 should say this.
- Pitangui 2021 also states that "PCN-silencing does not affect the adherence or phagocytosis" of yeasts. Add it beside the laminin result.

**H3. D7 option A contradicts the fixed report text.**

- spec:256 (fixed text for every `lectin_domain` row): "It does not show that the protein binds a host glycan."
- spec:282 and 292: label PF07938 proteins "host-glycan-binding lectin, role disputed".
- spec:139: family-table class text "host-interaction lectin".

The H1 evidence holds for FleA and SapL1. It does not transfer to every PF07938 hit. AOL (*Arthrobotrys*) has no phenotype (Balogh 2003). The pool S members differ from FleA (spec:265). With 997 whole-proteome hits (H1), most labelled proteins would have no binding data. The label must name the family, not the function. Example: "AAL-like (PF07938) lectin domain; host-glycan binding shown for FleA and SapL1 only; role disputed".

### Medium

**M1. The model list comes from a name search, so carbohydrate-binding effectors not named "lectin" are outside the design.** F4 (spec:43) selects Pfam models by the words lectin, agglutinin, galectin or glycan binding. The rule of 3.1 is then applied only to those models. Example missed: *Cladosporium fulvum* Avr4. van den Burg 2006, Mol Plant Microbe Interact 19(12):1420-30, PMID 17153926, [doi:10.1094/MPMI-19-1420](https://doi.org/10.1094/MPMI-19-1420), abstract: Avr4 "is a chitin-binding lectin containing an invertebrate chitin-binding domain (CBM14)", binds fungal cell-wall chitin during infection and protects against plant chitinases. Its evidence is as strong as that for MoCVNH. Its Pfam model is probably PF01607 ("Chitin binding Peritrophin-A domain"); I did not scan Avr4. Other chitin- or glucan-binding modules (CBM18 / PF00187 "Chitin recognition protein"; WSC, cited at F21 d) are in the same position. Revision 3 must say that v1 is limited to models with "lectin" names and give the reason, or extend the H/U table.

**M2. H measures how much a family has been studied, not only its biology.** H1 needs a paper. Families from lung pathogens and model plant pathogens get papers. Mushroom families and rare Ascomycota families do not. H0 is "not found", which the spec states (spec:70). But the consequence is a bias: the rule favours families with a famous member. Say this in Risk 10. Also state the mixed meaning of H1: FleA and SapL1 bind host glycans; Ecp6 binds the fungus's own chitin to hide it. Both fit the H1 wording ("an immune effect that depends on glycan binding"), but they are different mechanisms. MoCVNH binds the fungus's own cell-wall chitin too (Koharudin 2015 title: "Fungal Cell Wall Recognition") and gets H2 only because no gene-loss phenotype is published. State that reason in the CVNH row (spec:79).

**M3. Is a one-model call useful?** Even before H1, M1 = PF07938 is a FleA-family detector: 2 proteins in six Phase C proteomes, both *A. fumigatus*, one a 108-aa fragment (Q4WMH4); 3 of 25 structure clusters. D1 option B adds four mushroom models with no host evidence to a host-interaction tier, and FB_lectin is above 10% on whole proteomes. Neither option gives a useful general "lectin" call. If a call is kept, its name should say what it is (for example `aal_like_lectin_domain`), not `lectin_domain`.

**M4. Missing *Aspergillus* papers (title and abstract checked).**
- Kuboi 2013, J Infect Chemother, PMID 23695231, [doi:10.1007/s10156-013-0614-9](https://doi.org/10.1007/s10156-013-0614-9): AfuFleA purified as a 32-kDa L-fucose hemagglutinin, gene cloned. Independent identification of FleA.
- Houser 2015, Acta Cryst D, PMID 25760594, [doi:10.1107/S1399004714026595](https://doi.org/10.1107/S1399004714026595): AFL binding sites are non-equivalent; affinity for oligosaccharides found on human epithelia.
- Lehot 2018, Chemistry, PMID 30277619, [doi:10.1002/chem.201803602](https://doi.org/10.1002/chem.201803602): multivalent fucosides "highly efficient in preventing fungal spore adhesion to bronchoepithelial c[ells]". A second inhibitor study, same type of evidence as Duca 2022.
- Tronchin 2002, Infect Immun, PMID 12438366, [doi:10.1128/IAI.70.12.6891-6895.2002](https://doi.org/10.1128/IAI.70.12.6891-6895.2002): a 32-kDa sialic acid-specific lectin purified from conidia, framed around adherence to extracellular matrix. The abstract does not name a gene. Whether it is FleA or another protein is not known to me.
- Powers-Fletcher 2011, PLoS One, PMID 22163332, [doi:10.1371/journal.pone.0028865](https://doi.org/10.1371/journal.pone.0028865): calnexin ClxA of *A. fumigatus* is an ER lectin chaperone; a paper source for the M3 calnexin row in a Phase C species.

F14's conclusion ("role disputed") does not change with these papers. Lehot 2018 adds weight to the inhibitor line of evidence. It still is an inhibitor result.

**M5. F23 already contains look-alike enzyme candidates that L3 says are "not yet identified" (spec:177).** The *A. fumigatus* scan shows Ricin_B_lectin on Q4WVZ3 (Swiss-Prot AGALA_ASPFU, an alpha-galactosidase) and LysM on Q4WER9 and Q4WRR0 (Swiss-Prot LYSMA, LYSMB). These are M2 models, so they do not test M1. They are ready inputs for L3 if a sub-call is ever specified.

### Low

**L1.** spec:12 says the PF07938-only counts were reproduced "with the committed `cluster_check.py`". The script has the five-model set hard-coded (`M1 = {...}`, line 14) and prints only those counts. The PF07938 numbers are correct (my `cov.py`), but the script does not print them. Add a model-set argument.

**L2.** The F23 proteome FASTAs and the paracoccin domtbl are not committed (spec:355). The scan is cheap to repeat, so this is acceptable. Record the UniProt proteome release date.

**L3.** Section 4.4 heading "the honest answer" is a judgement word. Use "Cluster counts and the expected status".

**L4.** F20 says the C1G842 signal peptide comes "by SignalP in UniProt". The UniProt JSON feature has no evidence code. Say "UniProt Signal feature (evidence not shown)".

## 4. Specific questions of the task

**Q3. The M1 rule.** The rule is applied in the same way to every model in the table. It is not circular with the discovery pass: H comes from papers about *Aspergillus*, *Scedosporium* and *Cladosporium*, and U from all 653 genera; the Onygenales signal did not select PF07938. But the rule has two defects: U is measured on the wrong population (H1), and H depends on research attention (M2). The 10% threshold is not the deciding input. Only two models have H1 (PF07938, LysM). Any threshold between 4.3% and 78% gives the same M1 on pool S, and on whole proteomes any threshold below 60% removes PF07938. A one-model call is of low use (M3). D1 option B is not a better default: it puts models without host evidence in a host-interaction tier, and FB_lectin is common on whole proteomes. I recommend option (c) of H1 as the default for revision 3. This is a judgement without performance data.

**Q4. Validation design.**
- Circularity and leakage: handled correctly (4.1, statement 2, leakage `in_reference`).
- 20-cluster arithmetic: correct. Most positive clusters in one species is 3 (*Coprinopsis*) or 1 (*A. fumigatus*). N1 most per species 16 (*S. cerevisiae*).
- D4b default (no status file): consistent with the code and with `docs/CLASSES.md` (a status holds only where a status file is installed; otherwise `unvalidated`). I agree with the default.
- D5 and double counting: the spec is consistent for D5 = No. If the owner chooses Yes, a secreted PF07938 protein is counted by two composites (`secreted_lectin_candidate` and `surface_attachment_candidate`). The spec should state how the report shows that overlap.
- The N1 scan result (0 hits) is correctly called uninformative.

**Q5. Biology.**
- FleA "role disputed": supported by the papers. Dussouy 2021 uses the word "controversial". Fair.
- Paracoccin: see H2. The handling is not yet correct.
- LysM and CVNH: the ubiquity reason is right. The H-level reason for CVNH needs one sentence (M2).
- Signal peptide claims: correct for UniProt features. SignalP 6 not run (stated).
- M3: PF07915 is supported by the scan (OS-9 hits). Calnexin has a paper (Powers-Fletcher 2011). Malectin is still "from memory"; it has 1 pool S protein, so the risk is small.

**Q6. Literature gaps.** H2 (Coltri 2006), M1 (Avr4), M4 (four *Aspergillus* papers). Claims that rest on abstracts only: all of F14 except Houser 2013; F21 in full. Of these, the strongest claim from an abstract is Koharudin 2015 ("function is associated with protecting M. oryzae from the host immune response"). That is an abstract framing sentence, not a shown result. The spec already marks it as abstract-level and uses it only for H2. Acceptable.

**Q7. Work plan and owner questions.** The tasks are small and realistic. Changes:
- Run L2 (SignalP 6 on the 56 sequences) before D1, not after. D1 depends on whether FleA has an R0 signal peptide (H1 option a).
- Add the whole-proteome U to F5 for every model, with the score distribution of PF07938.
- Owner question 1 must be rewritten after H1. Questions 2 to 6 are acceptable. Question 4 needs the label fix of H3.

## 5. Required changes for revision 3

1. Measure U on the population each call reads. Report whole-proteome U for every model in the 3.1 table. Restate the consequence for M1 and choose a default among options (a), (b), (c) of H1. Rewrite D1, owner question 1, section 1 point 1, 4.4 and Risk 10.
2. Paracoccin: cite Coltri 2006 (laminin binding, 70 kDa). State the size conflict with C1G842 and the partial GH18 match (HMM 2-115 of 317). State that the predicted lectin sites are inside the GH18 region. Add to D9 that a paracoccin model would be a GH18 model. Add the "silencing does not affect adherence" sentence.
3. Remove the conflict between D7 option A and the fixed report text. Label the family, not the function.
4. State that the model list is limited to Pfam models with lectin-type names. Name Avr4 (CBM14) as a known gap, or add CBM14, CBM18 and WSC to the H/U table.
5. Add the H-bias statement and the H1 mixed-mechanism statement (M2).
6. Add Kuboi 2013, Houser 2015, Lehot 2018, Tronchin 2002 and Powers-Fletcher 2011, or state why not.
7. Move L2 before D1 in the work plan.
8. Low items L1 to L4.

## 6. Ready for the owner?

No. Revision 3 is needed. Item 1 changes the main recommendation. Items 2 and 3 change text that the owner will read as fact. After these changes, a short third review of items 1 to 3 is enough. The rest of the spec is in good order.

## 7. What I did not do

- I did not re-run `pool_recount.py` or the six-proteome scan. I read their outputs.
- I did not run SignalP. I did not check Pfam seed membership.
- My whole-proteome scan used five models only (the revision-1 M1 set). I did not measure whole-proteome U for the M2 models. They are already common on pool S, and whole-proteome U can only be equal or higher.
- I did not check whether many of the 997 PF07938 hits are fragments. The score distribution in section 2 is the only check.
- Papers in M4 and H2 were read at abstract level, except Pitangui 2021 (full text).
- My scan commands are in `analysis/lectin_class_spec/review2/scan653.sh` (login node, `hmmsearch --cut_ga --cpu 2`, HMMER 3.4, Pfam 38.2). The script writes to a review scratch directory that is not committed. The combined hits are in `whole_proteome_hits.tsv`.
