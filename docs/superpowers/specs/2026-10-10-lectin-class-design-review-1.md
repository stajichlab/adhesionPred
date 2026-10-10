# Lectin class design: review 1

*2026-10-10. Independent review of `docs/superpowers/specs/2026-10-10-lectin-class-design.md` (revision 1, commit `d0b8465` on branch `lectin-class-spec`). Reviewer: Claude Code (claude-opus-5-5). Scope: facts, biology, validation design, completeness against task 09, work plan. The spec is not edited. Line numbers refer to the spec file at `d0b8465`.*

## Verdict

**Pass with changes.** The numbers are accurate. I reproduced every number I checked (section 1). The main conclusion holds: the status is `smoke` at most and is `unvalidated` in most species. The structure of the design (M1/M2/M3, ungated evidence call plus gated composite, report-only M2) is sound.

The spec must not go to the owner as it is. The reason is the biology and literature, not the arithmetic. Three findings change the owner questions: the FleA role argument, the missed *Paracoccidioides* lectin, and the LysM/CVNH placement. A fourth finding changes the truth-set statement: M1 hits only 7 of the 25 pooled clusters, and the spec does not say so. Revision 2 must fix the items in section 4.

## 1. Fact checks (what I re-ran or re-read)

| Fact | Check | Result |
|---|---|---|
| F1 | `docs/reports/data/sorting_hat/family_tests/family_tests.tsv` rows PF10528, PF07691 | Counts match (12/54, 1/41, 3/14, 0, 1, 42, 9; 9/54, 1/41, 8/14, 0, 4, 128, 32). Note: that file says GLEYA `active=no`. `data/sorting_hat/family_table.tsv` says `active=yes` (2026-10-09). The spec's "active" claim is correct, but it cites the wrong file for it. |
| F2 | `_workdir/family_tests/lectin_*.domtbl` | 16 GLEYA and 18 PA14 protein-model pairs on positives. Hard negatives: G2WSY1 (GLEYA) and Q4WA69, Q4WL79, Q4WLY1, Q4WU49 (PA14). Match. |
| F3 | `data/controls/hard-negatives-shared/controls.tsv` | Five rows are N2. Glucosidases are HN0069, G2WSY1 is HN0018, CNE1 J9VLH0 is N1 (HN0033). 137 N1 rows in 76 clusters. N1 clusters per species: *S. cerevisiae* 16, *A. fumigatus* 15. Match. |
| F4 | `_workdir/discovery/pfam_names.tsv`, case-insensitive grep | 74. Match. |
| F5 | Independent recount over all 653 `pools.tsv.gz` files with `analysis/fungi_scan/selection_discovery.tsv` for groups (my own script) | Pool totals: S 436,120 (U 422,071 + S-only 14,049). Every protein and genus count in F5 matches, including LysM 1,797 / 1,698 and 26 / 7 / 478. M1 union: 68 in pool S, 63 in pool U, no overlap between M1 models. Match. |
| F6 | `docs/reports/data/sorting_hat/secretome_discovery/families_U_onygenales.tsv` | PF21691 6/27 q 0.01347; PF07938 6/27 q 0.01628; PF18647 15/27 q 0.02608. Match. |
| F7 | `GA` lines in `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm` | All 16 values match. Extra fact: PF07938 Fungal_lectin has `NSEQ 3` (three seed sequences). PF21691 LDL has 14, PF07367 FB_lectin 12, PF21595 CCL2-like 17. |
| F9 | `analysis/lectin_class_spec/unilectin_fungal_sequences.tsv` | 57 rows, sum of `n_pdb_structures` 197, 19 classes, 36 species names. Match. |
| F10 | Re-fetched the 56 UniProt sequences (rest.uniprot.org, 2026-10-10) and re-ran `mmseqs easy-cluster --min-seq-id 0.3 -c 0.5` (MMseqs2 17-b804f) | 33 clusters. The partition is identical to the `cluster30` column. After the three exclusions: 43 sequences, 25 clusters. Ascomycota: 12 sequences, 8 clusters. Match. |
| F11 | Re-ran `hmmscan --cut_ga` (HMMER 3.4, Pfam 38.2) on the same 56 sequences | Identical `pfam_ga_hits` for all 56 rows; 46 with a hit, 10 without. Match. |
| F12 | UniProt TSV `ft_signal` for the 56 accessions | Match. PA14: 9 of 10 (A0A1A9TAD0, 220 aa, has none). No AAL-like entry has a Signal feature. |
| F13 | `data/curated/adhesins/adhesins.tsv` line 191 | Q4WW81, `indirect_regulator`, EXP, PMID 27058347, Pfam PF07938. Match. |
| F14-F18 PMIDs | PubMed metadata for 27058347, 29370403, 33885578, 21149680, 22349222, 26105055, 18331471, 34200153, 33174598, 30239928, 23853587, 26396244, 27546806, 22341093, 34566931, 25504801, 42171338, 32358020, 29610256, 24391211, 24379281, 35259948 | All 22 PMIDs exist. Title, journal, year and DOI match the claimed paper. The abstract statements in F14 match the abstracts. |
| F16 discrepancy | Lebreton 2021 full text (PMC8227253) | **Resolved.** Section 3.3 of the paper: "The MycoLec database contains 33,518 putative lectin sequences (28,691 with filter at 0.25 quality score)". The abstract gives 33,485. The task file number (28,691) is the filtered count. Section 3.4 of the paper names calnexin, calreticulin, the P-type lectin and p58/ERGIC-53 as house-keeping. This confirms the M3 statement at line 61-62, which the spec marks as unread. |
| PMID 25504801 | PubMed | van Eerde 2014, Glycobiology, "Atomic-resolution structure of the α-galactosyl binding *Lyophyllum decastes* lectin ...". The claim at line 49 is correct. |

Not re-checked: F8 (UniLectin licence page and API behaviour) and the "157 of 1,740" count from the UniLectin3D 2019 full text.

Conclusion: no numeric error found. The fact table is reliable.

## 2. Findings, ranked by severity

### High

**H1. The argument for "FleA is not an adhesin" is not valid, and the spec leaves out papers that disagree.** Lines 98, 239, 249 base the label on the knockout outcome: "both knockout studies found worse disease without it". A disease outcome does not show whether a protein mediates attachment. Kerr 2016 (PMID 27058347), which the spec cites, shows that FleA and conidia bind lung mucin and macrophage surface proteins in a fucose-dependent way. That is host-glycan binding by a conidial lectin. Literature that the spec does not cite:

- Houser 2013, PLoS One 8(12):e83077, PMID 24340081, [doi:10.1371/journal.pone.0083077](https://doi.org/10.1371/journal.pone.0083077). AFL (= FleA) structure with a fucoside, glycan array, antibodies "suggested the presence of AFL on the Aspergillus' conidia", IL-8 induction in bronchial cells. The authors conclude that AFL "might be an important virulence factor". This is the opposite reading to Kerr 2016.
- Martínez-Alarcón 2021, Sci Rep 11:16109, PMID 34373510, [doi:10.1038/s41598-021-95008-4](https://doi.org/10.1038/s41598-021-95008-4). SapL1 (*Scedosporium*, the A0A084FYP2 row of the truth table) binds bronchial epithelial cells in a fucose-dependent manner. This is a second host-binding AAL-like lectin. The spec says FleA is "the one protein with host-binding evidence" (line 81). That is not correct.
- Duca 2022, ACS Chem Biol 17(12):3515-3526, PMID 36414265, [doi:10.1021/acschembio.2c00708](https://doi.org/10.1021/acschembio.2c00708). Multivalent fucosides "reduced the adhesion of [*A. fumigatus*] spores to lung epithelial cells". The paper frames FleA and SapL1 as anti-adhesion targets.
- Balogh 2003, Fungal Genet Biol 39(2):128-35, PMID 12781671, [doi:10.1016/s1087-1845(03)00023-9](https://doi.org/10.1016/s1087-1845(03)00023-9). Deletion of AOL (*Arthrobotrys*, G1XA82 row) gave no phenotype, including nematode pathogenicity.

The host-interaction tier is still a defensible place for FleA. The reason must change: the papers show host-glycan binding and immune recognition; the attachment role is disputed. D7 and owner question 4 must state the disagreement and not decide it from disease direction.

**H2. A fungal-encoded lectin in Onygenales is missed.** Lines 7-12, 155 and 217 say there is no fungal lectin literature for the endemic fungi, and the 4.4 table gives "Onygenales 0". The searches were limited to *Coccidioides*, *Histoplasma* and *Blastomyces*. *Paracoccidioides* is in the same family as *Histoplasma* and *Blastomyces* (Ajellomycetaceae, Onygenales). Its paracoccin is a GlcNAc-binding lectin with NAGase/chitinase activity. Recombinant paracoccin binds N-glycans of TLR2 and TLR4, and the binding depends on carbohydrate recognition:

- Alegre-Maller 2014, PLoS Negl Trop Dis 8(12):e3317, PMID 25474158, [doi:10.1371/journal.pntd.0003317](https://doi.org/10.1371/journal.pntd.0003317).
- Oliveira 2017, PLoS One 12(8):e0184010, PMID 28846733, [doi:10.1371/journal.pone.0184010](https://doi.org/10.1371/journal.pone.0184010) (cell-wall and hyphal-tip localisation).
- Pitangui 2021 review, Front Mol Biosci 8:700797, PMID 34532342, [doi:10.3389/fmolb.2021.700797](https://doi.org/10.3389/fmolb.2021.700797); Pitangui 2023, Virulence 14:2150455, PMID 36599817, [doi:10.1080/21505594.2022.2150455](https://doi.org/10.1080/21505594.2022.2150455).

Paracoccin probably has no M1 Pfam model; the reviews describe a chitinase region. I did not check its Pfam hits, so that is an unverified expectation. The spec must add it: as a host-interaction lectin in Onygenales with a published host-glycan binding result, and as an example of a lectin that the Pfam-defined class cannot see (as BAD-1 at line 274). Statement 4 of section 1 and the Onygenales row of 4.4 must change from "none" to "one, in *Paracoccidioides*, not reachable by M1".

**H3. The M1 rule is applied unevenly. LysM and CVNH have the strongest pathogen host-interaction evidence among the 56 structures, but they are in M2.** The tier is "host interaction and immune evasion". The best-documented host-interaction function of a fungal carbohydrate-binding domain in Ascomycota is the LysM effector:

- de Jonge 2010, Science 329:953-5, PMID 20724636, [doi:10.1126/science.1190859](https://doi.org/10.1126/science.1190859). Ecp6 (B3VBK9, a truth-table row) "sequesters chitin oligosaccharides ... to prevent elicitation of host immunity".
- Bolton 2008, Mol Microbiol 69(1):119-36, PMID 18452583, [doi:10.1111/j.1365-2958.2008.06270.x](https://doi.org/10.1111/j.1365-2958.2008.06270.x). Ecp6 is a virulence factor; secreted (UniProt B3VBK9, ECO:0000269).
- Koharudin 2011, Structure 19(5):662-74, PMID 21565701, [doi:10.1016/j.str.2011.03.004](https://doi.org/10.1016/j.str.2011.03.004). MoCVNH (G4N991, a truth-table row) is a CVNH-LysM lectin with measured carbohydrate binding, located in the appressorium during rice infection.

In the spec, line 60 says LysM is "chitin binding, not a glycan-lectin in the narrow sense". Chitin is a glycan, and UniLectin3D lists LysM-like as a lectin class. Line 59 says CVNH "function in fungi not shown". The MoCVNH paper shows binding and infection-stage location. In contrast, four of the five M1 models (LDL, Gal-bind_lectin, FB_lectin, CCL2-like) are mushroom lectins with "host relevance not shown" (lines 49-52). So M1 is selected by "has a fungal structure and is rare", not by host relevance.

Exclusion of LysM from the v1 call can still be right, because LysM is in 511 of 653 genera (my recount agrees). But the spec must state the reason correctly (ubiquity; LysM also occurs in chitinases), give the effector evidence, and offer a defined sub-call as a D2 option. An example: secreted, LysM-only architecture, no catalytic domain. Owner question 1 says "every model has fungal structure evidence". It must also say that host evidence exists for only one M1 class (AAL-like) and for two M2 classes (LysM, CVNH).

**H4. The spec does not state how much of the pooled truth M1 can recall.** Line 162 says the 25-cluster table "measures how many fungal lectin classes the M1 models cover", but it gives no number. I computed it from the committed table: M1 hits **16 of 43 sequences, in 7 of 25 clusters** (AAL-like 2 fungal clusters plus Q2UNX8, LDL 1, galectin 1, FB_lectin 1, CCL2 1; *Ralstonia* excluded). In Ascomycota, M1 hits 3 of 8 clusters (all AAL-like). Every number reported from this table will be close to 0.28 at the cluster level, and the hits are `in_reference`-type (H5). The spec must give this number in 4.4 and in owner question 3.

### Medium

**M1. The only species status that can be measured is sensitivity on one in-reference protein.** Line 160 says "the status rules decide whether one cluster gives `smoke`; this spec does not assume it". The code decides it: `status_from_measure` in `src/cellsurface_sorting_hat/status.py:105-133` returns `smoke` for `n_pos >= 1` with a sensitivity. Leakage other than `none` also caps at `smoke` (`call_status.py:99`). So *A. fumigatus* gets `smoke` from FleA, 1 of 1. PF07938 has 3 seed sequences (F7 check) and its Pfam/InterPro description names AAL. Seed membership of FleA is not checked. The spec must state the consequence plainly: the *A. fumigatus* `smoke` number is 1/1 on a protein that is probably in the model's reference set, so it carries no information about sensitivity. Option for the owner: do not write a status file for this call until there is a positive that is not structure-derived.

**M2. No negative tests the specificity of M1.** Line 143 says this for the PA14 glucosidases. The same is true for every negative in 4.3. CNE1 (calnexin), Emp46/47 (Lectin_leg-like) and the glucosidases hit models that are not in M1. The shared N1 set was scanned only with the 11-model `lectins.hmm` (F2), which does not contain PF00337, PF07367 or PF21595. A specificity on the 15 *A. fumigatus* N1 clusters will be near 1.0 for a reason that says nothing about M1. The spec must say this in 4.3 and 4.4. It must also add a task to scan the N1 set with all five M1 models.

**M3. Emp46/47 have a paper; the spec says none was read.** Lines 142, 235 and 271 say the localisation needs a paper. The UniProt record of Q12396 cites PMID 12134087: Sato and Nakano 2002, Mol Biol Cell 13(7):2518-32, [doi:10.1091/mbc.e02-01-0027](https://doi.org/10.1091/mbc.e02-01-0027). The abstract states Golgi localisation, ER-Golgi cycling, and a role "for the export of specific glycoprotein cargo from the endoplasmic reticulum". This supports the M3 exclusion and N1 status. Both proteins have a signal peptide (UniProt, experimental evidence) and a TM helix. They would pass the R0 gate. They are relevant negatives for `secreted_lectin_candidate` only if Lectin_leg-like ever joins the call.

**M4. "FleA-orthology" for the Onygenales PF07938 hits is not supported.** Line 223 calls a Pfam GA hit a "FleA-orthology ... sequence statement". A Pfam hit is family membership, not orthology. My recount shows something that the spec does not report: PF07938 has **0 of 7 Aspergillaceae genera** in the pools, while 6 of 27 Onygenales genera have it. The pools need a signal peptide. FleA has no UniProt signal peptide (F12). So the Onygenales PF07938 proteins in pool U carry a predicted signal peptide, unlike FleA. They may be a different subgroup. Risk 3 (line 256) touches this. Section 6.3 must say it and must drop "orthology".

**M5. D5 departs from precedent without saying so.** `docs/CLASSES.md:17` shows that `surface_attachment_candidate` already reads hydrophobin and HsbA. Neither is an attachment mechanism in a strict sense. The spec cites the hydrophobin precedent for D2 (line 91) but chooses the opposite for D5. The default may be right. The spec must state the inconsistency so the owner decides it with that knowledge.

**M6. The ungated call is defensible, but the spec does not say what it lets in.** The reasoning at line 81 is sound. Many fungal lectins leave the cell without a classical signal peptide. Data in the truth table supports this. CGL2 (Q9P4R8) has no UniProt Signal feature, but UniProt records it as secreted to the extracellular matrix and the cell wall (ECO:0000269, PMID 10931889; I did not read that paper). The cost: `lectin_domain` will call cytoplasmic lectins in a tool for cell-surface proteins. In mushrooms, that includes fruit-body defence lectins. The report text at line 212 must say that the ungated call does not imply surface location.

### Low

**L1.** F1 cites `family_tests.tsv` for "active", but that file has `active=no` for GLEYA. Cite `data/sorting_hat/family_table.tsv` (GLEYA active since 2026-10-09; PA14 has `second_condition=signal_peptide`).

**L2.** Lines 266 and 272: the Lebreton discrepancy is resolved (section 1). Update F16 and section 10. The house-keeping statement can now cite the paper's section 3.4.

**L3.** Line 106 and L4 (line 287) say that the per-module family digest (H1c) "must exist first". It exists: `_family_digest(families, module)` in `src/cellsurface_sorting_hat/modules/cli.py:225-232` hashes only the rows of one module. Change the dependency to "done".

**L4.** Reproducibility. The recount script for F5 is not committed (line 272). The sequence FASTA, the MMseqs2 output and the domtbl are not committed (line 296). I re-fetched the sequences from UniProt and got the same result, so the table is reproducible today. UniProt can change. Commit `sequences.faa` (UniProt data is CC BY 4.0) and the recount script, as COMMON-RULES asks (`sequences.faa`, `clusters.tsv`).

**L5.** Line 49: LDL is in M1 partly "for the discovery-pass signal". That selects a model by the result it will then be reported beside. Its structure evidence (PMID 25504801) is enough for the rule; remove the discovery reason, or label it as a selection effect.

**L6.** Missed house-keeping class: Lebreton 2021 names the P-type lectin (ERAD) as house-keeping. In yeast that is the MRH-domain protein Yos9 (from memory, not checked). M3 should list it, or say that it was not checked.

**L7.** Language. Line 70 "M1 is small on purpose" and line 246 "every model has fungal structure evidence" read as strengths. Line 248 ("More truth would need new experiments in the endemic fungi") is an opinion without data. With H2, it is also not accurate: paracoccin and SapL1 show that published host-binding lectins exist outside the six Phase C species.

## 3. Completeness against task 09

| Item | Covered | Gap |
|---|---|---|
| 1 Models | yes | H3 (LysM/CVNH evidence), L6 (P-type) |
| 2 Gate | yes | M6 (what the ungated call admits) |
| 3 Truth set | yes; "smoke only" stated | H4 (M1 coverage 7 of 25), M1 (1/1 in-reference), M2 (no M1 specificity negatives) |
| 4 Curation | yes; rubric follows the hydrophobin pattern | Rule 1 (`lectin (T4)` if line C holds) lets one paper with "host binding" make a T4 lectin call, without glycan binding. Say "glycan binding" for C, or record host binding only as the H flag. |
| 5 Relation to calls | yes | M5 (D5 precedent) |
| 6 Measurement, discovery families | yes | M4 (orthology wording, Aspergillaceae 0) |
| 7 Work plan, decisions | yes | see section 5 |

Literature named in the task: FleA (found), Epa/Flo PA14 (found), endemic fungi carbohydrate-binding proteins (partly; H2). Lebreton "references ... read them" (line 20 of the task): the inflammation papers behind Lebreton's sentence were not identified. Houser 2013 (H1) is the obvious candidate. I did not confirm that it is Lebreton's reference.

Other searches I ran (PubMed, 2026-10-10): "Histoplasma capsulatum fungal lectin carbohydrate-binding protein yeast surface" gave 0 hits. That is consistent with F17. I did not search *Candida albicans* beyond Als/Epa or *Coccidioides* further. The absence of a fungal-side lectin paper for those two remains "not found by these searches".

## 4. Required changes before the owner sees the spec

1. Rewrite the FleA role (lines 81, 98, 239, 249; F14). Cite Houser 2013, SapL1 2021 and Duca 2022. State that the attachment role is disputed. Base the host-interaction label on host-glycan binding and immune recognition, not on the direction of the disease outcome. Remove "the one protein with host-binding evidence".
2. Add paracoccin (*Paracoccidioides*, Onygenales) to F17, section 1 point 4, the 4.4 Onygenales row and section 6.2 item 6. State whether it has an M1 hit, after a check, or mark it unchecked.
3. Restate the M1 rule and the M2 rows for LysM and CVNH with the effector evidence (PMIDs 20724636, 18452583, 21565701). Keep or change the default, but give ubiquity as the reason. Add to D2 the option of an architecture-defined LysM-effector sub-call as future work.
4. Add the M1 coverage number to 4.4 and owner question 3: 16 of 43 sequences, 7 of 25 clusters; Ascomycota 3 of 8.
5. State in 4.4 statement 2 that the code gives `smoke` for one positive, and that the *A. fumigatus* number is 1/1 on a probable in-reference protein (PF07938 NSEQ 3). Give the owner the option of no status file.
6. State in 4.3 that no current negative tests M1 specificity. Add a task (L3) to scan the N1 set with all five M1 models.
7. Cite PMID 12134087 for Emp46/47 and remove "no paper read".
8. Replace "FleA-orthology" (line 223). Report PF07938 Aspergillaceae 0 of 7 in the pools, and that the pool members have a signal peptide, unlike FleA.
9. State the D5 departure from the hydrophobin/HsbA precedent in `surface_attachment_candidate`.
10. Report text: the ungated `lectin_domain` does not imply surface location.
11. Low items L1 to L7.

## 5. Work plan and owner-decision defaults

- The order is realistic. L1 to L3 are small and need no SLURM. L2 (SignalP 6 on 56 sequences) is small. Running it before D6 is final is correct.
- Remove the H1c dependency from L4 (L3 above).
- L3b must include the papers of this review: Houser 2013, SapL1, Duca 2022, paracoccin, LysM effectors, MoCVNH. Add an explicit stop if L3b finds a host-binding lectin family that M1 cannot detect. Paracoccin already meets that condition.
- Add a task: scan the shared N1 set and the six Phase C proteomes with the five M1 models. It gives the real counts of positives and look-alikes per species before L6.
- Defaults: D1 (M1 inactive until review) is safe. D2 (M2 report only) is acceptable if H3 is stated. D3 is supported, now with a paper for Emp46/47. D4 follows the hydrophobin D8 precedent. D6 is sound. D7 must be rewritten (H1). D8, D9, D10 are acceptable.
- None of the defaults is data-driven about performance. The spec says this for D1 (line 233). It should say it for D2 and D7 as well.

## 6. What I did not do

- I did not read full texts other than Lebreton 2021. All other papers were read at abstract level.
- I did not check Pfam seed membership for PF07938 (no `Pfam-A.seed` in `/srv/projects/db/pfam/2026-01-27-Pfam38.2/`).
- I did not run SignalP, and I did not scan proteomes with the M1 models.
- I did not check the UniLectin licence page or the API.
- My recount, clustering and hmmscan inputs and outputs were in a scratch folder in the worktree and are not committed. Commands: `mmseqs easy-cluster seq.faa clu tmp --min-seq-id 0.3 -c 0.5`; `hmmscan --cut_ga --cpu 4 --domtblout scan.domtbl Pfam-A.hmm seq.faa`; the sequences came from `rest.uniprot.org/uniprotkb/accessions?accessions=<56 IDs>&format=fasta`. No SLURM jobs were run.
