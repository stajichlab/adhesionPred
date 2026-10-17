Errata applied 2026-10-17 from reports 05 to 07; items still unverified are marked.

# Cell wall remodeling and cell fusion: what the repo covers and where the gaps are

*2026-10-17. Read-only inspection by Claude Code (claude-sonnet-5-5). Nothing in the repo was changed. This note is one of four in `docs/reports/cell-wall-fusion/`. It describes the repo, not the literature.*

## 0. How to read this note

**Two states of the repo.** `main` and several unmerged remote branches differ. I read both. Every citation says which one.

- `main` (checked out): the family table, `categories.yaml`, the hydrophobin and surface-attachment work (PR #77).
- Unmerged remote branches (read with `git show`, not checked out): `origin/lectin-reviews-and-full-run` (2026-10-17) and `origin/allergen-predictor-wip` (2026-10-16, adds the secretome discovery reports). Cited as `origin/<branch>:path`. These hold `docs/CLASSES.md`, the M2 spec of the v1 class list, the lectin class spec, the lectin report script and the secretome discovery reports. **None of these exist on `main`.** `docs/CLASSES.md` is absent from the working tree (checked with `find`).

**Verification level of Pfam IDs.** I did not run `hmmfetch` or open InterPro [corrected 2026-10-17, see 05: all 53 Pfam IDs in reports 01 to 04 and the family table were later fetched from InterPro and resolve; names match, except that which proteins sit in GH76, GH17, GH16 and GH5 is only partly settled, see 05 Summary item 3 and 07]. "In repo" means the ID and name appear in text I read. The 2026-10-04 report says it checked each accession against a local `Pfam-A.hmm` of 30,134 models, and that it did not check that a family covers every gene named in the papers (`docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md:20-23`). I did not repeat that check. An ID seen only as an annotation on a data row is marked as such. Anything else is marked "unverified: check InterPro/Pfam before use".

**Discrepancy to resolve before any use.** The task brief names `PF03142` as Chitin_synth_1 and `PF01644` as Chitin_synth_2. The repo says the reverse: `PF01644` Chitin_synth_1, `PF03142` Chitin_synth_2, `PF08407` Chitin_synth_1N (`docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md:56`). The orchestrator spec lists the two IDs without names (`docs/superpowers/specs/2026-10-04-orchestrator-design.md:305`). Check both against Pfam.

## 1. What the repo already covers

### 1.1 Sorting-hat family table and calls (these produce calls on a proteome)

The family table is `data/sorting_hat/family_table.tsv` (24 rows plus header). A row is read only if `active=yes` (`src/cellsurface_sorting_hat/modules/pfam.py:135`). `active` means the family passed a specificity test and an owner sign-off (`docs/reports/2026-10-07-pfam-family-proposals.md:7`; `docs/superpowers/specs/2026-10-04-orchestrator-design.md:259-260`). Four module names exist: `pfam_adhesion`, `pfam_allergen`, `pfam_hydrophobin`, `pfam_hsba` (`modules/pfam.py:22`).

Status column: state on `main` / state on `origin/lectin-reviews-and-full-run` (that branch turns on 8 more rows; read from `git show origin/lectin-reviews-and-full-run:data/sorting_hat/family_table.tsv`).

| Family (Pfam) | Wall or fusion role | File:line | Status main / branch | Call or control it feeds |
|---|---|---|---|---|
| CFEM PF05730 | GPI wall protein family (Rbt5, Pga10, Ag2/PRA) | `data/sorting_hat/family_table.tsv:2` | active (`no_tm`) / active | `wall_family_domain` via `pfam_adhesion` (`src/cellsurface_sorting_hat/categories.yaml:50-51`), so also `cell_wall_adhesion_candidate` (`categories.yaml:56-61`) |
| PA14 PF07691 | flocculin and Epa lectin domain; also in beta-glucosidases | `family_table.tsv:16` | active (`signal_peptide`) / active | `wall_family_domain` |
| Hydrophobin PF01185, PF06766, PF28987, PF22354, PF29785, PF29802, PF29465 | surface-active wall proteins (RodA, RodB, RodD) | `family_table.tsv:4,5,6,7,17,18,19` | active / active | `hydrophobin_domain` (`categories.yaml:52-53`); in `surface_attachment_candidate` (`categories.yaml:66-72`) and in the `other_*` mechanism list (`categories.yaml:97,103`) |
| HsbA PF12296 | hydrophobic surface binding protein | `family_table.tsv:20` | active / active | `hsba_domain` (`categories.yaml:54-55`), same two uses |
| Bys1 PF04681 | yeast-phase surface protein (CalA) | `family_table.tsv:3` | inactive / inactive | none. The secretome discovery lists it as a cell wall family in 27 of 27 Onygenales (`origin/allergen-predictor-wip:docs/reports/2026-10-10-secretome-discovery-pass.md:41`) |
| PIR PF00399 | covalently bound Pir wall protein repeat (structural, not adhesin) | `family_table.tsv:15` | inactive / inactive | none. The note says "counted as a wall family only" |
| PIR1-like_C PF22799 | cell wall mannoprotein | `family_table.tsv:24` | inactive / active | `wall_family_domain` on the branch |
| Hyr1 PF15789, Hyphal_reg_CWP PF11765 | hyphally regulated GPI wall proteins (Hyr/Iff) | `family_table.tsv:23,14` | inactive / PF15789 active, PF11765 inactive | `wall_family_domain` on the branch (PF15789) |
| Flocculin PF00624, Flocculin_t3 PF13928, Flo11 PF10182, GLEYA PF10528 | flocculin and Epa-type adhesin domains | `family_table.tsv:21,22,12,13` | inactive / active [corrected 2026-10-17, see 07 section 2 tables: Flocculin, Flocculin_t3, Flo11 and GLEYA are active=yes on main now (owner 2026-10-09, commit d94a8ab); the "inactive" main state is the pre-10-09 state] | `wall_family_domain` on the branch |
| Candida_ALS_N PF11766, Candida_ALS PF05792, ALS_M PF30910 | Als adhesins | `family_table.tsv:8,9,25` | inactive / PF11766 and PF05792 active, PF30910 inactive [corrected 2026-10-17, see 07 row 6: PF11766 and PF05792 are active=yes on main since commit d94a8ab (2026-10-09); only PF30910 is still inactive; the "inactive" main entry describes the pre-10-09 state] | `wall_family_domain` on the branch |
| AltA1 PF16541, Allergen_Asp_f_4 PF25312 | allergen families, not wall enzyme families | `family_table.tsv:10,11` | inactive / inactive | none |

The wall families in the table are adhesion and structural wall protein families. **No row is an enzyme family** (GH72, GH16, GH17, GH18, chitin synthase, glucan synthase), and there is no enzyme module in `MODULES` (`modules/pfam.py:22`).

Rules that matter for remodeling enzymes:

- "GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get no finer label in version 1. There is no `cell_wall_protein` call." (`docs/superpowers/specs/2026-10-04-orchestrator-design.md:229-230`, known limit 3.)
- "Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, non-protein adhesins, moonlighting proteins and biofilm are not categories." (same file, `:246-247`.)
- Enzyme classes are deferred: D2 (`:449`), D5 (`:452`), and the sketch `wall_remodeling_enzyme`, `wall_degradation_enzyme`, `wall_synthesis_enzyme` (`:297-308`; rows at `:303-305`). The sketch names PF03198 and PF00704 for remodeling, and PF01644, PF03142, PF02364 plus GPI biosynthesis for synthesis. Synthesis needs a validated TM module.
- Large families (GH18, GH5, Asp, Cu-oxidase) need a second condition such as a signal peptide before a row can be active (`:261-262`).
- The roadmap lists "Enzyme classes (spec §3.7)" as owner-later work (`docs/ROADMAP-2026-10-08-DRAFT.md:273`, row 47; also `:79`).

### 1.2 Wall and fusion-related classes that are mapped but not implemented

`docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md` is a map, not a test: "No class below was run through any tool for this note" (`:7-8`). Owner scope decision: synthase, glucanase and degradation-enzyme classes are later work (`:10`). Rows (`:49-65`):

| Class | Pfam named by the report | Line | Repo status |
|---|---|---|---|
| Transglycosidases (GEL, PHR, GAS, Crh, Sps2, Dfg) | Glyco_hydro_72 PF03198 ("to verify"), Glyco_hydro_16 PF00722 ("to verify") [corrected 2026-10-17, see 05 Table A: both IDs and names verified in InterPro (Glyco_hydro_72, Glyco_hydro_16)] | `:53` | not in any table; no label for enzymes |
| Yapsins and wall proteases | Asp PF00026 (generic) | `:54` | none |
| Chitinases, glucanases, chitin deacetylases | Glyco_hydro_18 PF00704, Cellulase (GH5) PF00150, Polysacc_deac_1 PF01522 | `:55` | none |
| Chitin synthases | PF01644, PF03142, PF08407 (see section 0) | `:56` | none |
| beta-1,3-glucan synthase | Glucan_synthase PF02364 | `:57` | none |
| beta-1,6-glucan (Kre9) | Kre9_KNH PF10342; KRE2, KRE5, KRE6, AGS "not looked up" | `:58` | none |
| GPI anchor biosynthesis | GWT1 PF06423, Gpi1 PF05024, Gpi16 PF04113, Gaa1 PF04114, GPI2 PF06432 | `:60` | none |
| Wall integrity sensors | WSC PF01822 (also in secreted enzymes) | `:62` | none; out of scope for v1 |
| Septation, polarized growth | not looked up | `:63` | out of scope |

The same report states the design finding: GPI remodeling enzymes and hydrolases sit inside the step 1 population but are not adhesins, and "A sorter that has 'adhesin' as its only cell wall label would call these 'other'" (`:37-41`).

### 1.3 Other analyses that touch wall families

- **CBM18 / Chitin_bind_1 PF00187** in chytrids: `analysis/chytrid_batrach/NOTES.md:22` (CAZy overview has no CBM18) with the call table at `:25-30`; census row at `docs/model-review/2026-09-27-review-and-framework-plan.md:237`.
- **GH18 chitinase CTS1 as a structural control** for the CFEM fold work: `analysis/class2b_structure/candidates.tsv:11` (also `analysis/class2b_pra3_fulllength/candidates_fulllength.tsv:11`), PF00704, class `enzyme_control`; self hit at `docs/reports/2026-09-29-class2b-structure.md:178`.
- **Pfam annotations on protein rows** (labels on proteins, not families in any call): PF01522 (chitin deacetylase) at `data/curated/surface/surface.tsv:362,959`; PF01476 (LysM) at `surface.tsv:151,481`; PF00332 (on BGL2 rows, name not verified) at `surface.tsv:930,3051`; PF03663 (on mannan endo-1,6-alpha-mannosidase rows, name not verified) at `surface.tsv:172`; PF03198 with PF07983 on GEL rows at `data/curated/antigens/coccidioides_candidates.tsv:57`.
- **Secretome discovery pass** (branch only): wall families in the enriched sets are Tos1 PF10287 and PF10290, Ecm33 PF12454, Bys1; for Aspergillaceae by copy number SUN PF03856, Kre9/KNH PF10342, X8 PF07983, LysM PF01476, Ecm33 PF12454 (`origin/allergen-predictor-wip:docs/reports/2026-10-10-secretome-discovery-pass.md:41` and `:62`). Accessions and names are as quoted in that report. They are clade associations and not tool classes.
- **Lectin-like families, report only** (branch only): LysM PF01476, CBM_14 PF01607 (Avr4) and Chitin_bind_1 PF00187 are reported-only models (`origin/lectin-reviews-and-full-run:analysis/lectin_class_report/lectin_report.py:69-84`, `MODELS`). The lectin spec records them as H1b, binding the fungus's own chitin to shield it from host or plant chitinases (`origin/lectin-reviews-and-full-run:docs/superpowers/specs/2026-10-10-lectin-class-design.md`, section 3.1 table). This is the nearest existing coverage of wall chitin binding.

### 1.4 Cell fusion

**Nothing in the repo targets cell fusion.** I searched for Fus1, Prm1, Fig1, Fus2, Kar5, Hap2/Gcs1, "cell fusion" and "mating" across `docs`, `src`, `analysis/curation`, `data/curated/adhesins`, `data/sorting_hat` and `TODO.md`. The only hits are mating adhesins:

- AGA1 and SAG1 are E1 adhesin rows (`data/curated/adhesins/adhesins.tsv:38,45`). AGA1, AGA2, SAG1 and FIG2 are among the 9 "known adhesins" of the S288C check (`docs/model-review/2026-09-27-review-and-framework-plan.md:152`).
- SAG1 hits Candida_ALS_N PF11766 (`docs/reports/2026-10-07-pfam-family-proposals.md:41`). [corrected 2026-10-17, see 07 rows 4, 5, 8: confirmed, score 244.9; PF05792 (Candida_ALS) does not hit SAG1 (none at GA, none at E 10)]
- A CGD gene headline says a Candida protein is "repressed during the mating process" (`docs/reports/data/sorting_hat/specificity/Calb_SC5314/cgd_gene_names.json:12`). It is a gene annotation and not a class.

Fusion membrane proteins (Fus1, Prm1, Fig1, Hap2) appear in no table.

## 2. Wall housekeeping proteins as hard negatives and false-positive sources

| Finding (quoted or tightly paraphrased) | Where |
|---|---|
| "On S. cerevisiae S288C the shipped 8M model makes 77 calls ... precision is ≈12%: the other calls are GPI cell-wall mannoproteins (SED1, TIR4, DAN1/4, CCW12...), mucin-like sensors (MSB2, HKR1, WSC2–4) and dubious ORFs. It also calls 12/20 curated hard negatives." | `docs/model-review/2026-09-27-review-and-framework-plan.md:22-23` |
| "The top CAZy family among the calls is AA1 multicopper oxidases/laccases (12.5k calls), which are not adhesins." | `docs/model-review/2026-09-27-review-and-framework-plan.md:20` |
| False-positive list includes CRH1, UTR2, SUN4, SCW11, KRE1/9, NCW2, CTS1; sensors MSB2, HKR1, WSC2/3/4, MID2, SLG1, MTL1. | `docs/model-review/2026-09-27-review-and-framework-plan.md:163-164` |
| "Every fungal GPI cell-wall protein, mucin-like sensor (Msb2, Hkr1), secreted glycosidase with an S/T-rich linker, and heavily O-glycosylated extracellular enzyme lands on the positive side." | `docs/model-review/2026-09-27-review-and-framework-plan.md:371-374` |
| Hard negatives to curate: Gas1, Crh1, Cwp1/2, Ccw12, Ccw22, Sed1, Tir1/3/4, Dan1, Ecm33, Yps1, Ncw2, the Pir family; glucanases (Exg1), chitinases, laccases. | `docs/model-review/2026-09-27-review-and-framework-plan.md:440-443` |
| Sorting hat on S288C: "R0 removes the intracellular proteins and keeps nearly every secreted wall protein ... Hard negatives: 27 of 28 called in S288C. R0 is a filter for surface candidates. It does not separate adhesins from other wall proteins." | `docs/reports/2026-10-07-classifier-runs-and-scale.md:56` |
| "The mechanism calls do not separate them either": the repeat call fires for 5 of 27 S288C hard negatives (PIR1, HSP150, TIR1, TIR4, FIT1) and for 0 of 3 strong-evidence adhesins. | `docs/reports/2026-10-07-classifier-runs-and-scale.md:57` |
| Gel1 (beta-1,3-glucanosyltransferase, Coccidioides) scored 0.000: "missed (an enzyme — arguably correct)". | `docs/model-review/STATUS.md:193` |
| In the named panel, MSB2, HKR1 and SAP9 are labelled negative and most candidates call them positive; "This was not analysed." | `docs/model-review/STATUS.md:87`; `docs/paper/04-limits-and-open-questions.md:43` |
| The antigen ranking top was pan-fungal wall enzymes: "Gel1 glucanosyltransferases, pepsins, subtilisins, chitinases". | `docs/reports/2026-09-27-coccidioides-antigen-findings.md:53` |
| A CFEM hit sets `cell_wall_adhesion_candidate` (Ag2/PRA in RS); 18 of the 22 domain-route hits are CFEM; CFEM is a wall family not shown to mediate adhesion. | `origin/lectin-reviews-and-full-run:docs/superpowers/specs/2026-10-09-v1-class-list-design.md:13-14` |

**Curated wall housekeeping rows** (hard negatives, evidence level N1 or N2), `data/curated/adhesins/hard_negative_seeds.tsv`:

- GH72 GAS1 (`:4`); GH16 CRH1 (`:5`) and UTR2 (`:29`); GH5 EXG1 (`:16`); GH18 CTS1 (`:28`); SUN4 (`:27`); KRE9 (`:30`); yapsin YPS1 (`:15`); ECM33 (`:14`).
- Wsc and Mid2 sensors: WSC1, WSC2 (`:19,20`), MID2 (`:21`). Pir: PIR1, CIS3, HSP150 (`:22-24`). Mannoproteins CWP1, CWP2, TIR1, TIR3, TIR4, DAN1, DAN4 (`:6,7,9-13`). Laccase ABR1 (`:32`).
- The seed list has 31 rows and "Hard negatives are almost all S. cerevisiae" (`data/curated/adhesins/README.md:36`). Issue #14 is open (`docs/superpowers/specs/2026-10-04-orchestrator-design.md:370`).
- `data/curated/adhesins/adhesins.tsv` carries the same N1 rows with accessions (gel1 `:156`, EXG1 `:164`, GAS1 `:167`, KRE9 `:170`, MID2 `:171`, SLG1 `:174`, WSC2 `:177`), plus 36 rows of class `surface_other_adhesion_phenotype`: surface enzymes whose mutant affects adhesion, "likely indirect" (`data/curated/adhesins/README.md:20`), for example BGL2 `:114`, PHR1 `:125`, SUN41 `:134`, XOG1 `:136`. These are the best existing seeds for a remodeling gene set.
- `data/curated/surface/surface.tsv` labels 1,257 proteins `non_adhesin_putative` ("proposed: annotated catalytic activity, no adhesion evidence", `data/curated/README.md:29`). They include Aspergillus chitinases, chitin deacetylases and glucanases (for example `surface.tsv:51,68,151,346,362`). The label is derived from annotation and is not reviewed.
- Task 07 defines shared hard-negative classes `gpi_wall_enzyme`, `wall_hydrolase`, `st_linker_enzyme`, `laccase_etc`, `domain_non_member` (`docs/agent-tasks/07-hard-negatives-shared.md:24-30`). It is a plan, not a table.

**Limit of this evidence.** The numbers above come from the earlier ESM classifier (`docs/model-review/`) and from R0 plus the active mechanism calls (`docs/reports/2026-10-07-classifier-runs-and-scale.md`). No measurement shows how a proposed enzyme or fusion family would behave in the sorting hat, because none is implemented.

## 3. Gaps (candidates only)

"In repo" means the ID appears in text I read.

### 3.1 Cell wall remodeling and synthesis

| Candidate family or gene set | Pfam ID | Verification | Note |
|---|---|---|---|
| GH72 (GEL, GAS, PHR; Gel1, Gas1, Phr1) | PF03198 | In repo (`2026-10-04 report:53`, itself "to verify"; `orchestrator spec:303`; annotation on `coccidioides_candidates.tsv:57`) | Often with the X8 domain PF07983 (same row; name X8 per the discovery report) [corrected 2026-10-17, see 05 Table A (PF03198 and PF07983 names verified) and 07 row 17 (Gas1 hits PF03198 504.4 and X8 102.4)]. GPI-anchored. Hard negatives exist. |
| GH16 (Crh, Utr2) | PF00722 | In repo (`2026-10-04 report:53`, "to verify") [corrected 2026-10-17, see 05 Table A: ID and name verified; 07 row 18: Crh1 and Utr2 hit it, 170.8 and 159.5] | Also holds non-wall glucanases; needs a signal peptide and GPI second condition. |
| GH76 (Dfg5, Dcw1 alpha-1,6-mannanases) | PF03663 | ID in repo only as an annotation on mannan endo-1,6-alpha-mannosidase rows (`surface.tsv:172`); name and the Dfg5/Dcw1 link are my recall [corrected 2026-10-17, see 05 Table A: PF03663 = Glyco_hydro_76 verified (Pfam type is repeat); 07 row 19: Dfg5 530.8 and Dcw1 547.4 hit it by hmmsearch --cut_ga] | verified ID and name (05); Dfg5/Dcw1 hits measured (07); CAZy records alone do not name Dfg5/Dcw1 (05 Table B) |
| GH17 (Bgl2, Exg-like) | PF00332 | ID only as annotation on BGL2 rows (`surface.tsv:930,3051`); name unverified [corrected 2026-10-17, see 05 Table A: PF00332 = Glyco_hydro_17 verified; 07 row 20: Bgl2 hits it at 353.2, Scw4 45.6] | verified ID and name (05); Bgl2 hit measured (07) |
| GH5 glucanases (Exg1) | PF00150 | In repo (`2026-10-04 report:55`) | large family with intracellular members |
| GH18 chitinases (Cts1, Cts2) | PF00704 | In repo (`orchestrator spec:303`; `2026-10-04 report:55`) | very large; needs a second condition |
| Chitin deacetylase (CE4) | PF01522 | In repo (`2026-10-04 report:55`) | Cda2 of Cryptococcus is already a hard-negative candidate (`data/curated/biofilm/basidiomycete_seeds.tsv:63`) |
| Chitin synthase | PF01644, PF03142, PF08407 | In repo; names conflict with the brief (section 0) | Membrane protein. Class I to VII needs phylogeny (`orchestrator spec:263-264`). Chs3 regulators Chs4, Skt5: unverified: check InterPro/Pfam before use. |
| beta-1,3-glucan synthase (Fks, Gsc) | PF02364 | In repo | membrane protein; other Fks domains: unverified: check InterPro/Pfam before use |
| alpha-1,3-glucan synthase (Ags) | none | not in repo ("not looked up", `2026-10-04 report:58`) | unverified: check InterPro/Pfam before use |
| Kre6/Skn1 (beta-1,6-glucan) | none | not in repo | unverified: check InterPro/Pfam before use |
| Kre9/Knh1 | PF10342 | In repo (`2026-10-04 report:58`) | The same domain is on uncharacterized Aspergillus GPI proteins (`surface.tsv:216,327`), so a hit is not Kre9 by itself. |
| SUN family (Sun4, Sun41, Sim1, Uth1, Nca3) | PF03856 | In repo only as quoted in the discovery report | SUN41 is a curated `surface_other_adhesion_phenotype` row (`adhesins.tsv:134`) |
| Ecm33 / Sps2 | PF12454 | In repo only as quoted in the discovery report [corrected 2026-10-17, see 05 Table A (PF12454 = Ecm33 verified) and 07 row 22: S. cerevisiae Ecm33 has NO hit to PF12454, none even at E 10, and no Pfam-A hit at GA, so "ECM33 a member" is refuted for the S288C protein] | ECM33 is an N1 negative (`hard_negative_seeds.tsv:14`) |
| Tos1 | PF10287, PF10290 | In repo only as quoted in the discovery report | |
| Scw4/Scw10, Ssg1 and other soluble wall glucanases | none | not in repo | unverified: check InterPro/Pfam before use |
| Yapsin (Yps), Sap9/Sap10 | PF00026 | In repo, generic | needs curated examples (`2026-10-04 report:54`) |
| Wall integrity sensors (Wsc) | PF01822 | In repo | WSC also in secreted enzymes (`2026-10-04 report:62`). Mid2, Mtl1 Pfam: unverified: check InterPro/Pfam before use |
| Expansin-like (swollenin), cerato-platanin | none | Expansin-like named as an unverified look-alike (`docs/superpowers/specs/2026-10-08-hydrophobin-custom-hmm-design.md:147`). Cerato-platanin is out of scope by owner statement (`origin/lectin-reviews-and-full-run:docs/superpowers/specs/2026-10-09-v1-class-list-design.md:26`). | unverified: check InterPro/Pfam before use |
| GPI anchor biosynthesis (Gpi1, Gpi2, Gwt1, Gaa1, Gpi16, Pig-X) | PF06423, PF05024, PF04113, PF04114, PF06432 | In repo (`2026-10-04 report:60`). PF10333 appears as a Pig-X annotation (`surface.tsv:1534`) | ER membrane; no GPI module; "GPI as a class" is not in v1 (`v1-class-list-design.md:118` on the branch) |
| Chitin and glucan binding modules LysM, CBM_14, Chitin_bind_1 | PF01476, PF01607, PF00187 | In repo (`lectin_report.py:69-84` on the branch) | reported-only in the lectin set; no remodeling context |
| Pir/Cwp structural proteins | PF00399, PF22799 | In repo (table rows 15, 24) | inactive on main |
| Mannosyltransferases, galactofuran enzymes (Och1, Mnn9, Ugm1, GfsA) | none | not in repo (`2026-10-04 report:59`, "not looked up") | unverified: check InterPro/Pfam before use |
| Rho1, Pkc1, Bck1, MAPK wall integrity pathway | none | cytosolic, out of scope (`2026-10-04 report:33,62`) | |

### 3.2 Cell fusion

No family in the repo. Every ID is unverified. I give none. [corrected 2026-10-17, see 07 row 26: Pfam models exist for Prm1 and Fig1 (below); they are not in the family table]

| Candidate | Role | Pfam ID |
|---|---|---|
| Fus1 | fusion zone wall protein | unverified: check InterPro/Pfam before use |
| Prm1, Fig1, Kar5, Prm2 | membrane and cell fusion | unverified: check InterPro/Pfam before use [corrected 2026-10-17, see 07 rows 26, 27 and Fusion table: Prm1 = PF27758 (hit 475.5 on S. cerevisiae PRM1) and Fig1 = PF12351 (hit 192.3 on FIG1) verified by hmmsearch; neither is in the family table; Kar5 and Prm2 remain unverified] |
| Hap2/Gcs1 | gamete fusion fusogen | unverified: check InterPro/Pfam before use |
| Ham-2, Ham-5, Ham-6, Ham-7 (Neurospora anastomosis) | cell-cell fusion | unverified: check InterPro/Pfam before use |
| Mating agglutinins Aga1, Aga2, Sag1 | mating adhesion before fusion | Sag1 hits PF11766 in repo (above); Aga1 and Aga2: unverified: check InterPro/Pfam before use [corrected 2026-10-17, see 07 row 7: Aga1 has no Pfam-A hit; Aga2 hits PF17366 (AGA2, 95.8), which is not in the family table; neither hits a table row]. FIG2 is a curated adhesin row with no family [corrected 2026-10-17, see 07 row 9: FIG2 hits Flocculin_t3 PF13928 (64.8, five weak domains of 24-28 each), an active row; the "no family" statement is refuted]. |
| Fus2, Rvs161/Rvs167, Spa2 | fusion zone scaffolds | cytosolic, out of scope for a surface tool |

### 3.3 Structural gaps in the tool

1. No enzyme module in `MODULES` (`modules/pfam.py:22`). `load_family_table` rejects any other module name (`:60-61`).
2. No call or report line for remodeling enzymes or fusion proteins. v1 has five categories (`orchestrator spec:21-22`, D2 at `:449`).
3. `SECOND_CONDITIONS` has only `""`, `signal_peptide`, `no_tm` (`modules/pfam.py:23`). There is no GPI condition and no GPI module.
4. Membrane proteins (synthases, fusion proteins) fall outside the signal peptide gate (`categories.yaml:44-46`). The spec says synthases would need a TM segment and a validated TM module (`orchestrator spec:305`). The `tm` module is listed as an evidence field only (`origin/lectin-reviews-and-full-run:docs/CLASSES.md:35`).
5. The shared hard-negative table (issue #14) is not built. Without it no remodeling family can pass the specificity test required for `active` (`orchestrator spec:259-262`).

## 4. Suggestion: report-only context rows for "cell wall remodeling" and "cell fusion"

Description only. Nothing here is implemented.

### 4.1 The precedent (found on branches, not on main)

The owner's decision of 2026-10-10: D1 option (c), no lectin call, "All lectin-like models are reported only"; D12 option **(c-i)**, "a script outside the tool. No change to `src/`, `family_table.tsv` or `categories.yaml`" (`origin/lectin-reviews-and-full-run:docs/superpowers/specs/2026-10-10-lectin-class-design.md:5-16`). Section 3.4 of that spec explains why the in-tool option (c-ii) is worse: rows must be `active=yes` to produce output, so `active` stops meaning "passed a specificity test", and an `evidence:` edit changes `config_sha256` and makes existing call status files stale.

The M2 spec states the same staleness rule: "any edit of `categories.yaml` makes the repeat call files stale" (`origin/lectin-reviews-and-full-run:docs/superpowers/specs/2026-10-09-v1-class-list-design.md:120`; also `:59`).

### 4.2 Recommended pattern: copy the lectin script, one script per gene set

Copy `analysis/lectin_class_report/` (branch `origin/lectin-reviews-and-full-run`):

- `lectin_report.py`: runs `hmmsearch --cut_ga` per proteome against a fixed `MODELS` list, writes counts, hits and a markdown report. It makes no call, writes no status file and does not touch `src/`. Optional `--signalp` and `--name-list` inputs (header text, `lectin_report.py:1-35`).
- `submit_lectin_report.sbatch`, `compare_committed.py`, `README.md` (versions, run commands).
- Tests: `tests/lectin_class_report/test_lectin_report.py` (21 test functions; names include `test_counts_hits_and_label`, `test_report_text`, `test_a2_group_only_with_signalp`, `test_domain_below_ga_is_flagged`, `test_zst_input`, `test_needs_tmpdir_or_scratch`, `test_missing_model_stops`, `test_duplicate_names_stop`, `test_duplicate_id_inside_one_proteome`, `test_empty_fasta_does_not_abort`, `test_name_list_bare_accession_and_crlf`).

Proposed new files (names are suggestions):

| New file | Content |
|---|---|
| `analysis/cellwall_context_report/cellwall_report.py` | Same CLI. `MODELS` for the remodeling set, limited at first to IDs marked "In repo" in section 3.1. Per hit: protein ID, accession, name, score, E-value, below-GA flag, optional R0 result and optional TM count. |
| second `MODELS` block or `fusion_report.py` | Fusion set. Cannot be written until the IDs of section 3.2 are verified. |
| `analysis/cellwall_context_report/README.md`, `submit_*.sbatch` | Versions and commands, as the lectin README. |
| `tests/cellwall_context_report/test_cellwall_report.py` | See 4.4. |
| A model table (label, role, paper, whole-proteome genus share) | The role the lectin spec's section 3.1 table plays. |

**No change** to `src/cellsurface_sorting_hat/`, `data/sorting_hat/family_table.tsv`, `categories.yaml`, `docs/class_descriptions.yaml` or `docs/CLASSES.md` (on the branch a test keeps `docs/CLASSES.md` and the call names in sync, `v1-class-list-design.md:59-60`). Every existing status file stays valid.

### 4.3 What each report row should say

1. A family-only label with fixed text, in the style of `PF07938_LABEL` and `FIXED_TEXT` (`lectin_report.py:53-60`): "This is a hit to the Pfam model X. It does not show that the protein remodels the wall, is on the cell surface, or is expressed. No call is made."
2. Optional columns from outputs the tool already makes: the R0 result and `tm.n_tm_mature` (listed as evidence fields at `categories.yaml:40-41`). This separates GPI and secreted enzymes from membrane synthases by display only.
3. A note drawn from section 2 that housekeeping wall enzymes are not adhesins, and no use of the word "adhesin" as a claim (M2 rule: "Evidence is not a claim", `v1-class-list-design.md:61`).
4. A genus count per family over the Fungi_5k scan, as the lectin report and the secretome discovery do.
5. For the fusion set: a note that fusion proteins are membrane proteins outside the step 1 population.

### 4.4 Tests needed

- Counting: a stub `hmmsearch` or toy domain table gives the expected per-model counts and one row per protein and model (pattern of `test_counts_hits_and_label`).
- Report text: the fixed label sits beside each hit and no row claims adhesion (pattern of `test_report_text`).
- `--cut_ga` enforced: a domain below GA is flagged (`test_domain_below_ga_is_flagged`; `modules/pfam.py:94-95` already refuses a table without `--cut_ga`).
- Failure paths: missing model, duplicate names, empty FASTA, bad SignalP file, duplicate IDs, `.zst` input, no `--tmpdir` and no `$SCRATCH` (the existing lectin tests).
- A `MODELS` drift test: each accession carries the name recorded in the script, checked against a small committed expectation, so an unverified ID cannot enter silently. This test is new and closes the verification gap of section 3.
- A guard test that the script writes no status file and changes no file under `src/`, `data/sorting_hat/` or `docs/CLASSES.md`.

### 4.5 If the owner later wants these rows inside the tool (option c-ii)

This is the route the lectin spec decided against. Files that would change, from code I read on `main`:

| File | Change | Consequence |
|---|---|---|
| `src/cellsurface_sorting_hat/modules/pfam.py:22` | add module names to `MODULES`; `:60-61` otherwise rejects the rows | module-set tests must change |
| `src/cellsurface_sorting_hat/modules/cli.py:71` | help text lists the module names | |
| `data/sorting_hat/family_table.tsv` | rows with `active=yes`, `active_by`, `active_date` (`pfam.py:68-69` requires them), and a `specificity_note` of "report only, no specificity test" | `active` no longer means a passed specificity test |
| `src/cellsurface_sorting_hat/categories.yaml` | only if fields should reach `evidence.tsv.gz` | changes the config hash; existing call status files go stale and the repeat files must be re-measured |
| `tests/cellsurface_sorting_hat/test_data_and_scripts.py` | asserts row count, the signed-off set and the module set (`:30-60` on the branch) | edit all three |
| `tests/cellsurface_sorting_hat/test_engine.py`, `test_cli.py`, `modules/test_module_cli.py`, `golden/*` | module lists in fixtures and the golden report (for `pfam_hsba`: `test_engine.py:52,309,369,383`; `test_cli.py:44`; `golden/report.expected.md:22,39,57,58`) | regenerate golden files and read every changed line |

**Recommendation.** Use (c-i) for both sets. The enzyme families here are large with intracellular members (`orchestrator spec:261-262`), no hard-negative table exists (issue #14), and `active` already promises a passed specificity test. A script table gives the owner counts of remodeling and fusion families per proteome without a call.

## 5. Unknowns

- I ran no tool, test, script or `hmmfetch`. All Pfam names are from repo text.
- `PF00722` and `PF03198` are "to verify" in the 2026-10-04 report itself (`:53`). [corrected 2026-10-17, see 05 Table A: both IDs verified in InterPro]
- The GH76 = PF03663 and GH17 = PF00332 readings, and the Dfg5/Dcw1 link, are my recall. Treat as unverified. [corrected 2026-10-17, see 05 Table A and Summary item 3: the IDs and names are verified (Glyco_hydro_76 PF03663, Glyco_hydro_17 PF00332); see 07 rows 19, 20: Dfg5/Dcw1 and Bgl2 hit them by hmmsearch. 05 Table B notes CAZy records alone do not name Dfg5/Dcw1, Bgl2, Crh/Utr2 or Exg1]
- For branch files I read only the parts cited. The lectin spec has five revisions and an owner block that overrides earlier text (`:5-7`); I relied on that block and the sections cited.
- I did not read the sibling reports `01-cell-wall-reviews.md`, `02-cell-fusion.md` or `03-multicellularity-hyphae.md`.

## Errata log

Applied 2026-10-17 from reports 05, 06, 07. Line numbers are those of the file before this errata (the added header line shifts later lines by 2).

1. Line 12. Old: `I did not run `hmmfetch` or open InterPro.`. New: same text with an added marker. Source: 05 Method, Summary.
2. Line 34. Old: `| `family_table.tsv:21,22,12,13` | inactive / active |`. New: same text with an added marker. Source: 07 s2 (Active column).
3. Line 35. Old: `| inactive / PF11766 and PF05792 active, PF30910 inactive |`. New: same text with an added marker. Source: 07 s3 row 6.
4. Line 54. Old: `Glyco_hydro_72 PF03198 ("to verify"), Glyco_hydro_16 PF00722 ("to verify")`. New: same text with an added marker. Source: 05 Table A.
5. Line 79. Old: `- SAG1 hits Candida_ALS_N PF11766 (`docs/reports/2026-10-07-pfam-family-proposals.md:41`).`. New: same text with an added marker. Source: 07 s3 rows 4, 5, 8.
6. Line 119. Old: `Often with the X8 domain PF07983 (same row; name X8 per the discovery report).`. New: same text with an added marker. Source: 05 Table A; 07 row 17.
7. Line 120. Old: `| In repo (`2026-10-04 report:53`, "to verify") |`. New: same text with an added marker. Source: 05 Table A; 07 row 18.
8. Line 121. Old: `name and the Dfg5/Dcw1 link are my recall | unverified: check InterPro/Pfam before use |`. New: same text with an added marker. Source: 05 Table A/B; 07 row 19.
9. Line 122. Old: `name unverified | unverified: check InterPro/Pfam before use |`. New: same text with an added marker. Source: 05 Table A; 07 row 20.
10. Line 132. Old: `| In repo only as quoted in the discovery report | ECM33 is an N1 negative`. New: same text with an added marker. Source: 05 Table A; 07 row 22.
11. Line 146. Old: `No family in the repo. Every ID is unverified. I give none.`. New: same text with an added marker. Source: 07 row 26.
12. Line 151. Old: `| Prm1, Fig1, Kar5, Prm2 | membrane and cell fusion | unverified: check InterPro/Pfam before use |`. New: same text with an added marker. Source: 07 s2, rows 26-27.
13. Line 154. Old: `Aga1 and Aga2: unverified: check InterPro/Pfam before use. FIG2 is a curated adhesin row with no family. |`. New: same text with an added marker. Source: 07 rows 7, 9.
14. Line 230. Old: `are "to verify" in the 2026-10-04 report itself (`:53`).`. New: same text with an added marker. Source: 05 Table A.
15. Line 231. Old: `are my recall. Treat as unverified.`. New: same text with an added marker. Source: 05 Table A, Summary 3; 07 rows 19, 20.

Still unverified after 05 to 07: PF10333 as Pig-X; CE4 for PF01522; Pfam-to-CBM links; Avr4; Kar5/Prm2/Fus1 IDs; 12/20 hard-negative figure is not testable on the 07 set (07 row 15).
