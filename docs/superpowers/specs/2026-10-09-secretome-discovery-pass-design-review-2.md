# Review 2 of the secretome discovery pass spec (rev 2)

*2026-10-09. Independent second-round review of `docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md`, revision 2. First review: `...-design-review-1.md`. Reviewer: Claude Code (claude-opus-5-5). The spec was not edited. Numbers below come from the files and commands named in each finding. Genus counts were recomputed with the three-part rule of `analysis/fungi_scan/select_proteomes.py` (FASTA in `input/`, at least 4,000 proteins, taxid in `_workdir/sorting_hat/taxdump/nodes.dmp`). Tree checks used Biopython `Phylo` on `Phylogeny/fungi_tree/final_tree.nw`, rerooted on a Basidiomycota tip. Fisher p-values were computed with exact hypergeometric sums (`math.comb`).*

## Verdict

**PASS.** No blocker. Four major findings and eight minor findings.

The blocker of review 1 (B-1) is fixed. Most first-round findings are fixed in substance. Three of the fixes added in rev 2 have a defect of their own: the proteome choice rule (m-6 fix), the origin count (M-2 fix) and the Pfam-free cluster source (M-4 fix). One first-round finding (M-3, secretome size) is only partly addressed. Each major finding has a cheap fix that can go into the spec or the plan before B1. The owner can accept the spec with these fixes as conditions; they do not require a new design.

## Status of the first-round findings

| Finding | Status | Note |
|---|---|---|
| B-1 active family set | Fixed | Section 4: re-run converters and core command; B2 check of `run.json` and `pfam_*.json`. Checked: the scan's `pfam_hmmsearch.sbatch` fetched every row of `family_table.tsv` (24 rows at commit `3ff1d5e`, same as today), so the raw tables hold all 24 models and no new `hmmsearch` is needed. |
| M-1 empty genus, *Paracoccidioides* | Fixed | Counts correct (see "checked"). Recovery path works (see m-B below). |
| M-2 non-independence of genera | Fixed in wording; method defect | Label next to p and q is in. The origin count as specified does not count origins. See MAJOR M2-2. |
| M-3 copy number and secretome size | Partly fixed | Copy-number test added. The presence test has no pool-size adjustment. See MAJOR M2-3. |
| M-4 pools and Pfam-free proteins | Fixed in design; source defect | Two pools, TM and HDEL/KDEL exclusion, Pfam-free units. The named no-Pfam dump excludes small clusters. See MAJOR M2-4. |
| M-5 clade enrichment vs host interaction | Fixed | Section 1 sentence; pathogen contrast in 5.3; limit 6.2. |
| M-6 compute cost | Fixed | Job count is in section 4. Arithmetic checks (see "checked"); one count is off by 7 (m-A). |
| M-7 Fungi_5k annotations | Fixed | Reuse decisions and B2a gate stated. |
| M-8 no section column | Fixed in design; tree caveat | Tree clades replace sections. *Aspergillus* is not monophyletic on the tree (MAJOR M2-2, point c). |
| M-9 post-hoc role table | Fixed | Default "no stated role", per-row basis, T4 label, blind order, first-page sentence. |
| m-1 medians | Fixed | 412.5 of 427 and 784 of 808.5 recomputed and correct. |
| m-2 step names | Fixed | "Seven module steps" is gone. |
| m-3 BH filter | Fixed | BH over families with at least one present genus; the 3-genus rule limits only reporting. |
| m-4 sidedness, sister comparison | Fixed | Two-sided; sister comparison descriptive. See m-E on depletion reporting. |
| m-5 class allocation | Mostly fixed | Rule stated. Empty-class genera still not handled; the cap never binds (m-D). |
| m-6 proteome quality | Fixed in wording; rule defect | See MAJOR M2-1. |
| m-7 antigen module | Fixed | "*C. immitis* RS only". |
| m-8 HPCC rules | Fixed | `${SCRATCH:?}`, absolute paths, no `BASH_SOURCE[0]`, `.tsv.gz`, resolved Pfam 38.2 path. Full Pfam-A is now reused, so the chunking point no longer applies. |
| m-9 fixed levels | Fixed | Four levels fixed in 5.3; effect size and genus counts required. |

## Findings

### BLOCKER

None.

### MAJOR

**M2-1. The fallback rule "highest protein count" selects hybrids, possible mislabels and inflated annotations. BUSCO completeness is not available in Fungi_5k, so the fallback is the rule.**

- Evidence (BUSCO): `samples.csv` has `BUSCO_LINEAGE` only, not completeness. `scripts/prep_for_bigquery_load.py` has a `busco()` loader that reads `results/stats/busco`; that directory does not exist and the call is commented out (`#busco()`). `sql/schema.sql` has no BUSCO score column. A search of Fungi_5k to depth 4 (excluding `work/`) found no BUSCO or compleasm output. So the spec's "BUSCO completeness if the file in Fungi_5k gives it" always falls through to protein count.
- Evidence (what the maximum picks; passing proteomes, ratio of maximum to genus median protein count): *Epichloe* 19,583 vs median 7,396 (2.65x; pick `Epichloe_coenophiala_ATCC_90664`, a known interspecific hybrid); *Coniochaeta* 2.05x; *Aureobasidium* 1.94x; *Exophiala* 1.85x; *Verticillium* 1.72x (pick `Verticillium_longisporum_PD589`, a known allodiploid hybrid); *Penicillium* 17,770 vs 10,110 (1.76x; pick `Penicillium_cf._chrysogenum_S_N-308-OC-P1`); *Aspergillus* 14,965 vs 11,256 (pick `Aspergillus_latus_MM151978`, a known allodiploid hybrid); *Ascosphaera* 1.56x; *Blastomyces* 1.37x. On the rerooted tree, `Penicillium_cf._chrysogenum_S_N-308-OC-P1` sits inside the main *Aspergillus* clade, not with *Penicillium*.
- Why it matters: hybrids roughly double copy number, which feeds straight into the copy-number test (5.2). A high protein count is also a sign of split gene models, which lose N-termini and signal peptides; that is the problem m-6 raised. In the background, 149 genera have 3 or more passing proteomes, so the rule affects many genera.
- Fix: choose the proteome nearest the genus median protein count (or a RefSeq `GCF_` assembly where one exists), and exclude known hybrids by name. Optional: run compleasm or BUSCO on the candidate proteomes (protein mode is fast) and choose by completeness. In B1, check that each chosen proteome sits inside its labelled genus or family on the tree, and record the ones that do not. Change the sentence in section 3 to say that no BUSCO completeness is available in Fungi_5k.

**M2-2. The "number of independent origins by Fitch parsimony" is not an origin count. The tree is rooted arbitrarily, and the target clades are not clean on it.**

- (a) Fitch parsimony gives the minimum number of state changes. It does not separate gains from losses. A family that is ancestral to Pezizomycotina and lost many times gets a high Fitch count, and the spec would print it as many "independent origins". That reads as convergent gain, which is the opposite of the truth, and makes a lineage marker look like independent support.
- (b) `final_tree.nw` is rooted on one tip: the root splits `Ustilaginoidea_albicans_UA25` (Hypocreales) from the other 5,811 tips (last characters of the file: `...,Ustilaginoidea_albicans_UA25:0.00084):0.00000;`). Gains and losses cannot be read without an outgroup root.
- (c) After rerooting on a Basidiomycota tip: the 71 Onygenales tips are not one clade (the best clade has 66 Onygenales tips; all 5 *Ascosphaera* tips are outside it). *Aspergillus* has 2 outlier tips (`Aspergillus_subversicolor_LH_A137`, `Aspergillus_foetidus_LH_A311`) and its best clade contains 7 other tips, including `Alternaria_destruens_LH_A593`, `Alternaria_alstroemeriae_LH_A124` and `Bertia_moriformis_D4H25SX9FM`. `Mycosphaerella_populi_pn0402` sits among *Penicillium*. The genus labels and the tree disagree for a small number of proteomes, and the max-count rule picked one of these (M2-1).
- Fix: root the pruned genus tree on an outgroup (for example the Pezizomycetes/Orbiliomycetes genera, or a non-Pezizomycotina tip) and state the root. Report gains and losses separately: Dollo parsimony (one gain, count losses), or Fitch with ACCTRAN and DELTRAN reconstructions and the root state. Call the column "parsimony gains (losses)", not "origins". Count gains within and into the target clade separately from gains elsewhere. Say how the Aspergillaceae and Onygenales clades are defined on the tree when the labels disagree (labels, or tree clade), and say how *Ascosphaera* is placed. Apply the same rule to the *Aspergillus* subtree used for the supplement in 5.3.

**M2-3. The presence test still has no adjustment for secretome size. With a two-sided test, Onygenales results will be dominated by depletion from a smaller pool (M-3 partly unaddressed).**

- Evidence (`proteome_summary.tsv`, `signal_peptide_protein[R0]`): Onygenales median 427 (range 192 to 738); *Aspergillus* median 808.5; the three Pezizomycotina background-like scan proteomes are 628 (*Neurospora*), 1,085 (*Fusarium*) and 1,324 (*Pyricularia*). An Onygenales pool is about half the background size. A family with a low per-protein rate is then absent in more Onygenales genera by chance. A two-sided Fisher test will report this as depletion. Review 1 asked for log pool size as a covariate or a stratified Fisher test; rev 2 adds the copy-number test (which normalises by pool size) but leaves the presence test unadjusted.
- Fix: one of: (a) logistic regression of presence on group plus log pool size, with the Fisher test kept as the ranking statistic; (b) a stratified (Cochran-Mantel-Haenszel) test over pool-size bins; (c) rarefy each pool to a common size (for example the 10th percentile pool size, with a fixed seed) before presence is scored. Report pool size per genus beside presence in every table.

**M2-4. The named no-Pfam cluster dump excludes the small clade-restricted clusters that the Pfam-free unit is meant to find.**

- Evidence: `results/mmseqs_nopfam_genedump.csv` is written by `scripts/Rscripts/mmseq_cluster_profile.R`. That script keeps only clusters with `group_count > 100` and `null_ratio > 0.95`. Recomputed from the file: 12,743 clusters, minimum `group_count` 101, 8,678 with ratio 1. A family restricted to Onygenales (71 proteomes) with one copy per proteome has at most 71 members and cannot be in the dump. SOWgp-like proteins are the stated example (5.2). The spec says the dump holds "clusters with no Pfam hit", which is true only for clusters above 100 members.
- Other formats: `Fungi5K_cluster.tsv.gz` is a two-column representative/member table of protein IDs (for example `FBBC54B5_004720-T1  FBBC54B5_004723-T1`), which joins to the chosen proteomes by protein ID. `orthogroups_mmseqs.tsv.gz` is a wide table, one column per proteome named like `Aaosphaeria_arxii_CBS 175.79` (a space before the strain, not the input file stem), with cluster IDs `MMCLUST########`. The cluster IDs reach `MMCLUST04601075`, so millions of clusters exist.
- Fix: in the spec, name `Fungi5K_cluster.tsv.gz` (or `orthogroups_mmseqs.tsv.gz`) as the source, with no size filter. Define "no Pfam hit" for a cluster from the pass's own join (for example, no member of the chosen proteomes has a Pfam-A hit). State that the Pfam-free clusters get their own BH family, and state the number tested from the data.

### MINOR

**m-A. Background count.** Section 3 says the background includes the other Eurotiales genera and cites "the review counts 603". The 603 of review 1 excludes all Eurotiales. Recomputed: 603 passing named genera outside Onygenales and Eurotiales, plus 7 other Eurotiales genera = 610. The job estimate changes by about 0.6 h. Fix: write 610.

**m-B. *Paracoccidioides* recovery.** `merged.dmp` is not extracted in `_workdir/sorting_hat/taxdump/` (only `names.dmp`, `nodes.dmp`, `taxdump.tar.gz`). It is inside the tarball and maps 502780 to 121759 and 502779 to 1048829; both targets are in `nodes.dmp`. Fix: say that B1 reads `merged.dmp` from `taxdump.tar.gz`.

**m-C. Fisher example uses the old background.** The example "10% of 150 background genera gives p = 4.2e-7" refers to the rev 1 design. With 610 background genera and 61 present: p = 1.5e-7 (two-sided equals one-sided here). The point of the sentence stands. Fix: update the numbers.

**m-D. Subsample stratification.** Proportional allocation over 610 genera gives Sordariomycetes 217/610 = 35.6%, so the 40% cap never binds. Five genera have an empty CLASS and are still not handled (review 1, m-5). Classes with 1 to 5 genera get fractional allocations below 1; the rounding rule is not stated. The stability criterion (for example, reported in at least 4 of 5 subsamples) and the seeds are not stated. Fix: state all four.

**m-E. Depletion and the reporting rule.** The test is two-sided, but the reporting rule needs at least 3 target genera present (or above the background median). A depleted family has few target genera present and is mostly not reported. Fix: say whether depletion is reported, and with which rule.

**m-F. BH scope.** "BH over all families tested for each pool and target" does not say whether Fisher and Mann-Whitney share one BH family or have one each. Fix: state it.

**m-G. Mann-Whitney on copies per 1,000 includes the zeros.** For a family absent in most genera, the test mostly repeats the presence test, with heavy ties at zero (an exact test is not possible; a tie-corrected normal approximation is). "Hits per 1,000" can count several domains in one protein. Fix: count proteins with at least one hit, state the tie handling, and optionally add the test among genera where the family is present, so expansion is separate from presence.

**m-H. Joins and formats for reused tables.** Tree tips are the input file stems (5,812 of 5,813 match; the missing one, `Herpotrichiellaceae_sp._CTeuk-1634`, has an empty genus), not ASMID or genus; the spec says neither. Fungi_5k TMHMM is in `key=value` short format (`PredHel=0`), not the tool's `tmhmm.tsv` header format, so a converter is needed. MEROPS results are raw BLAST tabular (`*.blasttab.gz`, no filter applied; 5,815 entries because some are also uncompressed), so a hit threshold must be stated before a MEROPS class is assigned. Fix: state the tip key, the TMHMM conversion, and the MEROPS and CAZy (dbCAN `overview.tsv.gz`) thresholds in the plan.

## Items checked and correct

- `results/function/` contents: `pfam/` has 5,813 `*.domtblout.gz` (plus `tblout` and `log`); header shows `hmmsearch --cut_ga`, HMMER 3.4, run 2025-05-03; `PF10417.14` in the output, which supports Pfam 37.x as review 1 found. `predgpi/` has 5,813 `*.predgpi.gff3.gz`. `signalp/` has 5,813 `*.signalp.results.txt.gz` with header `SignalP-6.0 Organism: Eukarya`. `tmhmm/` has 5,814 `tmhmm_short` and `tmhmm_results` files. `cazy/` has 5,813 per-proteome directories (`cazymes.tsv.gz`, `overview.tsv.gz`, `substrates.tsv.gz`). `merops/` has BLAST tabular results. `targetP`, `wolfpsort`, `aiupred`, `repeatmasker` exist.
- B2a looks feasible. For *C. immitis* RS, the Fungi_5k SignalP 6 SP set and the scan's SignalP 6 (`--mode fast --organism eukarya`, 6.0h) SP set are identical (371 SP each, 0 differences). Fungi_5k TMHMM PredHel and topology equal the scan's `tmhmm.tsv` for all 7,630 proteins (0 differences). The scan FASTA has the same 7,630 IDs as the Fungi_5k input. This is one proteome; B2a must still test all 56.
- Named-genus counts in section 3: Onygenales 27 named, 26 pass (only *Paracoccidioides* fails, on part 3); Aspergillaceae 7 named, 7 pass (*Aspergillus*, *Evansstolkia*, *Hamigera*, *Monascus*, *Penicilliopsis*, *Penicillium*, *Xeromyces*); other Eurotiales 7 pass. 14 Onygenales genera and 3 Aspergillaceae genera have a second passing proteome (limit 6.7).
- *Paracoccidioides*: both proteomes have more than 4,000 proteins; taxids 502780 and 502779 are absent from `nodes.dmp`; `merged.dmp` in the tarball maps them to 121759 and 1048829, which are present.
- Every proteome that the max-count rule would choose for a passing genus (targets and background) is a tip of `final_tree.nw`.
- Job arithmetic: `sacct` for 29681841 (`fs_pfam_1`, 1:24:12), 29681843 (`fs_repeats_1`, 0:38:32), 29681850 (`fs_tmhmm_9`, 1:07:13), 29681858 (`fs_signalp_29`, 0:25:35), all COMPLETED. Per-proteome rates as in review 1. 5 min x about 660 to 684 proteomes = 55 to 57 h, which is 37 to 46 jobs at 1.25 to 1.5 h. The statement that running TMHMM too would exceed the job limit is correct (8 to 11 min x 660 = 88 to 121 h more). The number of reused scan proteomes is uncertain, because the new choice rule may not pick the scan proteome for a genus; the effect is at most about 4.7 h.
- Section 1 medians: Onygenales 412.5 of 427, *Aspergillus* 784 of 808.5; summed ratio 0.967.
- Two-sided Fisher with the full background: 7/7 vs 0/610 gives p = 1.5e-16; 3/7 vs 0/610 gives 9.0e-7. Both pass a BH threshold of 0.05/3000 = 1.7e-5 for the top rank, so the test has resolution for Aspergillaceae. This is consistent with the spec's statement that the p-value is a ranking statistic.
- `run_list.tsv`, `scan_step.sbatch` and `scan_converters.sh` exist and use `${SCRATCH:?}` through the called job scripts; the Pfam step checks that the resolved path is in a `Pfam38.2` directory.
- The spec states associations, not function, in sections 1 and 6, and the role-table rules in 5.4 meet M-9.
