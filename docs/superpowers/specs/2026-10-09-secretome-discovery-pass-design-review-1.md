# Review 1 of the secretome discovery pass spec (rev 1)

*2026-10-09. Independent review of `docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md`, revision 1. Reviewer: Claude Code (claude-opus-5-5). The spec was not edited. All numbers below come from the files and commands named in each finding. Fisher p-values were computed with an exact one-sided hypergeometric sum (`math.comb`), not taken from a library.*

## Verdict

**FAIL.** One blocker, nine major findings, nine minor findings.

The blocker is cheap to fix (re-run the converters and the core command on the 56 scan proteomes). Several major findings also have a cheap fix: Fungi_5k already holds full Pfam-A, SignalP 6, TMHMM and PredGPI results for all 5,813 proteomes, a species tree, and MMseqs2 clusters. The spec does not mention them.

## Findings

### BLOCKER

**B-1. The 56 scan outputs were made with a different active family set. The spec reuses them "as they are".**

- Evidence: the scan module tables are dated 2026-10-09 00:48 (`_workdir/fungi_scan/Coccidioides_immitis_RS/modules/`). `pfam_adhesion.json` lists `"families": ["PF05730", "PF07691"]`. At commit `3ff1d5e` (the scan report), `data/sorting_hat/family_table.tsv` had 10 active families. Today it has 18 active (commits `cbfdc45` and `d94a8ab` came after the scan). The 8 more active families are all `pfam_adhesion` models: Candida_ALS_N, Candida_ALS, Flo11, GLEYA, Flocculin, Flocculin_t3, Hyr1, PIR1-like_C.
- Why it blocks: the reused proteomes are mostly target genera (24 Onygenales, 12 *Aspergillus*, 4 Penicillium/Talaromyces). The background proteomes would be new runs with 18 active families. So the "unexplained" pool has a different definition in target and background. This difference lines up with the target/background split, which is the comparison under test. Spec section 6.2 says the active set must be recorded. Recording it does not remove the difference.
- Fix: delete the sentence "Their outputs stay as they are" (section 3). Re-run `analysis/fungi_scan/run/scan_converters.sh` and the core command on the 56 proteomes with the current family table. The raw domain tables already hold hits for all 24 models (for example, PF30910 and PF04681 are in `raw/pfam/domtbl.txt` of *C. immitis* RS), so no new hmmsearch job is needed. Add a check to B2: every proteome's `run.json` and `pfam_*.json` list the same active set.

### MAJOR

**M-1. Section 3 counts the empty genus as a genus. *Paracoccidioides* fails the selection rule for a reason the spec does not give.**

- Evidence (`samples.csv`, recomputed): Onygenales has 28 distinct GENUS values, and one of them is the empty string (3 proteomes: "Onygenales sp. PD_10", "PD_12", "PD_40"). So there are 27 named genera. Aspergillaceae has 8 values, and one is empty ("Aspergillaceae sp.", "Eurotiomycetes sp. MA 6039"), so there are 7 named genera. The rule of `select_proteomes.py` also requires `NCBI_TAXONID` in `_workdir/sorting_hat/taxdump/nodes.dmp`. Both *Paracoccidioides* taxids (502780, 502779) are not in that file, although both proteomes have more than 4,000 proteins (7,041 and 7,251). Result: at most 26 named Onygenales genera and 7 named Aspergillaceae genera pass. *Paracoccidioides*, a primary dimorphic pathogen, is lost.
- Fix: count named genera only, and say how rows with an empty genus are handled (exclude them). State all three parts of the selection rule (FASTA exists, 4,000 proteins, taxid in the taxdump). Recover *Paracoccidioides* through `merged.dmp` or the species taxid. Change "28" and "8" in sections 2, 3 and 5.1.

**M-2. Genus presence by Fisher's test treats genera as independent. For Aspergillaceae the test can only find families that are almost restricted to the clade.**

- Evidence: one-sided exact Fisher test, background 150 genera. The BH threshold for the top-ranked family with about 3,000 tests is 0.05/3000 = 1.7e-5.

  | Target present / genera | Background fraction 0 | 0.1 | 0.3 | 0.5 | 0.7 |
  |---|---|---|---|---|---|
  | 7 / 7 | 2.5e-12 | 4.2e-7 | 3.3e-4 | 9.3e-3 | 8.9e-2 |
  | 5 / 7 | 2.8e-8 | 3.7e-4 | 3.4e-2 | 0.24 | 0.65 |
  | 3 / 7 | 5.5e-5 | 3.3e-2 | 0.36 | 0.77 | 0.97 |
  | 26 / 26 | 1.2e-31 | 7.4e-21 | 2.1e-12 | 1.1e-7 | 2.3e-4 |
  | 15 / 26 | 3.9e-15 | 2.1e-7 | 6.7e-3 | 0.30 | 0.93 |

  For 7 Aspergillaceae genera, a family in all 7 passes only when 10% or less of the background has it. With 3 of 7 it passes only when no background genus has it. One gain on the Aspergillaceae stem gives 7 of 7 and p = 2.5e-12. That is one evolutionary event counted seven times. The "3 independent genera" rule (owner decision 1) does not make the genera independent.
- Where it misleads: the reported q-values will look like strong evidence. In fact they measure clade restriction plus shared ancestry. Section 6.1 names the confounding. Section 5.1 still reports BH q-values with no warning next to them.
- Fix: (a) In every output table, label the Fisher p and q as a ranking statistic under an independence assumption that does not hold. (b) For each reported family, count the independent origins on the genus tree. `Phylogeny/fungi_tree/final_tree.nw` in Fungi_5k exists, and parsimony on presence/absence is enough. Report that count next to the q-value. A family with one origin is "clade-specific", not "found in N genera". (c) Optional: use a phylogenetic logistic regression on the same tree as a second test. Pick one of these before B4, because each choice changes what B4 tests.

**M-3. Presence/absence ignores copy number and secretome size.**

- Evidence: medians of `signal_peptide_protein[R0]` in `proteome_summary.tsv` are 427 (Onygenales, 24) and 808.5 (*Aspergillus*, 12). A larger secretome raises the chance that any family is present, so the Aspergillaceae test is biased upward for common families. Lineage-specific expansions (copy-number gains with presence in all genera) give odds ratio 1 and are never found. Those expansions are the usual genomic signal for host-associated adaptation.
- Fix: add a count test per family (per-genus copy number in the pool, with log pool size as an offset, for example a negative binomial or Poisson GLM; or Mann-Whitney on copies per 1,000 pool proteins). Add log pool size as a covariate to the presence test (logistic regression), or report the Fisher result stratified by pool size. Report depletion as well as enrichment, or state that the test is one-sided (see m-4).

**M-4. The "unexplained" pool removes the protein classes that the three purposes care about most. The candidate lists cannot include proteins without a Pfam hit.**

- Evidence: `docs/CLASSES.md` defines `other_surface_no_mechanism` as signal peptide and no repeat, wall family, hydrophobin, HsbA or antigen-rank evidence. So every secreted protein with a tandem repeat call leaves the pool. Ser/Thr-rich repeat glycoproteins are the classic surface and immune-masking class. The repeat14 coverage feature in 5.2 is therefore close to zero by construction. The pool also keeps all proteins that have a signal peptide plus transmembrane helices, and all ER/Golgi/vacuolar residents. TMHMM is already run but the pool does not use it. Section 5.4 makes candidates only from 5.1 (Pfam families). Section 6.3 says proteins with no Pfam hit become a count and a composition summary only. Proteins of the SOWgp/Cocci-antigen kind, which have no Pfam domain, can never become candidates.
- Fix: (a) Run the family test on two pools: the unexplained pool, and the full signal-peptide pool. A family that is only enriched in the second pool is then visible. (b) Remove proteins with a mature TM helix (the CFEM rule already uses `n_tm_mature`) and C-terminal HDEL/KDEL from the pool, or report them separately. (c) Add a Pfam-free family unit: Fungi_5k has MMseqs2 clusters (`results/Fungi5K_cluster.tsv.gz`, `results/orthogroups_mmseqs.tsv.gz`) and a no-Pfam cluster dump (`results/mmseqs_nopfam_genedump.*`). Run the same genus test on clusters with no Pfam hit, and let them enter 5.4.

**M-5. The tests answer "what is clade-enriched", not "what serves host interaction".**

- Evidence: the 26 Onygenales genera include saprotrophs and keratinophiles (*Onygena*, *Amauroascus*, *Malbranchea*, *Uncinocarpus*), the bee pathogen *Ascosphaera*, dermatophytes, and the dimorphic pathogens. The 7 Aspergillaceae genera include *Monascus*, *Xeromyces*, *Hamigera* and *Penicilliopsis*, which are not known as pathogens of mammals. A family enriched in either clade is a lineage marker. The design gives no step that ties a family to host interaction or immune evasion.
- Fix: add a descriptive contrast inside Onygenales: genera with mammal-pathogen species against the other Onygenales genera. Write the pathogen list in the plan before the run, with a source for each genus. The n is small, so make it a descriptive table, not a test. Say in section 1 that clade enrichment cannot tell host interaction apart from other lineage traits.

**M-6. The 5-seed background design is not costed. It roughly triples the compute.**

- Evidence: 603 named Pezizomycotina genera outside Onygenales and Eurotiales pass the selection rule (recomputed). Five independent draws of 150 cover about 603 x (1 - (453/603)^5) = about 459 distinct genera. Scan rates from `sacct` (jobs 29681841 to 29681858): TMHMM about 8 to 11 min per proteome (8 per job, 0:53 to 1:30), repeats about 2 min (about 19 per job, 0:37 to 0:42), allergen about 1.5 min (28 per job, 0:41 to 0:45), SignalP about 0.9 min on GPU (28 per job, 0:23 to 0:26), family-subset Pfam about 1.5 min (56 in 1:24). That is about 15 to 18 min per proteome. One draw: about 34 target + 150 background + 40 supplement = about 224 proteomes, less the reused ones, so about 55 to 65 h, or about 40 to 50 jobs at 1.25 to 1.5 h. That is already at the 40 to 60 limit before the full Pfam-A and PredGPI steps. Five draws: about 460 + 74 proteomes, so well over 60 jobs.
- Fix: use the Fungi_5k precomputed results (M-7) for TMHMM, SignalP and PredGPI, which removes the heaviest steps. Then run the repeat and family-subset steps on all 603 pool genera once (about 3.5 min per proteome, so about 35 h or about 25 jobs). Draw the 5 backgrounds in the analysis from that fixed pool. A full background of 603 also removes the seed dependence and increases power. Write the job count in the spec, not only in the plan, because it decides the design.

**M-7. The spec does not use the Fungi_5k annotations that already exist.**

- Evidence: `/bigdata/stajichlab/shared/projects/Fungi_5k/results/function/` holds `pfam/*.domtblout.gz` for 5,813 proteomes (`hmmsearch --cut_ga`, full Pfam-A, run 2025-05-03), `predgpi/*.predgpi.gff3.gz` for 5,813 proteomes, `signalp/` (SignalP 6.0) and `tmhmm/` results. The Pfam release is 37.x, not 38.2: PF10417 is version .14 in the Fungi_5k output and in Pfam 37.0 and 37.3, and version .15 in 38.2. `Phylogeny/fungi_tree/final_tree.nw` is a species tree.
- Fix: decide in the spec whether to (a) reuse the Fungi_5k full Pfam-A 37.x and PredGPI results, and run only the 38.2 family subset (needed for the mechanism exclusion, because PF28987, PF29785, PF29802, PF29465 and PF30910 are newer families), or (b) run Pfam-A 38.2 again with a reason. If SignalP or TMHMM is reused, show that the Fungi_5k run gives the same R0 calls as `signalp_gpu.sbatch` on the 56 scan proteomes before use (same mode, same organism setting).

**M-8. The supplement groups by "species group as given in `samples.csv`". That column does not exist.**

- Evidence: the `samples.csv` header is `ASMID,SPECIESIN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,LOCUSTAG`. No section or species-group column exists. The 228 *Aspergillus* proteomes have 228 distinct SPECIES strings.
- Fix: name the source of the sections (a published section list, with a citation, joined by species name), or use clades of the Fungi_5k tree. Otherwise drop the grouping. State how the 40 species are chosen (for example, a fixed number per section).

**M-9. The T4 role table is written after the results are known. Nothing limits post-hoc stories.**

- Evidence: section 5.4 says the role table is "in this spec's plan" and maps "Pfam descriptions to the three purposes". The families exist only after B4. Pfam descriptions seldom state host interaction, and almost never immune evasion. The T4 label is stated once in the spec, but the spec does not require it on each output row.
- Fix: (a) Allow and default to "no stated role" for each family. (b) Every role row cites its basis: the quoted Pfam/InterPro text, or a paper (PMID). A row that only has assistant reasoning gets the T4 label in the row itself. (c) Write the role table blind to the enrichment statistics where possible: write roles for the family list sorted by accession before the q-values are attached. (d) Put a sentence on the report's first page: the role column is assistant judgement, not evidence.

### MINOR

**m-1. Medians in section 1 do not match the cited report table.** The report table (`docs/reports/2026-10-09-fungi-scan-descriptive.md`) gives signal-peptide medians 438 (Onygenales) and 822 (Eurotiales), which are upper medians. The true medians from `proteome_summary.tsv` are 427 and 816.5. The spec's 427 is correct, but 413 is 412.5 and 809 is 808.5 (rounded). The "*Aspergillus*" medians (12 proteomes) are not in the report; they come from `proteome_summary.tsv`. Fix: cite `proteome_summary.tsv` and the median convention. Also note the report's convention.

**m-2. "The seven module steps" is not correct.** `scan_step.sbatch` has five steps (pfam, signalp, tmhmm, repeats, allergen). `scan_converters.sh` then writes eight module tables. Fix: name the steps as in the script.

**m-3. The filter before BH is not independent of the test.** Filtering by "3 target genera present" uses the group label, so BH on the filtered set is not valid. BH over all families and then a reporting filter is valid. Fix: say BH runs over all families with at least one present genus (or filter on pooled presence, which does not use the label), and the 3-genus rule only limits what is reported. State the number of families tested from the data. "About 3,000" has no source.

**m-4. The sidedness and the sister comparison are not defined.** State one-sided (enrichment) or two-sided Fisher. "Onygenales against the other Eurotiales" is ambiguous: 7 non-Aspergillaceae genera or all 14 named Eurotiales genera. With 7 background genera no test can pass BH. Call the sister comparison descriptive.

**m-5. The class stratification is not defined.** Passing named genera per class: Sordariomycetes 217, Dothideomycetes 176, Leotiomycetes 76, Lecanoromycetes 60, Pezizomycetes 24, Eurotiomycetes 19 (non-target), and 10 classes with 1 to 5. "No class dominates" can mean equal or capped allocation. Fix: state the allocation rule (for example, proportional with a cap) and how the 5 genera with an empty class are handled.

**m-6. Proteome quality is not controlled.** One proteome per genus is chosen as the first by species name. Gene-model errors at the N-terminus remove signal peptides and so change the pool. Fix: choose by BUSCO completeness or by protein count within the genus, or record BUSCO for each chosen proteome and report it next to presence/absence. Optional: repeat the analysis with a second proteome per genus where one exists.

**m-7. The antigen module is absent in the scan runs.** `run.json` of *C. immitis* RS lists `"absent_modules": ["antigen_lookup"]`. The antigen call works only for *C. immitis* RS. Section 5.4 promises "any hit to the antigen ... modules". Fix: say that this column applies to *C. immitis* RS only.

**m-8. The HPCC rules for new steps are not stated.** The existing job scripts use `${SCRATCH:?...}` and do not use `BASH_SOURCE[0]` (checked: `scripts/sorting_hat/*.sbatch`, `analysis/fungi_scan/run/*`). The spec does not require the same for the new full Pfam-A, PredGPI and composition jobs. Running full Pfam-A per proteome on about 500 sequences is inefficient. Fix: state `$SCRATCH` and absolute paths for new jobs. Pool the unexplained sequences of all proteomes into one FASTA and split it into 1 to 1.5 h chunks. Write per-proteome tables as `.tsv.gz`. Pin the resolved Pfam path (`2026-01-27-Pfam38.2/Pfam-A.hmm`), not `current`.

**m-9. The 5.2 "fixed level" is chosen in the plan.** Fix the levels before the data are seen, and write them in the spec or the plan before B3. The Mann-Whitney test with 7 Aspergillaceae genera has enough resolution, but the effect size must be given with the genus counts.

## Items checked and correct

- `samples.csv` has 5,813 data rows.
- Section 3 table, as raw distinct GENUS strings: Onygenales 71 proteomes and 28 genera; Aspergillaceae 406 and 8; Eurotiales 483 and 15; Pezizomycotina 2,802 and 654; other Ascomycota 1,255 and 125; Basidiomycota 1,294 and 387. (The genus counts include the empty string; see M-1.)
- *Aspergillus*: 228 proteomes in one genus.
- Eurotiales genera outside Aspergillaceae: 7 named (*Talaromyces*, *Paecilomyces*, *Elaphomyces*, *Pseudotulostoma*, *Rasamsonia*, *Thermoascus*, *Thermomyces*).
- 97%: the ratio of summed `other_surface_no_mechanism[R0]` to summed `signal_peptide_protein[R0]` over the 56 proteomes is 0.967. The median ratio is 0.970.
- 18 of 24 families are active in `data/sorting_hat/family_table.tsv` today.
- `other_surface_no_mechanism` in `docs/CLASSES.md` is "Signal peptide and no mechanism evidence". It already requires the signal peptide, so the AND in section 4.2 is redundant but harmless.
- Scan timings: TMHMM 0:53 to 1:30 for 8 proteomes per job; SignalP 0:23 to 0:26 for 28 per job (`sacct`, all `COMPLETED`).
- Pfam-A 38.2 is the file the family scan used (`raw/pfam/provenance.json` resolves to `2026-01-27-Pfam38.2/Pfam-A.hmm`; `hmmsearch --cut_ga`, HMMER 3.4).
- The module `predgpi/202001` exists on the cluster.
- The existing job scripts use `$SCRATCH` with a fail-fast check and do not use `BASH_SOURCE[0]`.
- The spec states associations, not function. It names phylogenetic confounding (6.1), the dependence on the active set (6.2), the Pfam coverage limit (6.3) and the lack of validation of PredGPI (6.5). These statements are correct. Findings M-2, M-4 and B-1 say where the design does not yet act on them.
- The T4 label for the role table is stated in the spec, and the spec says the table is not a truth set.
