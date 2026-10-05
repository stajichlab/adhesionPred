# Review 1 of the orchestrator design spec (revision 2), with dispositions

*Reviewer: an independent Opus subagent, 2026-10-04. It read the spec, the plan, the architecture and
status documents, the handoff, the analysis code and outputs, `data/curated/`, `src/surface_glyco/` and
issues #14, #16 and #19. It changed no file and submitted no job. Verdict: **needs rework before a plan.**
Count: 2 blocker, 13 major, 10 minor.*

*I spot-checked these claims before accepting the review: SOWgp has prevalence 0.9201 and no
`fungi5k_id` in `cocci_antigen_ranking.tsv`; PRA2 fails the acceptance test (3 of 4 anchors pass);
`src/surface_glyco/models/` holds only a README; the plan separates two meanings of "other"; TMHMM
(`01b_tmhmm.sh`) and the six-model HMM script (`01_known_family_hmm.sh`) exist; issue #16 gives 60 to
150 GPU-hours for ESM C 300M. The other findings I did not re-derive; the dispositions below say
where revision 3 relies on them.*

Disposition codes: **A** accepted and applied in revision 3. **O** needs an owner decision
(listed as D10 to D12 in the spec). **N** noted, no change.

| # | Sev | Finding (short) | Disposition |
|---|---|---|---|
| 1 | blocker | Three-valued logic covers AND only. "Other" is ill-defined and ignores that the plan has two meanings. | A. Kleene logic with truth tables (spec 3.3). `other` split into `other_not_surface` and `other_surface_no_mechanism`. |
| 2 | blocker | Neither step 1 module can run: no ML model ships, R1/R2 have unfrozen parameters. Choosing a module in a category rule is choosing step 1. Calls on RS differ much: R0 460, R2 94, M8 V-go 758, M8 V-kw 5,706. | A in design (every gated call is computed per step 1 variant; a variant with no frozen artefact is `unavailable`). O for the default gate (D10). |
| 3 | major | The `run(fasta)` interface does not fit the antigen and expression layers, which read Fungi_5k DuckDB, ID maps, kallisto TPM. | A. New module kind K (lookup by ID, with an ID-mapping step and `not_in_reference`). |
| 4 | major | `antigen_candidate` cannot call SOWgp (prevalence 0.9201) or Ag2/PRA; "top tier" is undefined; R2 keeps 5 of 14 Tier 1 and 19 of 45 Tier 2. | A. Antigen call from the ranking percentile, not the tiers. No step 1 gate in the base call; gated call is a separate column. |
| 5 | major | The R2 gate costs recall (R2 misses 58%, 77%, 88% of surface proteins in S1, Eurotiomycetes, Basidiomycota) and risk 2 is worded wrongly. | A. Risk 2 rewritten. Gated and ungated calls are both written. O for the default gate (D10). |
| 6 | major | "Validated in Saccharomycotina only" belongs to the ESM classifiers, not the repeat detectors. The detectors have a synthetic benchmark plus SOWgp only. | A. Detector status is `unvalidated`. |
| 7 | major | Status is per clade, not per module. Clade labels are four class-level strings; `--clade` is not checked; one clade per run cannot score a multi-clade panel. | A. Status is a function of (module, version, clade). NCBI taxon IDs with lineage matching. Per-protein taxon map. Scope rule stated. |
| 8 | major | "Agreement not yet measured" is out of date: Phase C measured R2 against ML per proteome. | A. Cite the Phase C tables. |
| 9 | major | `cys_rich_secreted` contradicts its README ("no accuracy claim", "does not mean PRA3-like") and omits SignalP. | A. Moved out of the adhesion category into an evidence column; scope *Coccidioides*; SignalP required. |
| 10 | major | Cache key misses parameters, model hashes, the categories hash. The funannotate Pfam copy has no release file. | A. Cache key and provenance rules. Central Pfam path with a recorded release. Atomic writes. |
| 11 | major | The test panel leaks: the panel proteins tuned the modules, and homologs cross the split. | A. Per (module, protein) tuning flags; homology-cluster split. |
| 12 | major | Acceptance is not testable. No category outside step 1 can reach the `estimated` rule. A golden test needs GPU SignalP. | A. v1 acceptance is software correctness on stored fixtures. Per-category numbers are report-only. |
| 13 | major | Gaps against the five categories: structural wall proteins (Cwp1, Ccw12, Sed1, Pir) get only `surface_glycoprotein`; allergen was a placeholder; allergen surface policy unstated. | A for allergen (IUIS homology module, no surface gate). O for a `cell_wall_protein` sub-call (D11). |
| 14 | major | Compute claims wrong: 488 is the *Coccidioides* isolate count; issue #16 is about 58 M proteins, ESM C 300M. | A. |
| 15 | major | Failure modes missing; `not_assessable` mixes out-of-scope with failure. | A. Failure-mode table (spec 3.8). Separate run states. |
| 16 | minor | Curated table counts and `needs_review` shares are loose. | A. Counts from the review, marked as not rechecked. |
| 17 | minor | "3 of 4 anchors" reads as 3 failures. | A. |
| 18 | minor | `01_known_family_hmm.sh` already runs the HMMs; PF28987 and PF22354 missing from the spec. | A. |
| 19 | minor | `validated` reads as good; "called (unvalidated)" breaks parsing; `status` is ambiguous. | A. `estimated`; separate `_status` columns. |
| 20 | minor | TMHMM was also run; v1 needs TM evidence (R2 calls MSB2 and HKR1). | A. TM column noted as evidence in v1. |
| 21 | minor | Panel: Hsp60 has no expected value; hard negatives rest on a seed list. | A. |
| 22 | minor | Expression conflicts with the non-goal. | A. Evidence column only. |
| 23 | minor | CFEM filed as adhesion conflicts with the hemophore fold; `surface_glycoprotein` truth is GO-based. | A. Report header notes. |
| 24 | minor | Kappa is unstable at low prevalence. | A. Also report 2x2 counts. |
| 25 | minor | Section order, revision label, timing source. | A. |
