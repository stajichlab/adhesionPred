# Review 2 of the orchestrator design spec (revision 4), with dispositions

*Reviewer: an independent Fable subagent, 2026-10-04 (a different model from the author and from
reviewer 1). It checked the review 1 fixes, the new content and the owner's decisions, and ran one
BLASTP search for the allergen cutoff. It changed no repository file and ran no SLURM job. Verdict:
**needs rework before a plan.** Count: 3 blocker, 9 major, 10 minor.*

*I did not re-derive the BLASTP numbers (108 Af293 proteins at 35%/80 aa; 79 not known allergens).
Counts of proteomes and files in findings 7 were checked by me on 2026-10-04: RS RefSeq FASTA 9,910
proteins; Fungi_5k RS file 7,630; Fungi_5k Af293 file 9,161.*

Disposition codes: **A** applied in revision 5. **O** owner answered a question in this round.
**P** left for the plan (stated in the spec).

| # | Sev | Finding (short) | Disposition |
|---|---|---|---|
| 1 | blocker | Scope and status are merged. Unvalidated modules would give U for every protein. | A, O. Applicability split from status (spec 3.2). D7 revised: `--taxon` required; not applicable gives U; applicable but unmeasured gives a value with status `unvalidated`. |
| 2 | blocker | `other_*` is never T outside *Coccidioides*, because antigen is U. Text does not match formulas. | A, O. `other_*` defined over the assessable mechanism calls, with an `other_basis` column. Formulas written (spec 3.4). |
| 3 | blocker | The 35%/80 aa allergen rule would call about 79 housekeeping paralogs in Af293; the scoping note advised against it. | A, O. D6 revised to two tiers: evidence `allergen_homolog_hit` (35%/80) and call `allergen_candidate` (>= 70% identity, >= 80% coverage, or allergen Pfam). Cutoffs are config items chosen without a non-allergen set. |
| 4 | major | Lineage matching not defined; "Eurotiomycetes" would pass `estimated` to Onygenales. | A. Status applies only to tested taxa and their descendants; most specific wins; taxonomy dump version recorded (spec 3.2). Unit tests added (spec 4). |
| 5 | major | `--taxon` and `--taxon-map` interplay undefined. | A. Map overrides `--taxon`; at least one required; protein with no taxon is an error. |
| 6 | major | Protein key and identical sequences. | A. Key is (ID, sha256); per-sha256 caching; test added. |
| 7 | major | D13 proteome facts incomplete (RS has two FASTA files; Af293 IDs differ between files; no S288C path). | A. Paths and counts in spec 4.1. S288C protein count not checked. |
| 8 | major | Kind K match thresholds missing. | A, P. Exact hash, then ID map, then fuzzy match at >= 95% identity and >= 90% mutual coverage as proposed defaults. To be fixed in the plan. |
| 9 | major | Cache and `status_source` check version only. | A. Check version, `params_hash`, `artefact_hash`; tool versions in key; batch settings for ML (spec 3.6). |
| 10 | major | Column count and unavailable variants. | A. Long format primary; wide derived; only available variants get columns. |
| 11 | major | Failure-mode gaps (trailing `*`, ambiguous residues, partial output, GPU OOM, wrong taxon, empty FASTA). | A. Rows added (spec 3.8). |
| 12 | major | "27 cases" wrong (9 AND + 9 OR + 3 NOT). | A. |
| 13 | major | No run-level acceptance on a real proteome. | A. Af293 end-to-end check on HPCC, not CI (spec 4). |
| 14 | minor | R0 means "has a signal peptide". | A. Sentence added to the known limits. |
| 15 | minor | P = 15 puts all four anchors inside the cut; weak label; `NOT CALIBRATED` note. | A. Known limit 5; report prints share called. |
| 16 | minor | Antigen outside *Coccidioides* should be "not applicable", not `not_in_reference`. | A. |
| 17 | minor | Reason for `percentile` over `percentile_dedup`. | A. |
| 18 | minor | `na_window` count per variant. | A. |
| 19 | minor | Denominator for "about one third have a signal peptide" (30 of 97). | A. |
| 20 | minor | Stale revision labels. | A. |
| 21 | minor | Driver waiting, resources, copy-back not specified. | P. Listed as required plan content (spec 3.6). |
| 22 | minor | `other_surface_no_mechanism` lists four mechanism calls; gated forms implied. | A. |

Missing for a plan, per the reviewer, and where they stand: scope/status redesign (done); lineage
algorithm (done in outline; the taxonomy source is an NCBI dump); allergen definitions (done);
K-module thresholds (proposed); output schema (done in outline); `categories.yaml` example,
module resource table, fixture list, golden file (plan); the five FASTA paths and counts (done,
except two proteomes not yet downloaded and the S288C count).
