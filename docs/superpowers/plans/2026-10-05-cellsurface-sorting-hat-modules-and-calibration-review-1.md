# Review 1 of Plan 2 (modules and calibration), with dispositions

*Five independent reviews on 2026-10-05, each on the first draft: general and replay (Fable), bioinformatics
(Opus), computational biology and statistics (Opus), immunology (Opus), microbiology (Opus). Counts per
reviewer: general 1 blocker, 6 major, 17 minor; bioinformatics 2, 9, 10; statistics 3, 7, 5; immunology 4, 12, 7;
microbiology 1, 9, 9. All five verdicts: needs rework. The owner decided four questions (below). The plan was
then revised. I checked the revision myself: all 37 file blocks and the Plan 1 patch were extracted into a
clean directory (256 tests pass, 7 shellcheck tests skipped), ruff 0.3.5 is clean, twelve mutation checks make
the named tests fail, and the Pfam key extraction was run against the real database (15 models fetched; the names
and accessions match the family table). The revised plan has **not** had a second independent review. Literature
claims that the reviewers gave from memory were not verified.*

## Owner decisions taken in this round

| Question | Decision |
|---|---|
| Output names that claim more than the tools measure | Rename all (table "Changes to the spec" in the plan). |
| Gene models for the Af293 calibration and run | UniProt UP000002530 first; Fungi_5k as a second annotation. |
| Use of the antigen list | Serodiagnostic marker candidates; the gated column is the headline. |
| Allergen call | Applied the reviewers' common recommendation (two evidence flags, not a mechanism category). The owner was not asked separately; reverse in the config if not wanted. |

## Dispositions (grouped; A = applied, P = left for a later plan or owner decision, N = noted)

| Finding (reviewer) | Disposition |
|---|---|
| `hmmfetch -f` fails on unversioned accessions (general, bioinformatics) | A. Fetch by model name; names and accessions checked against the database; tested on the real database. |
| Pfam and allergen modules write all-zero tables when the tool output belongs to another proteome, is empty or truncated (bioinformatics, general) | A. Domain table needs `# [ok]` and `--cut_ga`; all targets and queries must be FASTA IDs; a table with no usable result is refused; `partial` and `unavailable` run states. |
| While no family is active the domain test is "done" with zeros (general) | A. `unavailable` modules. |
| Allergen leave-cluster-out recall is zero by construction; cutoffs do not match the engine's rules; denominator from BLAST (statistics, bioinformatics, immunology) | A. Replaced by leave-species-out recall with the engine's own rules and the FASTA denominator; refuses a sequence with no BLAST line. |
| 35%/80 aa is one local alignment, not the sliding window; "FAO/WHO rule" label (bioinformatics, immunology) | A for naming (`iuis_allergen_similarity`, text says single alignment). P for a sliding-window implementation. |
| Allergen FASTA builder: whitespace, free-text entry, fragments, wrong count (bioinformatics, statistics) | A. Cleaning, skip list, meta table with species, exposure, evidence. Real result: 111 written, 5 skipped. |
| Names claim more than the tools measure (immunology, microbiology, bioinformatics) | A. Renames; research-use header and limits in the report. |
| Pooled R0 estimate attached to species it does not describe; Phase C converter overwrites the computed status (statistics, microbiology, general) | A. One entry per species; status is the weaker of the Phase C label and the rule; per-stratum specificity kept; zero-width intervals widened. |
| `estimated` passes zero-width intervals and counts proteins, not clusters (statistics, general) | A. Hybrid interval (bootstrap or Wilson on clusters, whichever is wider), cluster counts recorded and required (when recorded). |
| Calibration on a different annotation than the truth (microbiology, bioinformatics) | A. UniProt first; Fungi_5k second with R0 agreement reported; annotation named in the claims. |
| Stale status entries kept when the module identity changes (general) | A. Dropped with a message. |
| Tool identity could not tell builds, scripts or options apart (bioinformatics) | A in part: tool version, `--seg`, detector script hash, hmmer version are recorded. P: pin exact module builds and hash model weights. |
| `no_tm` drops proteins whose signal peptide TMHMM reads as a helix (microbiology, bioinformatics) | A. `n_tm_mature` ignores helices starting at or before residue 35. The 35 is a choice, not a measurement. |
| Repeat calls are mostly intracellular; proteins under 80 aa become `error` (bioinformatics) | A. Repeat calls are evidence; the adhesion call needs a signal peptide; short proteins are `not_called`. P: a Pfam gate for ubiquitin, EF-hand, ankyrin, WD40, zinc-finger repeats. |
| Ungated antigen call mostly not secreted (microbiology, immunology) | A. Gated headline column; renamed; header text. |
| Antigen lookup uses gene then best transcript (bioinformatics) | N and gap. `idmap_method` column added; exact sequence match not implemented (126 RS isoforms affected). |
| Antigen cut of 15 unsupported by four anchors (immunology) | A in the wording (table, report). P: leave-species-out antigen test (owner decision C3). |
| Self-recognition in the allergen panel and Af293 counts (immunology) | A. `tuning = in_reference` excludes such rows; counts split into same-species and other-species. |
| Allergen negatives could come from IEDB negative IgE assays (immunology) | P. Owner decision C2; the number of such records was not checked. |
| Specificity of the Pfam test is computed from unlabelled proteins (statistics, bioinformatics, microbiology) | A. The command prints the numbers and the plan says not to report them as specificity; reviewed non-member hits are classified. |
| Strain stability tests only determinism (statistics, microbiology) | A in the plan (best-reciprocal-hit comparison by hand, named cases). P: code for reciprocal best hits. |
| No PPV or expected false calls (statistics) | A in the report steps (by hand, assumed prevalences). P: store in `measure`. |
| Leakage is self-declared; SignalP 6 training overlap unmeasured (statistics, microbiology) | A in part: `in_reference` value; R0 recorded as `leakage = unknown`. P: automatic overlap measurement. |
| Sn bound with not-assessable counted as missed (statistics) | A. |
| Panel: missing proteins counted as disagree; wrong sentence (general) | A. `not_in_run`, `tuning` column. |
| Plan text: packaging order, counts, Step 2 text, `$SCRATCH` shell, PYTHONPATH, S288C dubious ORFs, venv (general, microbiology) | A. |
| Family table lacks yeast adhesin domains (microbiology) | A. Five inactive rows added (Flo11, GLEYA, Hyphal_reg_CWP, PIR, PA14 with a signal-peptide condition). |
| Allergen evidence: exposure route, species, IUIS evidence text (immunology) | A. Carried through the module table. |
| Scale to 5,000 proteomes; per-sha256 cache not used; Kind K sha256 match; spec 3.8 retry (general, bioinformatics) | Listed as gaps in the plan. |
| Literature points from memory (microbiology: Vaknin 2014, Liu 2016, Balajee 2007, Fedorova 2008, Gravelat 2013, A1163 lineage, gp43) | N. Not verified; not relied on in the code. Check before citing. |
| `status_from_measure` and Plan 1 definition of `estimated` differ (statistics) | A in the plan text; the spec sentence is updated in the pull request. |

## Not done

- A second independent review of the revised plan.
- The sliding-window allergen search, exact sequence-hash ID mapping, a repeat-family gate, SignalP training-set overlap, IEDB negative assays and the leave-species-out antigen test.
- Verification of the taxon IDs 330879 (Af293), 451804 (A1163) and 746128 (W72310) against a taxdump; Task 11 checks them.
