# Design spec: the orchestrator (one table, one column group per tool)

*Drafted 2026-10-04. Revision 1. DRAFT. No independent review has run. No code, data or job exists
for this spec. Owner decisions are in section 9.*

Inputs: `docs/PLAN-2026-09-30-pipeline-and-decisions.md` section 6 (the proposal this spec
expands); `docs/TOOL-ARCHITECTURE.md`; `docs/model-review/STATUS.md` (2026-10-02);
`docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md` (the mapping from three cell wall
publications to the tools); `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`
(step 1 label definitions); `docs/HANDOFF-2026-10-03.md`; the code under `analysis/` named in section 2.

## 1. Purpose and non-goals

**Purpose.** The owner wants one entry point. It takes a proteome and writes one row per protein.
Each row says, per tool, what the tool found, how sure it is, and whether that tool has been
validated for this clade. A second layer turns the evidence into category calls. The categories
the owner named are: cell wall and adhesion candidate, allergen, antigen, glycoprotein on the cell
surface, and other. The cell wall publications add a further need: cell wall gene classes that are
defined by function (synthases, remodeling enzymes) and not by location (section 3.4).

**Non-goals.**
- It does not choose the step 1 method (rule, ML or hybrid). That is the owner's open decision. The
  orchestrator carries both candidates until the decision is made (section 3.2).
- It does not train a model or set gate values.
- It does not create truth data. Allergen curation (issue #19) and Onygenales curation come first
  (section 8).
- It does not replace any single tool. Every tool still runs alone.
- It does not predict cell wall chemistry, expression changes or signaling (section 3.4).

## 2. Facts checked

What exists on 2026-10-04 and what does not.

| Item | State | Where |
|---|---|---|
| Step 1, ESM-2 plus logistic regression | CLI exists. **No model ships.** Output columns `surface_glycoprotein_score`, `surface_glycoprotein`. | `src/surface_glyco/` |
| Step 1, SignalP 6 and PredGPI | Run as jobs on a protein FASTA (GPU job, 9,910 proteins in 2 min 54 s). Not wrapped. | `analysis/step1_compare/jobs/j1_features.sh`, `predgpi_scores.py`; `analysis/cocci_repeats/01_signalp.sh` |
| Step 1 rules R0, R1, R2 | Defined in the Phase C evaluation. Not exposed as a tool. | `analysis/step1_compare/phasec/` |
| Step 2a repeat detectors | Two scripts. Neither replaces the other. Validated in Saccharomycotina only. | `analysis/cocci_repeats/02_repeat_profile.py`, `14_repeat_detect_general.py` |
| Step 2b/2c HMMs | Models exist (PF05730, PF04681, PF01185, PF06766). No scan wrapper and no specificity test. | `docs/TOOL-ARCHITECTURE.md` section 2 |
| Pfam scan | `hmmsearch` over Pfam-A on HPCC (`/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm`, 30,134 models) is scripted for one project. Not a general tool. | `analysis/cocci_repeats/22_pfam_hmmsearch.sh` |
| Cys-rich candidates | A candidate finder with unit tests. Thresholds from 7 controls. | `analysis/cys_candidates/` |
| Antigen layer | *Coccidioides* only. Prints `NOT CALIBRATED` at 3 of 4 anchors. | `analysis/cocci_antigens/` |
| Expression layer | One data set (*C. immitis* spherule against mycelium, 2 replicates per condition). | `analysis/cocci_spherule/` |
| Allergen layer | **Does not exist.** No data directory. Issue #19 is parked. | issue #19 |
| Biofilm layer | **Does not exist.** Blocked on phenotype data. | `docs/TOOL-ARCHITECTURE.md` |
| Orchestrator | **Does not exist.** Proposal only. | plan section 6 |
| Curated tables | surface 3,540 rows; adhesins 269 lines; biofilm 207; antigens 87 lines. All drafts (`needs_review=yes` on nearly all). | `data/curated/` |
| Step 1 label classes | P-ext, P-gpi, N-int, N-sec, PM-TM, ambiguous. | step 1 spec 2.2 |

## 3. Design

### 3.1 Three layers

```
proteome FASTA
   |
   v
[ evidence modules ]  one run per module, cached by input hash
   L  localisation      SignalP 6, GPI call, transmembrane segments, step 1 rule and ML
   D  family / domain   hmmsearch against Pfam-A and a curated family table
   A  architecture      repeat detectors, composition (Ser/Thr, Pro, Cys), Cys-rich finder
   C  comparative       antigen ranking, allergen homology, expression
   |
   v
evidence.tsv.gz         long format: protein, module, field, value, status
   |
   v
[ category rules ]      config file, three-valued, multi-label
   |
   v
calls.tsv.gz            one row per protein, one column per category
report.md               counts, per-module status, what could not be assessed
```

### 3.2 Evidence modules

Each module is a plugin with the same small interface.

| Field | Meaning |
|---|---|
| `name`, `version` | for example `signalp6`, `6.0` |
| `kind` | L, D, A or C |
| `clade_scope` | clades where the module has been measured. Free text plus a machine list of taxon names. |
| `status` | `validated` (an estimate set exists: recall interval half-width at most 0.10 and at least 20 direct positives, the step 1 rule), `smoke` (fewer), `unvalidated` (no truth) |
| `status_source` | path to the measurement that supports `status` (a JSON written by the tier harness) |
| `run(fasta, workdir)` | writes a per-protein table |
| `fields` | names and types of the output columns |

Rules:
1. A module cannot set its own `status`. The status comes only from a measurement file. Without a
   file the status is `unvalidated`.
2. A module run on a clade outside its `clade_scope` writes the value and sets `in_scope = false`.
3. Step 1 appears as **two** L modules, `step1_rule` (R0, R1, R2 columns) and `step1_ml`, until the
   owner decides. The agreement of the two is reported (the plan lists it as not yet measured).

### 3.3 Category rules

A category is a boolean expression over evidence fields, stored in a config file
(`categories.yaml`), not in code. Each category lists the modules it requires.

The result for one protein and one category has three values:

| Value | When |
|---|---|
| `called` | the rule is true and every required module is in scope |
| `not_called` | the rule is false and every required module is in scope |
| `not_assessable` | at least one required module is out of scope or missing |

"Other" is `not_called` for every category that was assessable. It is not the same as
`not_assessable`. A protein from a clade where the adhesin modules are out of scope gets
`not_assessable` for adhesin, not "other".

The category status is the weakest status among its required modules. A `called` value with status
`unvalidated` is printed as `called (unvalidated)` in every output. The architecture rule "hold every
predictor to an acceptance test and let it fail loudly" applies: a module with an acceptance test
that fails writes a warning into the report header.

Categories are not exclusive. A protein can be an adhesin, an antigen and a surface glycoprotein.
Calls are kept as separate columns. The orchestrator does not merge them into one score
(architecture rule 4).

### 3.4 Categories in version 1

The owner named five. The cell wall publications (mapping note, section 2) add three more that are
defined by function.

| Category | Rule sketch | Required modules | Today |
|---|---|---|---|
| `surface_glycoprotein` | step 1 call | `step1_rule` or `step1_ml` | measured, undecided |
| `adhesion_repeat` (2a) | repeat detector call and `surface_glycoprotein` | repeat detectors, step 1 | validated in Saccharomycotina only |
| `adhesion_domain` (2b-i, 2b-iii, 2c) | Pfam hit (PF05730, PF04681, PF01185, PF06766, ALS families) and `surface_glycoprotein` | Pfam scan, step 1 | HMMs exist, no wrapper, no specificity test |
| `cys_rich_secreted` | Cys-rich finder call | `cys_candidates` | thresholds from 7 controls |
| `antigen_candidate` | antigen rank in top tier and `surface_glycoprotein` | antigen layer | *Coccidioides* only, not calibrated |
| `allergen_candidate` | homology to a curated allergen set, or an allergen Pfam hit | allergen module | does not exist |
| `wall_remodeling_enzyme` | family table hit (for example PF03198, PF00704) and `surface_glycoprotein` | Pfam scan, family table, step 1 | does not exist |
| `wall_synthesis_enzyme` | family table hit (PF01644, PF03142, PF02364, GPI biosynthesis) and transmembrane segments | Pfam scan, family table, TM call | does not exist |
| `other` | all assessable categories `not_called` | all above | derived |

Decision D2 (section 9) asks which of these are in version 1.

Not categories in version 1: cell wall integrity signaling (cytosolic kinases and sensors),
septation and polarized growth, polysaccharide chemistry, moonlighting proteins, biofilm. Reasons are
in the mapping note, sections 2 and 5.

### 3.5 The family table

`wall_remodeling_enzyme`, `wall_synthesis_enzyme` and the domain adhesin classes need a table that
maps Pfam accessions to a class. Rules for the table:

1. Each row has: Pfam accession and name, class, source publication with PMID, the Pfam release in
   which the accession was verified, and a `specificity_note`.
2. A row is `active` only after a specificity test: run the HMM on the proteomes of two or more
   clades with known members and known negatives, and report hits that are not members. The plan
   already requires this for CFEM, Bys1 and hydrophobin.
3. Large families (GH18, GH5, Asp, Cu-oxidase) have intracellular and non-wall members. Such a row
   needs a second condition (signal peptide, GPI call or transmembrane count) before it is `active`.
4. Sub-class assignment inside a family (chitin synthase classes I to VII) is **not** done by this
   table. A row can carry a pointer to a separate phylogenetic placement module, which is out of
   scope for version 1.

### 3.6 Interfaces and files

- Command: `cellsurface_sort --fasta P.faa --clade <name> --workdir W --out O` (name: decision D3).
- Heavy modules run as SLURM jobs with `$SCRATCH` (`${SCRATCH:?}`), never with
  `BASH_SOURCE`. The driver submits the jobs, waits, then reads their outputs. A module output that
  exists and whose input hash matches is reused.
- Output tables are written as `.tsv.gz` (portable). Large intermediates inside `$WORKDIR` may use
  `.zst`. Readers accept plain, `.gz` and `.zst`.
- `--clade` is required. It selects the `clade_scope` check. There is no default clade, so a user
  cannot get a silent out-of-scope call.
- Clade names use the labels in `analysis/step1_compare/species.tsv`.

## 4. Test panel and acceptance

A panel of proteins with expected category sets, built from the curated tables and the literature
rows that already carry a PMID. Rules:

1. Each panel protein has a source (PMID or table row) and an expected value for each category it
   is known for, plus categories it is known **not** to be (hard negatives, issue #14).
2. Thresholds are set on a development part of the panel. A held-out part is used once for the
   report. A protein used to set a threshold cannot be in the held-out part (leakage control, as in
   the step 1 spec).
3. Report recall and false-positive rate per category and per clade only where truth exists. Where
   truth does not exist, print `unvalidated`, not a number.
4. Unit tests (CI): rule evaluation, three-valued logic, scope checks, status propagation, input
   hashing and cache reuse, reader of `.gz`/`.zst`. Slow tests (HPCC): a tier harness writes the
   JSON that each module's `status_source` points to.
5. Golden-output test: a fixed small FASTA with a stored `calls.tsv.gz`.

Candidate panel members and the categories they test, to be confirmed against the curated tables:
SOWgp (adhesion_repeat, antigen), Als1 and Flo11 (adhesion_repeat), Ag2/PRA and Rbt5 (CFEM),
RodA (hydrophobin), CalA (Bys1), Gel1 (wall_remodeling_enzyme, not adhesion), a chitin synthase
and Fks1 (wall_synthesis_enzyme, not surface_glycoprotein), Hsp60 (excluded by design). Whether
these have labels in `data/curated/` is not checked here.

## 5. Statistics

- Intervals are 95% cluster bootstrap by homology group (as in Phase C).
- Per-category thresholds are fixed before the held-out part is read.
- Agreement between two modules that make the same kind of call (`step1_rule` and `step1_ml`; the
  two repeat detectors) is reported as Cohen's kappa with a bootstrap interval.
- No pooled number across clades unless each clade passes the estimate-set rule.

## 6. Compute and effort

Estimates, not measurements, except where stated.

| Item | Size |
|---|---|
| SignalP 6 plus PredGPI, 9,910 proteins | measured: 2 min 54 s on one GPU |
| Pfam scan of one proteome | not measured. The Pfam-A file is 2.2 GB. |
| ESM-2 embedding of one proteome | not stated here; see the step 1 spec and issue #16 |
| Fungi_5k (488+ proteomes) | per issue #16, one embedding run to S3. A Pfam-only scan of all proteomes is much smaller. |
| Code | driver, module interface, config reader, three wrappers (SignalP/GPI, Pfam, repeat), report: not estimated |

## 7. Risks

1. The step 1 decision changes what `surface_glycoprotein` means as a rule input. Mitigation:
   the category rule names the module (`step1_rule` or `step1_ml`), not "step 1".
2. Rules that depend on `surface_glycoprotein` inherit its false-positive rate. R2 has FPR 0.006 and
   recall 0.42. Half of the surface proteins would fail the first condition.
3. Large Pfam families give many hits that are not cell wall genes (section 3.5 rule 3).
4. Clade scope hides most of the genome. In Onygenales and Basidiomycota many categories will be
   `not_assessable`. This is correct but may look like a failure. The report must show the share.
5. Truth is thin. Allergen has none. Antigen has *Coccidioides* only.
6. Cached module output can go stale when a model or database changes. Input hash must include the
   module version and the database release (for example the Pfam release).
7. The checkout can be shared. Check `git branch --show-current` before every commit.

## 8. Order of work

1. This spec, then an independent review by a different model, then a plan.
2. Allergen scoping note: answer the open questions in issue #19 (how many fungal families, how
   many surface versus intracellular, whether homology to known allergens is the only signal).
3. Family table, first version: CFEM, Bys1, hydrophobin, Als, GEL/GAS, GH18, chitin synthase, Fks,
   GPI biosynthesis. Specificity tests.
4. Module wrappers: SignalP/GPI, Pfam, repeat detector. Then the driver and report.
5. Panel and tier harness.

## 9. Decisions for the owner

Each has a recommendation. One question at a time, as the owner prefers.

| # | Decision | Recommendation |
|---|---|---|
| D1 | Multi-label or exclusive categories | Multi-label (section 3.3). Matches issue #19 and architecture rule 4. |
| D2 | Categories in version 1 | The five the owner named, plus `wall_remodeling_enzyme`, so that GPI-anchored enzymes are not "other". Add `wall_synthesis_enzyme` in version 2. |
| D3 | Name and place of the package and command | Do not add it to `surface_glyco`: that package is step 1 only. Use a new package `src/cellsurface_sort/` and a command `cellsurface_sort`. The owner may prefer another name. |
| D4 | Does the orchestrator wait for the step 1 decision | No. Carry `step1_rule` and `step1_ml` as two modules (section 3.2). |
| D5 | Include cell wall biosynthesis and remodeling families in version 1 | Remodeling yes (see D2). Synthases version 2. Signaling no. |
| D6 | Allergen scope in version 1 | Homology to a curated set plus allergen Pfam models, flagged `unvalidated`. Needs the scoping note first. |
| D7 | Out-of-scope clade handling | `not_assessable` (section 3.3). `--clade` required. |
| D8 | Execution engine | Python driver with SLURM job scripts for version 1. Nextflow later, if the module count grows. |
| D9 | Review model for this spec | A different model from the author (Fable or Opus), as in earlier specs. |

## 10. Deliverables

1. This spec after review (revision 2).
2. `docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md` (written).
3. A plan (not written). Code only after plan review.
