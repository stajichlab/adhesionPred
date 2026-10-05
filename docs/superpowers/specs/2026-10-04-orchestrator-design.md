# Design spec: `cellsurface_sorting_hat` (the orchestrator)

*Drafted 2026-10-04. Revision 3, the same day, after independent review 1
(`2026-10-04-orchestrator-design-review-1.md`: 2 blocker, 13 major, 10 minor findings, all
dispositioned). DRAFT. Revision 3 has not been reviewed. No code, data or job exists for this spec.
Owner decisions are in section 9; D10 to D12 are new.*

Inputs: `docs/PLAN-2026-09-30-pipeline-and-decisions.md` section 6 (the proposal this spec
expands); `docs/TOOL-ARCHITECTURE.md`; `docs/model-review/STATUS.md` (2026-10-02);
`docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md`;
`docs/reports/2026-10-04-fungal-allergen-scoping.md`;
`docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` (step 1 labels);
`docs/HANDOFF-2026-10-03.md`; the code under `analysis/` named in section 2.

## 1. Purpose and non-goals

**Purpose.** One entry point. It takes a proteome and writes one row per protein. Each row says,
per tool, what the tool found, and for which clades that tool has been measured. A second layer turns
the evidence into category calls. The main goal is to find cell surface proteins for five tasks:
surface glycoprotein, cell wall and adhesion candidate, antigen, allergen, and other. Enzyme classes
(synthases, glucanases, degradation enzymes) are later work (section 3.7).

**Non-goals.**
- It does not choose the step 1 method (rule, ML or hybrid) and does not set gate values. Every call
  that depends on step 1 is computed once per step 1 variant (section 3.4), so the choice can be made
  later without changing the design.
- It does not train a model.
- It does not create truth data. Allergen and Onygenales curation are separate work.
- It does not replace any single tool. Every tool still runs alone.
- It does not predict cell wall chemistry, signaling, or expression changes. Expression appears only
  as an evidence column.

## 2. Facts checked

State on 2026-10-04. Items marked (review) were checked by the independent reviewer and not re-derived.

| Item | State | Where |
|---|---|---|
| Step 1, ESM-2 plus logistic regression | CLI exists. **No model ships**: `src/surface_glyco/models/` holds only a README. | `src/surface_glyco/` |
| Step 1 rules R0, R1, R2 | R0 has no parameter ("SignalP calls a signal peptide"). R1 and R2 have fitted parameters (GPI cut g, Ser+Thr cut t) chosen per training fold by Youden's J. No frozen values. Not exposed as a tool. | `analysis/step1_compare/phasec/rules.py` (review) |
| Step 1 calls on *C. immitis* RS (9,910 proteins) | R0 460, R2 94, M8 V-go 758, M8 V-kw 5,706. The choice changes results a lot. | `_workdir/step1_compare/phasec/proteome_calls.tsv.gz` (review) |
| Step 1 rule against ML agreement | Measured in Phase C per proteome (example: *A. fumigatus*, M8 V-go: 194 in both, 9 rule-only, 1,115 ML-only). | `_workdir/step1_compare/phasec/report.md` (review) |
| Step 1 features | SignalP 6 and PredGPI run as jobs on a protein FASTA. 9,910 proteins in 2 min 54 s on one GPU. Not wrapped. | `analysis/step1_compare/jobs/j1_features.sh`; time from `docs/HANDOFF-2026-10-03.md` |
| TM evidence | Phobius (`docs/reports/2026-10-03-cocci-spherule-followup.md`) and TMHMM 2.0c (`analysis/cocci_spherule/01b_tmhmm.sh`) were run once. Not wrapped, no validation here. R2 calls the PM-TM negatives MSB2 and HKR1. | (review) |
| Repeat detectors (2a) | Two scripts. Neither replaces the other. Measured on a synthetic divergence series plus one real family (SOWgp). **No clade truth set.** The PR-AUC 0.94 to 0.98 in Saccharomycotina belongs to the ESM classifiers, not to the detectors. | `analysis/cocci_repeats/02_repeat_profile.py`, `14_repeat_detect_general.py`; `REPORT_2026-09-29_repeat_detector_divergence.md` section 7 |
| Family HMM scans (CFEM, Bys1, hydrophobin) | A script already runs `hmmsearch --cut_ga` for PF05730, PF04681, PF01185, PF06766, PF28987 (DewD) and PF28404, and records provenance. Hydrophobin model Eas PF22354 exists in Pfam but is not in the script. No specificity test. | `analysis/cys_candidates/01_known_family_hmm.sh` |
| Pfam database | Central link `/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm`, resolved with release recorded in `pfam_provenance.json` (release `2026-01-27-Pfam38.2` in the last run). The funannotate copy has no release file. | `analysis/cys_candidates/README.md` |
| Cys-rich finder | No model, no training, "no accuracy claim". Thresholds from 7 controls. Calibrated on *Coccidioides* proteomes. Reads SignalP output. | `analysis/cys_candidates/` |
| Antigen layer | *Coccidioides* only. 3 of 4 anchors pass the acceptance test (PRA3, Ag2/PRA, SOWgp pass; PRA2 fails). The script prints `NOT CALIBRATED` because one fails. Reads the Fungi_5k DuckDB (IDs `FA2214EC`), an ID map to `CIMG_` IDs, RS1 kallisto TPM and prevalence over 488 isolates. SOWgp (prevalence 0.9201) is in neither Tier 1 nor Tier 2. | `analysis/cocci_antigens/` |
| Expression layer | One data set (*C. immitis* spherule against mycelium, 2 replicates). 21 spherule genes are absent from the Fungi_5k RS annotation. | `analysis/cocci_spherule/`, `docs/HANDOFF-2026-10-03.md` |
| Allergen layer | **No module.** Scoping done: WHO/IUIS lists 120 fungal allergen molecules (4 Onygenales, none from *Coccidioides*). About one third have a signal peptide. Extracts exist. | `docs/reports/2026-10-04-fungal-allergen-scoping.md`; `analysis/allergen_scoping/` |
| Biofilm layer | Does not exist. Blocked on phenotype data. | `docs/TOOL-ARCHITECTURE.md` |
| Orchestrator | Does not exist. | plan section 6 |
| Curated tables | All drafts. Rows with `needs_review=yes` (review): adhesins 224 of 268, antigens 86 of 86, biofilm 207 of 207. `surface.tsv` has about 3,573 rows and no `needs_review` column. `antigens/coccidioides_candidates.tsv` has 1,069 rows. Hard-negative seeds: 31 rows. | `data/curated/` |
| Step 1 label classes | P-ext, P-gpi, N-int, N-sec, PM-TM, ambiguous. | step 1 spec 2.2 |

## 3. Design

### 3.1 Layers

```
proteome FASTA  (+ optional per-protein taxon map)
   |
   v
[ evidence modules ]  cached; one run state per module and per protein
   L  localisation      SignalP 6, GPI call, TM call (evidence), step 1 variants
   D  family / domain   hmmsearch against Pfam-A, curated family table
   A  architecture      repeat detectors, composition, Cys-rich finder
   H  homology          allergen set (WHO/IUIS), later AllergenOnline and COMPARE
   K  lookup by ID      antigen ranking, expression (precomputed tables)
   |
   v
evidence.tsv.gz         long format: protein, module, variant, field, value
   |
   v
[ category rules ]      config file, Kleene three-valued logic
   |
   v
calls.tsv.gz            one row per protein; per category: value and _status columns
report.md               counts, run states, scope, the share not assessable
```

### 3.2 Evidence modules

Each module is a plugin with this interface.

| Field | Meaning |
|---|---|
| `name`, `version`, `kind` | for example `signalp6`, `6.0`, L |
| `params_hash` | hash of all parameters (cutoffs, SignalP mode, Cys-rich L/K/W, R2 g and t) |
| `artefact_hash` | hash of the database or model it uses (Pfam release file, model weights, ranking table) |
| `scope(taxon)` | the taxa for which the module has been measured. Free text plus NCBI taxon IDs. |
| `status(module, version, taxon)` | one of `estimated`, `smoke`, `unvalidated` (below) |
| `status_source` | path to the measurement JSON that supports the status. Its recorded module version must equal `version`, else the status is `unvalidated`. |
| `run(inputs, workdir)` | writes a per-protein table and a run record |
| `fields` | names and types of the output columns |

**Status** is a function of module, version and clade, because one module can be an estimate in one
clade and a smoke test in another (step 1: estimate in S1 and Eurotiomycetes, smoke test in
Basidiomycota). Values: `estimated` (recall interval half-width at most 0.10 and at least 20 direct
positives: the Phase C rule; the word replaces `validated`, which reads as "good"), `smoke` (fewer),
`unvalidated` (no truth). **In scope** means status is `estimated` or `smoke` for the protein's taxon.
A module run outside its scope writes values and `in_scope = false`.

**Taxa.** Scope uses NCBI taxon IDs and lineage matching (*Coccidioides* lies inside Eurotiomycetes).
`--taxon` sets one taxon for the whole run. A per-protein taxon map is also accepted, which a
multi-clade test panel needs. A taxon given by the user is recorded in the report and is not
checked against the sequences.

**Run states** (per module, per run): `ok`, `unavailable` (no frozen artefact, for example step 1 ML
with no model), `not_run`, `error`. A protein-level value can also be `na_too_short`, `na_window`,
`na_invalid` (see section 3.8).

**Rules.**
1. A module cannot set its own status. Without a `status_source` it is `unvalidated`.
2. Step 1 appears as separate **variants**: `step1_rule@R0`, `step1_rule@R1`, `step1_rule@R2` and
   `step1_ml@<card>`. R1 and R2 are `unavailable` until the owner freezes g and t. An ML variant is
   `unavailable` until a model card exists. Agreement between variants is cited from the Phase C
   tables and re-measured only if a definition changes.
3. Detectors that make the same kind of call (the two repeat detectors) run as separate modules.

**Kind K (lookup by ID).** Antigen and expression read precomputed tables keyed on database IDs, not
on the user's sequence. Such a module needs an ID-mapping step (sequence hash first, then mmseqs)
with a reported match rate. A protein without a match gets `not_in_reference`. This is a run state
of the protein for that module and makes dependent categories `not_assessable`.

### 3.3 Category logic

Values: `called` (T), `not_called` (F), `not_assessable` (U, unknown). Kleene three-valued logic.

| AND | T | F | U |
|---|---|---|---|
| T | T | F | U |
| F | F | F | F |
| U | U | F | U |

| OR | T | F | U |
|---|---|---|---|
| T | T | T | T |
| F | T | F | U |
| U | T | U | U |

NOT T = F, NOT F = T, NOT U = U. A module that is out of scope, `unavailable`, `error` or
`not_in_reference` for a protein gives U for the inputs it provides. A rule is evaluated with these
tables. So a positive call is kept when another input is U (T OR U = T), and a false input makes an
AND false whatever the other input is.

**Status columns.** The value column holds only `called`, `not_called`, `not_assessable`. The
category status is a separate column `<category>_status`: the weakest status among the modules that
set the result (`estimated` > `smoke` > `unvalidated`). For a `called` result the contributing
modules are the ones that were true. For a `not_called` result they are all required modules. Reports
print `called` together with its status in two columns. They never merge the two into one string.

### 3.4 Categories in version 1

Calls that depend on step 1 are written once per step 1 variant, as `<call>[<variant>]`, and also
**ungated** (without the step 1 condition). The owner chooses the default gate (decision D10).

| Category / call | Rule | Inputs | Today |
|---|---|---|---|
| `surface_glycoprotein[v]` | step 1 variant `v` calls the protein | step 1 variant | R0 measured; R1, R2, ML not frozen |
| `adhesion_repeat` (ungated) | a repeat detector calls | `repeat02` OR `repeat14` | `unvalidated` in every clade |
| `adhesion_domain` (ungated) | Pfam hit in the adhesion table (PF05730, PF04681, PF01185, PF06766, PF28987, PF22354, ALS families) | Pfam scan | HMMs run once; no specificity test |
| `cell_wall_adhesion_candidate[v]` | (`adhesion_repeat` OR `adhesion_domain`) AND `surface_glycoprotein[v]`; ungated form without the AND | the above | see rows |
| `antigen_candidate` (ungated) | antigen ranking `percentile_dedup` at most P | antigen lookup | *Coccidioides* only; P default 10 (the "top decile" of the acceptance test) |
| `antigen_candidate_surface[v]` | `antigen_candidate` AND `surface_glycoprotein[v]` | antigen lookup, step 1 | see above |
| `allergen_candidate` | best hit to the WHO/IUIS fungal allergen set with identity and coverage at or above the report cutoffs (default 35% over 80 aa, the FAO/WHO rule), or an allergen-specific Pfam hit (PF16541, PF25312). **No surface gate**: most fungal allergens are intracellular. | allergen homology module, Pfam scan | `unvalidated`; module to be written |
| `other_not_surface[v]` | `surface_glycoprotein[v]` is F, and no other category is T | all above | derived |
| `other_surface_no_mechanism[v]` | `surface_glycoprotein[v]` is T, and every mechanism category (`adhesion_repeat`, `adhesion_domain`, `antigen_candidate`, `allergen_candidate`) is F | all above | derived |

`other_*` is U whenever one of its inputs is U and none is T. So a protein is "other" only when every
category that can be assessed is F and none is U. The report prints how many proteins are U for each
reason. The two `other` values keep the two meanings that the plan requires.

Evidence columns that are **not** categories in version 1:
- `cys_rich_sp_unassigned` (Cys-rich finder, tier named explicitly). No functional claim: the finder's
  README says "Cys-rich does not mean PRA3-like". Scope: *Coccidioides*. Inputs: SignalP, Pfam scan.
  `cys_rich_sp_known_family` proteins (which include CFEM proteins) are reported as such and not
  counted here.
- TM evidence (Phobius, TMHMM). Shown beside `surface_glycoprotein` so that a reader can see PM-TM
  proteins such as MSB2 and HKR1.
- Expression (spherule against mycelium). Evidence only.

Known limits, printed in the report header:
1. GPI-anchored and secreted enzymes (GEL/GAS, chitinases, glucanases, proteases) get only
   `surface_glycoprotein`. Non-adhesive structural wall proteins (Cwp1, Ccw12, Sed1, Pir) get the same.
   Decision D11 asks whether to add a `cell_wall_protein` call.
2. `surface_glycoprotein` is defined by GO cell wall and extracellular region evidence. It is not
   evidence of glycosylation.
3. CFEM is filed under adhesion because the 2b-i class is, but its confirmed fold is a hemophore
   (`docs/TOOL-ARCHITECTURE.md`). Binding to a host receptor is not shown.
4. The repeat detectors have no clade truth set. Their calls are hypotheses.
5. Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, moonlighting
   proteins and biofilm are not categories.

### 3.5 The family table

The `adhesion_domain` call (and, later, the enzyme classes of section 3.7) needs a table that maps
Pfam accessions to a class. Rules:

1. Each row has: Pfam accession and name, class, source publication with PMID, the Pfam release in
   which the accession was verified, and a `specificity_note`.
2. A row is `active` only after a specificity test: run the HMM on the proteomes of two or more
   clades with known members and known negatives, and report hits that are not members.
3. Large families (GH18, GH5, Asp, Cu-oxidase) have intracellular and non-wall members. Such a row
   needs a second condition (signal peptide, GPI call) before it is `active`.
4. Sub-class assignment inside a family (chitin synthase classes I to VII) is not done by this
   table.
5. Start from `analysis/cys_candidates/01_known_family_hmm.sh`. Use the central Pfam path and record
   the release (section 3.6).

### 3.6 Provenance, caching and files

- **Cache key** for a module output: sha256 of the input, module `version`, `params_hash`, and
  `artefact_hash` (Pfam release file, model weights, ranking table). The hash of `categories.yaml`
  is written into `calls.tsv.gz`. A `status_source` whose module version differs from the running
  version is refused.
- **Writes are atomic** (write to a temporary name, then rename), with a sha256 sidecar, as
  `j1_features.sh` does. Two runs on the same workdir do not share a partial file.
- **Protein key** is the sequence sha256, with the user's ID kept as a label (Phase C uses
  `seq_sha256`). Duplicate IDs are an error.
- **Pfam** comes from the central link, resolved, with the release and sha256 recorded. The
  funannotate copy is not used.
- Command: `cellsurface_sorting_hat --fasta P.faa --taxon <NCBI id> --workdir W --out O`
  (also `--taxon-map FILE`). Package `src/cellsurface_sorting_hat/`.
- Heavy modules run as SLURM jobs with `$SCRATCH` (`${SCRATCH:?}`), never with `BASH_SOURCE`.
  The driver submits jobs and reads their outputs. It sets a timeout and treats a killed or
  preempted job as `error` for that module. It does not trust queue-time estimates.
- Tables are `.tsv.gz`. Large intermediates may use `.zst`. Readers accept plain, `.gz`, `.zst`.

### 3.7 Later work: enzyme classes (not version 1)

The owner set these aside on 2026-10-04. No change to the module interface or the logic is needed.

| Later category | Rule sketch | Extra need |
|---|---|---|
| `wall_remodeling_enzyme` | family table hit (PF03198 GEL/GAS, PF00704 GH18) AND `surface_glycoprotein[v]` | family table with specificity tests |
| `wall_degradation_enzyme` | family table hit; large families need a second condition (3.5 rule 3) | same |
| `wall_synthesis_enzyme` | family table hit (PF01644, PF03142, PF02364, GPI biosynthesis) AND at least one TM segment | a validated TM module |

Mapping of the cell wall publications to these classes:
`docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md`.

### 3.8 Failure modes

| Case | Per-protein value | Report |
|---|---|---|
| Module unavailable (no frozen artefact) | U for its inputs; run state `unavailable` | listed in the header with the reason |
| Module job fails, times out, is preempted | run state `error` | header warning; dependent categories U |
| Protein out of the module's scope | value written, `in_scope = false` | share per module |
| PredGPI `too_short` | `na_too_short`; GPI input U | count |
| Protein longer than the ESM window (1,022 aa) | `na_window` unless a window rule exists (issue #10) | count |
| Internal `*`, invalid residues, empty sequence | `na_invalid`; protein excluded from all modules | count, list |
| FASTA header differs from the SignalP ID (Cys-rich finder needs an exact match) | `error` for that module | count |
| Duplicate IDs | run refused | error |
| No ID match in a lookup table | `not_in_reference` | match rate |

## 4. Test panel and acceptance

**Acceptance for version 1 is software correctness**, tested on stored module outputs (fixtures), so
CI needs no GPU and no SignalP:
1. Unit tests: the Kleene tables (all 27 AND and OR cases), `other_*` rules, status propagation,
   scope checks with lineage matching, cache key and refusal of a stale `status_source`, atomic
   writes, readers of `.gz` and `.zst`, every row of the failure-mode table.
2. Golden test: a small fixed set of module-output fixtures gives a stored `calls.tsv.gz`.
3. A slow tier harness on HPCC runs the real modules and writes the JSON that each `status_source`
   points to. It is not part of CI.

**Per-category accuracy is report-only** in version 1. Recall and false-positive rate are printed per
category and per clade only where truth exists. Elsewhere the report prints `unvalidated`. The
`estimated` rule (20 direct positives) cannot be met by the adhesion sub-calls (Bys1 has 1 positive,
PRA3 1), the antigen layer (4 anchors) or the allergen module (no *Coccidioides* truth). That is
expected, and the report says so.

**Thresholds the orchestrator owns:** the antigen percentile P and the allergen identity and coverage
cutoffs. Both are config items with the defaults above and are reported with every run. All other
cutoffs belong to the modules and come in through `params_hash`.

**Panel** (a description of expected calls, for the report; not an acceptance gate):
- Each protein has a source (PMID or table row) and an expected value per category: `called`,
  `not_called`, `known_miss`, or `excluded`.
- Record for each (module, protein) pair whether the module, or a homolog of the protein, was seen
  during tuning or training. Exclude those pairs from the metrics. Controls that tuned modules:
  PRA3, Ag2/PRA, PRA2, RodA, CalA, SOWgp, CTS1 (Cys-rich finder); SOWgp, PRA3, Ag2/PRA, PRA2, CF
  antigen (antigen acceptance test); SOWgp alleles (repeat benchmark); Als1 and Flo11 (step 1 S1
  training or out-of-fold). Homologs count (calB and calC with CalA; the Als family).
- Split development and held-out parts by homology cluster, not by protein.
- Candidates: SOWgp (`adhesion_repeat`, `antigen_candidate`), Als1 and Flo11 (`adhesion_repeat`),
  Ag2/PRA and Rbt5 (`adhesion_domain`), RodA (`adhesion_domain`), CalA (`adhesion_domain`), Gel1
  (`surface_glycoprotein` only), Asp f 1, Asp f 2, Asp f 34 (`allergen_candidate`), Hsp60
  (`known_miss`: moonlighting, no signal peptide; it is a known antigen and adhesin, so `not_called`
  would encode a false fact).
- The hard-negative source is the 31-row seed list. Issue #14 (curated hard negatives) is open.
- Whether these proteins have labels in `data/curated/` is not checked here.

## 5. Statistics

- Intervals are 95% cluster bootstrap by homology group (as in Phase C).
- Agreement between modules of the same kind is reported as Cohen's kappa with a bootstrap interval
  **and** as the 2x2 counts. Kappa is unstable at low prevalence (R2 calls 94 of 9,910 RS proteins).
- No pooled number across clades unless each clade passes the `estimated` rule.

## 6. Compute and effort

Estimates, not measurements, except where marked.

| Item | Size |
|---|---|
| SignalP 6 plus PredGPI, 9,910 proteins | measured: 2 min 54 s on one GPU (`docs/HANDOFF-2026-10-03.md`) |
| Pfam scan of one proteome | not measured. A *Coccidioides* run took about 25 s on a login node for six models (`analysis/cys_candidates/README.md`); a full Pfam-A scan is much larger and was not timed. |
| ESM embedding of one proteome | not measured here |
| Fungi_5k embedding (issue #16) | about 58 M proteins; 60 to 150 GPU-hours for ESM C 300M (estimate in the issue) |
| Code | driver, module interface, config reader, wrappers (SignalP/GPI, Pfam, repeat, allergen homology, lookup), report: not estimated |

(488 is the number of *Coccidioides* isolates, not the number of Fungi_5k proteomes.)

## 7. Risks

1. The step 1 decision changes what `surface_glycoprotein` means. Mitigation: every gated call is
   written per variant, and an ungated call is always written.
2. The rule gates lose recall. R2 recall is 0.418 (S1), 0.227 (Eurotiomycetes), 0.125 (Basidiomycota,
   smoke test), so R2 misses 58%, 77% and 88% of surface proteins there. On the *Coccidioides*
   antigen candidates R2 keeps 5 of 14 Tier 1 and 19 of 45 Tier 2 (review); R0 keeps all. The main
   cost of a gate is lost recall, not false positives. R0 has higher recall (0.603, 0.727, 0.938) and a higher
   false-positive rate (0.037, 0.010, 0.083).
3. Large Pfam families give hits that are not cell wall genes (section 3.5 rule 3).
4. Scope hides most of the genome. In Onygenales and Basidiomycota many categories will be
   `not_assessable`. The report must show the share and the reason.
5. Truth is thin. Allergen has no *Coccidioides* truth. Antigen has *Coccidioides* only.
6. Cached output can go stale. The cache key covers version, parameters and artefact hashes.
7. The checkout can be shared. Check `git branch --show-current` before every commit.

## 8. Order of work

1. This spec (revision 3): a second independent review if the owner wants one, then a plan.
2. Allergen: COMPARE and AllergenOnline downloads, SignalP/PredGPI on the IUIS sequences, the
   negative set for validation (`docs/reports/2026-10-04-fungal-allergen-scoping.md`).
3. Family table, first version: CFEM, Bys1, hydrophobin, Als. Specificity tests.
4. Module wrappers: SignalP/GPI, Pfam, repeat detector, allergen homology, lookup modules with ID
   mapping. Then the driver, the category engine and the report.
5. Fixtures, unit tests, golden test.
6. Later (section 3.7): enzyme classes, extra family table rows, a validated TM module.

## 9. Decisions for the owner

One question at a time.

| # | Decision | Status / recommendation |
|---|---|---|
| D1 | Multi-label or exclusive categories | Multi-label. Recommended; not yet answered. |
| D2 | Categories in version 1 | **Answered 2026-10-04:** the five named. Enzyme classes later. |
| D3 | Name and place | **Answered 2026-10-04:** `cellsurface_sorting_hat`, package `src/cellsurface_sorting_hat/`. |
| D4 | Wait for the step 1 decision | No. Variants are carried (section 3.2). |
| D5 | Biosynthesis and remodeling families in version 1 | **Answered by D2:** later. Signaling stays out. |
| D6 | Allergen scope in version 1 | Homology to the WHO/IUIS fungal set plus allergen-specific Pfam; no surface gate; `unvalidated`. |
| D7 | Out-of-scope handling | `not_assessable` with `--taxon` required. |
| D8 | Execution engine | Python driver with SLURM scripts; Nextflow later. |
| D9 | Review model | Different model from the author. Review 1 done (Opus). |
| **D10** | **Default gate for the gated calls** (`surface_glycoprotein[v]` used by `cell_wall_adhesion_candidate[v]` and `antigen_candidate_surface[v]`) | Recommend `step1_rule@R0` (no parameter, recall 0.603 / 0.727 / 0.938), with R2 and ML as extra columns once frozen. Ungated calls are always written. |
| **D11** | Add a `cell_wall_protein` call for non-adhesive structural wall proteins (Cwp1, Ccw12, Sed1, Pir)? | Not in version 1. It needs a GPI call and a curated family list that does not exist (issue #14). |
| **D12** | Antigen call definition | Percentile cut P of the ranking, default 10. The tiers are not used because they exclude SOWgp. |

## 10. Deliverables

1. This spec (revision 3) and `2026-10-04-orchestrator-design-review-1.md`.
2. `docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md` and
   `docs/reports/2026-10-04-fungal-allergen-scoping.md` (written).
3. A plan (not written). Code only after plan review.
