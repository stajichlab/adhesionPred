# Design spec: `cellsurface_sorting_hat` (the orchestrator)

*Drafted 2026-10-04. Revision 5, the same day, after independent review 2 (3 blocker, 9 major, 10 minor
findings; see `2026-10-04-orchestrator-design-review-2.md`). Revision 3 answered review 1
(`2026-10-04-orchestrator-design-review-1.md`). DRAFT. Revision 5 has not been reviewed. No code, data or
job exists for this spec. Owner decisions D1 to D13 are answered (section 9); D6 and D7 were changed
by review 2.*

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
| Allergen layer | **No module.** Scoping done: WHO/IUIS lists 120 fungal allergen molecules (4 Onygenales, none from *Coccidioides*). 30 of the 97 IUIS accessions found in UniProt (31%) have a signal peptide feature. Extracts exist. | `docs/reports/2026-10-04-fungal-allergen-scoping.md`; `analysis/allergen_scoping/` |
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
evidence.tsv.gz         long format: protein, module, field, value (fields listed in categories.yaml)
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
| `applicable(taxon, input)` | whether the module can run on this input and taxon (for example the antigen lookup: *Coccidioides* only; Pfam scan: any fungus). Not a measure of quality. |
| `status(module, version, taxon)` | one of `estimated`, `smoke`, `unvalidated` (below). Measured on the taxa listed in `status_source`. |
| `status_source` | path to the measurement JSON that supports the status. It records module `version`, `params_hash`, `artefact_hash` and the list of tested taxa. If any of the three differs from the running module, the status is `unvalidated`. |
| `run(inputs, workdir)` | writes a per-protein table and a run record |
| `fields` | names and types of the output columns |

**Applicability and status are two separate things** (review 2, finding 1).
- *Applicability* says whether the module can produce a value for this protein and taxon. If it
  cannot, the module gives U for its inputs. Examples: antigen and Cys-rich finder outside
  *Coccidioides*; a step 1 ML variant with no model.
- *Status* says how well the module has been measured. An applicable module with status
  `unvalidated` still gives `called` or `not_called`. The status goes in the `_status` column.

**Status** is a function of module, version and taxon, because one module can be an estimate in one
clade and a smoke test in another (step 1: estimate in S1 and Eurotiomycetes, smoke test in
Basidiomycota). Values: `estimated` (recall interval half-width at most 0.10 and at least 20 direct
positives: the Phase C rule; the word replaces `validated`, which reads as "good"), `smoke` (fewer),
`unvalidated` (no truth).

**Taxa.** The orchestrator reads an NCBI taxonomy dump (the version is recorded in the report).
A status applies to a protein's taxon only if that taxon is the tested taxon or a descendant of a
tested taxon listed in the `status_source`. It does **not** apply to a taxon only because both share a
broad label such as "Eurotiomycetes". Example: the step 1 estimate was measured on *A. fumigatus* and
*A. nidulans*, so it covers those species and their descendants; *Coccidioides* (Onygenales) is not
covered and gets `unvalidated` (or `smoke` where a smoke test lists it). When two entries match, the
most specific one wins. The status column also records `status_basis` (the tested taxon matched).
`--taxon` sets one taxon for the whole run. `--taxon-map FILE` sets a taxon per protein and overrides
`--taxon` for the proteins it lists. At least one of the two is required; a protein with no taxon is
an error. A taxon the user gives is recorded in the report and is not checked against the sequences.

**Run states** (per module, per run): `ok`, `partial` (output for some proteins only; the missing
proteins get `error`), `unavailable` (no frozen artefact, for example step 1 ML with no model),
`not_run`, `error`. A protein-level value can also be `na_too_short`, `na_window`,
`na_invalid` (see section 3.8).

**Rules.**
1. A module cannot set its own status. Without a `status_source` it is `unvalidated`.
2. Step 1 appears as separate **variants**: `step1_rule@R0`, `step1_rule@R1`, `step1_rule@R2` and
   `step1_ml@<card>`. R1 and R2 are `unavailable` until the owner freezes g and t. An ML variant is
   `unavailable` until a model card exists. Agreement between variants is cited from the Phase C
   tables and re-measured only if a definition changes.
3. Detectors that make the same kind of call (the two repeat detectors) run as separate modules.

**Kind K (lookup by ID).** Antigen and expression read precomputed tables keyed on database IDs, not
on the user's sequence. Such a module needs an ID-mapping step:
1. exact sequence sha256 match, then ID match through a stored map (for RS, `protein_map.tsv`; the
   RS ranking IDs `CIMG_*` equal the RefSeq GCF_000149335.2 annotation per `docs/HANDOFF-2026-10-03.md`);
2. a fuzzy match (diamond or mmseqs) is written as evidence (`idmap_method`, identity, mutual
   coverage) and is used for a call only at identity >= 95% and mutual coverage >= 90%
   (proposed defaults, config items, to be fixed in the plan);
3. if two ranking rows match one protein, the better identity wins and the tie is reported.
A protein without a usable match gets `not_in_reference`. The report prints the match rate. A
protein from a taxon where the table does not apply gets U with the reason "not applicable"
(*Coccidioides* only for antigen), not `not_in_reference`.

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

NOT T = F, NOT F = T, NOT U = U. A module that is not applicable, `unavailable`, `error` or
`not_in_reference` for a protein gives U for the inputs it provides. A module that is applicable but
`unvalidated` gives a value. A rule is evaluated with these
tables. So a positive call is kept when another input is U (T OR U = T), and a false input makes an
AND false whatever the other input is.

**Status columns.** The value column holds only `called`, `not_called`, `not_assessable`. The
category status is a separate column `<category>_status`: the weakest status among the modules that
set the result (`estimated` > `smoke` > `unvalidated`). The modules that set the result are the
ones that decided it: for a true OR, the true inputs; for a false AND, the false inputs; for a true
AND or a false OR, all inputs. A `not_assessable` result has no deciding module and gets status
`unvalidated` with an empty basis. Reports print the value and its status in two columns. They
never merge the two into one string.

### 3.4 Categories in version 1

Calls that depend on step 1 are written once per step 1 variant, as `<call>[<variant>]`, and also
**ungated** (without the step 1 condition). The default gate is `step1_rule@R0` (decision D10, answered 2026-10-04). It is a config item.

| Category / call | Rule | Inputs | Today |
|---|---|---|---|
| `surface_glycoprotein[v]` | step 1 variant `v` calls the protein | step 1 variant | R0 measured; R1, R2, ML not frozen |
| `adhesion_repeat` (ungated) | a repeat detector calls | `repeat02` OR `repeat14` | applicable to any fungus; `unvalidated` in every clade |
| `adhesion_domain` (ungated) | Pfam hit in the adhesion table (PF05730, PF04681, PF01185, PF06766, PF28987, PF22354, ALS families) | Pfam scan | HMMs run once; no specificity test |
| `cell_wall_adhesion_candidate[v]` | (`adhesion_repeat` OR `adhesion_domain`) AND `surface_glycoprotein[v]`; ungated form without the AND | the above | see rows |
| `antigen_candidate` (ungated) | antigen ranking combined `percentile` (all 9,139 proteins; `percentile_dedup` is empty for 597 non-representatives) at most P, P = 15. The columns antigenicity, specificity, prevalence and max cross-reaction identity are always written beside it. | antigen lookup | *Coccidioides* only |
| `antigen_candidate_surface[v]` | `antigen_candidate` AND `surface_glycoprotein[v]` | antigen lookup, step 1 | see above |
| `allergen_homolog_hit` (evidence, not a category) | best hit to the WHO/IUIS fungal allergen set at identity >= 35% over an aligned length >= 80 aa (the FAO/WHO rule). Written for every protein with such a hit: allergen name, identity, aligned length, coverage of the allergen, aligner. | allergen homology module | applicable to any fungus; `unvalidated` |
| `allergen_candidate` | `allergen_homolog_hit` at identity >= 70% and coverage >= 80% of the allergen length (starting values, config items, chosen without a fungal non-allergen set), **or** an allergen-specific Pfam hit (PF16541, PF25312). **No surface gate**: most fungal allergens are intracellular. | allergen homology module, Pfam scan | `unvalidated`; module to be written. Decision D6 (revised). |
| `other_not_surface[v]` | `surface_glycoprotein[v]` is F, and every mechanism call that is not U is F | `surface_glycoprotein[v]`, mechanism calls | derived |
| `other_surface_no_mechanism[v]` | `surface_glycoprotein[v]` is T, and every mechanism call that is not U is F | `surface_glycoprotein[v]`, mechanism calls | derived |

Mechanism calls are `adhesion_repeat`, `adhesion_domain`, `antigen_candidate` and `allergen_candidate`
(the gated forms are implied by them; do not add them to the formulas).

**In `categories.yaml`, the ungated form of `cell_wall_adhesion_candidate` is named `cell_wall_adhesion_ungated`, and `allergen_homolog_hit` is an ungated call (identity >= 35% and aligned length >= 80 aa) that is not a mechanism category.

`other_*` is defined over the assessable categories** (decision, review 2 finding 2). For one
protein, the mechanism calls that are U are left out of the test, and the column `other_basis` lists
them (for example `antigen`). `other_*` is U only when `surface_glycoprotein[v]` is U. So on
*A. fumigatus*, `other_not_surface` means "not surface, and none of adhesion, allergen was called;
antigen was not assessed", and `other_basis` says so. The report header prints the basis counts.
A protein can be T for `other_*` in one step 1 variant and U in another. The two `other` values
keep the two meanings the plan requires.

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
   Decision D11 (answered): no `cell_wall_protein` call in version 1.
2. `surface_glycoprotein` is defined by GO cell wall and extracellular region evidence. It is not
   evidence of glycosylation. With the default gate R0, the call means "SignalP calls a signal
   peptide" and nothing more (`analysis/step1_compare/phasec/rules.py`). Secreted enzymes and ER
   proteins with a signal peptide carry the name. The gated adhesion call is therefore close to the
   ungated one for proteins with a signal peptide.
3. CFEM is filed under adhesion because the 2b-i class is, but its confirmed fold is a hemophore
   (`docs/TOOL-ARCHITECTURE.md`). Binding to a host receptor is not shown.
4. The repeat detectors have no clade truth set. Their calls are hypotheses.
5. The antigen call is the top 15% of a fixed *Coccidioides* ranking (1,371 of 9,139 proteins). It is a
   weak label. The report prints the share called and carries the ranking's `NOT CALIBRATED` note
   (3 of 4 anchors pass the top-decile test; with P = 15 all four anchors are inside the cut, so
   they cannot test it).
6. Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, moonlighting
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

- **Cache key** for a module output: sha256 of each protein sequence, plus module `version`,
  `params_hash`, `artefact_hash`, and tool versions (SignalP mode, PredGPI, hmmer, torch, aligner).
  Results are cached per sequence sha256 and expanded to IDs at the end, so adding one protein does
  not invalidate the rest, and identical sequences with different IDs are computed once. For an ML
  step 1 variant, batch size and order are part of `params_hash`. The hash of `categories.yaml`
  is written into `calls`. A `status_source` is refused if its recorded `version`, `params_hash` or
  `artefact_hash` differs from the running module.
- **Writes are atomic** (write to a temporary name, then rename), with a sha256 sidecar, as
  `j1_features.sh` does. Two runs on the same workdir do not share a partial file.
- **Protein key** is the pair (ID, sequence sha256). Calls are per ID. Identical sequences with
  different IDs are allowed and get identical calls. Duplicate IDs are an error.
- **Pfam** comes from the central link, resolved, with the release and sha256 recorded. The
  funannotate copy is not used.
- Command: `cellsurface_sorting_hat --fasta P.faa --taxon <NCBI id> --workdir W --out O`
  (also `--taxon-map FILE`). Package `src/cellsurface_sorting_hat/`.
- Heavy modules run as SLURM jobs with `$SCRATCH` (`${SCRATCH:?}`), never with `BASH_SOURCE`.
  The driver submits jobs and reads their outputs. It sets a timeout and treats a killed or
  preempted job as `error` for that module. It does not trust queue-time estimates.
- **Output schema.** `calls.long.tsv.gz` (protein, call, variant, value, status, status_basis, other_basis) is
  the primary file. A wide `calls.wide.tsv.gz` is derived from it. `evidence.tsv.gz` carries the configured evidence fields (antigen axes, allergen hit fields, Cys-rich and expression evidence). `proteins.tsv.gz` lists ID, sequence sha256, taxon, state and note. `run.json` and `report.md` record the config and taxonomy hashes. Only available step 1 variants get
  columns; `unavailable` variants are listed in the report header with the reason, not written as
  all-U columns. The schema (call names, variants) is fixed by `categories.yaml`.
- **Waiting and resources** (to be fixed in the plan): poll interval, resubmission on preemption (a
  fixed maximum), resources per module (CPU, GPU, memory, time, partition), where SignalP weights
  live, and `$SCRATCH` copy-back to `/bigdata` before a job ends.
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
| Module not applicable to the protein or taxon | U for its inputs; no value written | share per module and reason |
| PredGPI `too_short` | `na_too_short`; GPI input U | count |
| Protein longer than the ESM window (1,022 aa) | `na_window` unless a window rule exists (issue #10) | count per step 1 variant (long cell wall proteins are affected) |
| Internal `*`, characters that are not residues, empty sequence | `na_invalid`; protein excluded from all modules | count, list |
| FASTA header differs from the SignalP ID (Cys-rich finder needs an exact match) | `error` for that module | count |
| Duplicate IDs | run refused | error |
| Empty FASTA or no valid proteins | run refused | error |
| Trailing `*` | stripped silently, counted | count |
| `X`, `B`, `Z`, `U`, `J` residues | allowed up to a module-specific fraction; above it `na_invalid` for that module only | count per module |
| Module finishes for some proteins only | run state `partial`; the missing proteins are `error` | header warning |
| GPU out of memory | one retry with half the batch size, then `error` | header warning |
| User `--taxon` does not match the proteome | not detected; recorded as given | stated in the report |
| No ID match in a lookup table | `not_in_reference` | match rate |

## 4. Test panel and acceptance

**Acceptance for version 1 is software correctness**, tested on stored module outputs (fixtures), so
CI needs no GPU and no SignalP:
1. Unit tests: the Kleene tables (9 AND cases, 9 OR cases, 3 NOT cases), the `other_*` formulas and `other_basis` (including a T in one variant and U in another), status propagation,
   scope checks with lineage matching, cache key and refusal of a stale `status_source`, atomic
   writes, readers of `.gz` and `.zst`, every row of the failure-mode table.
2. Golden test: a small fixed set of module-output fixtures gives a stored `calls.tsv.gz`.
3. Lineage matching: tested taxon, descendant, sibling clade with the same broad label (no match), most specific wins, `--taxon-map` overriding `--taxon`, protein with no taxon.
4. Identical sequences with different IDs; per-sha256 caching and expansion.
5. A slow tier harness on HPCC runs the real modules and writes the JSON that each `status_source`
   points to. It is not part of CI.

**Run-level check (not CI):** one real proteome (*A. fumigatus* Af293) end to end on HPCC. It must
finish, write the report, give correct run states, and record the wall time and the resources.

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

### 4.1 First target proteomes (decision D13)

| Proteome | Source | Why | State |
|---|---|---|---|
| *C. immitis* RS | RefSeq GCF_000149335.2 protein FASTA, `_workdir/cocci_spherule/ref/GCF_000149335.2_ASM14933v2_protein.faa.gz`, **9,910 proteins** (counted 2026-10-04). Not the Fungi_5k file `Coccidioides_immitis_RS.proteins.fa` (7,630 sequences, a different annotation). | the only proteome where antigen, expression, Cys-rich and repeat results all exist; SOWgp | available |
| *A. fumigatus* Af293 | Fungi_5k file `Aspergillus_fumigatus_Af293.proteins.fa`, 9,161 proteins, IDs `F85C5601_...`; the step 1 truth set uses UniProt UP000002530 (`_workdir/step1_compare/downloads/UP000002530.fasta.gz`), different IDs | WHO/IUIS *Aspergillus* allergens (38 entries), CalA, CspA, RodA; Eurotiomycetes step 1 truth | both files available; the two are joined by sequence sha256, so only identical sequences join |
| *S. cerevisiae* S288C | SGD `orf_trans_all.fasta.gz`, `_workdir/step1_compare/downloads/orf_trans_all.fasta.gz` (no S288C file with that name was found in Fungi_5k/input); protein count not checked here | step 1 estimate with direct positives; FLO11, AGA1 | available |
| *A. fumigatus* A1163 (CEA10, FGSC A1163, CBS 144.89) | UniProt UP000001699, 9,942 proteins (found 2026-10-04) | second strain of the allergen species; tests strain-to-strain stability of calls | **not downloaded**; not found in `Fungi_5k/samples.csv` |
| *A. fumigatus* W72310 | NCBI GCA_040167795.1 (`UCR_Afum_W72310_1.0`, chromosome level, UC Riverside, 2024-06-12); `..._protein.faa.gz` exists on the NCBI FTP | owner's UCR strain | **not downloaded**; who made the gene models was not checked |

Rules for the three *A. fumigatus* strains:
1. Gene model sources differ (UniProt, RefSeq, GenBank). Differences in calls between strains can
   come from annotation. The report states the source of each proteome.
2. A strain-to-strain comparison needs a protein key that is not the ID: use the sequence sha256
   and, for non-identical sequences, the best reciprocal hit.
3. Only Af293 has a step 1 truth set. A1163 and W72310 are for stability and allergen homology, not
   for accuracy numbers.
4. Before they enter the panel, check for duplicate proteins and for IDs with `*` or non-standard
   residues (failure-mode table).

Fungi_5k as a batch is not part of version 1.

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
4. Applicability hides most of the genome. Antigen and the Cys-rich finder apply to *Coccidioides* only, so on other proteomes many calls are `not_assessable`. The report must show the share and the reason.
5. Truth is thin. Allergen has no *Coccidioides* truth. Antigen has *Coccidioides* only.
6. Cached output can go stale. The cache key covers version, parameters and artefact hashes.
7. The checkout can be shared. Check `git branch --show-current` before every commit.

## 8. Order of work

1. This spec (revision 5): reviews 1 and 2 done. A third review only if the owner asks. Then a plan.
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
| D1 | Multi-label or exclusive categories | **Answered 2026-10-04: multi-label.** Each category is its own column. |
| D2 | Categories in version 1 | **Answered 2026-10-04:** the five named. Enzyme classes later. |
| D3 | Name and place | **Answered 2026-10-04:** `cellsurface_sorting_hat`, package `src/cellsurface_sorting_hat/`. |
| D4 | Wait for the step 1 decision | **Answered by D10:** no. Variants are carried; the default gate is R0. |
| D5 | Biosynthesis and remodeling families in version 1 | **Answered by D2:** later. Signaling stays out. |
| D6 | Allergen scope in version 1 | **Revised 2026-10-04 after review 2.** Two tiers: `allergen_homolog_hit` (evidence, 35% identity over 80 aa or more, always written) and `allergen_candidate` (>= 70% identity and >= 80% coverage of the allergen length, or an allergen-specific Pfam hit). No surface gate. `unvalidated`. Reason: the 35%/80 aa rule hits 108 of 9,161 Af293 proteins, of which 79 are not the known allergens (mostly housekeeping paralogs) (review 2, one BLAST run, not re-derived). Hits by identity: >= 50%: 74, >= 70%: 42, >= 95%: 29. |
| D7 | Out-of-scope handling | **Revised 2026-10-04 after review 2.** `--taxon` and/or `--taxon-map` is required. A module that is *not applicable* gives `not_assessable`. A module that is applicable but unmeasured gives a value with status `unvalidated`. |
| D8 | Execution engine | **Answered 2026-10-04:** Python driver with SLURM scripts. Nextflow later if the module count grows. |
| D9 | Review model | Different model from the author. Review 1 (Opus) and review 2 (Fable) done. |
| D10 | Default gate for the gated calls (`surface_glycoprotein[v]` used by `cell_wall_adhesion_candidate[v]` and `antigen_candidate_surface[v]`) | **Answered 2026-10-04: `step1_rule@R0`.** R0 has no fitted parameter and runs today (recall 0.603 / 0.727 / 0.938, FPR 0.037 / 0.010 / 0.083). The headline gated columns use R0. R1, R2 and ML variants are extra columns once frozen. Ungated calls are always written. The default is a config item; the owner can change it without a design change. |
| D11 | Add a `cell_wall_protein` call for non-adhesive structural wall proteins | **Answered 2026-10-04: not in version 1.** The report header states the limit. |
| D12 | Antigen call definition | **Answered 2026-10-04:** combined `percentile` (all 9,139 ranked proteins) at most **P = 15** (1,371 proteins); the separate axes are written, not required. The tiers are not used because they exclude SOWgp. Consequence: P = 15 includes PRA2 (percentile 10.79), which the top-decile acceptance test excluded. PRA2 can no longer test the cut, and the 3-of-4 anchor result is reported at the top decile as before. P is a config item. |
| D13 | First proteomes for building and demonstrating version 1 | **Answered 2026-10-04:** *C. immitis* RS, *A. fumigatus* Af293, *S. cerevisiae* S288C, *A. fumigatus* A1163, *A. fumigatus* W72310 (section 4.1). |

## 10. Deliverables

1. This spec (revision 5), `2026-10-04-orchestrator-design-review-1.md` and `-review-2.md`.
2. `docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md` and
   `docs/reports/2026-10-04-fungal-allergen-scoping.md` (written).
3. A plan (not written). Code only after plan review.
