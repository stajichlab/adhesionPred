# cellsurface_sorting_hat module wrappers and calibration: Implementation Plan (Plan 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Tasks 0 to 10 are code with tests that run anywhere. Tasks 11 to 15 run on the HPCC (not in CI) and write measured numbers to a report.

**Goal:** Produce the module tables that the core engine reads (SignalP rule R0, Pfam families, two repeat detectors, allergen homology, antigen ranking lookup, Cys-rich tiers, spherule expression, TM helices), and record for each module where it was calibrated and what sensitivity and specificity were measured there.

**Architecture:** One small wrapper per tool turns the tool's output into a module table plus a run record (the contract of Plan 1). Heavy tools run as SLURM jobs from scripts in `scripts/sorting_hat/`; parsing and decisions run in Python and are tested on small fixtures. A calibration package computes sensitivity and specificity with intervals that respect homology clusters, writes them into the `measure` object of the status sources (one entry per species), and refuses to call a module `estimated` unless both rates were measured on enough independent proteins. A first task patches Plan 1 files: it renames the output calls so that their names say what the tools measure, and adds the report wording.

**Tech Stack:** Python 3.11+ (`/usr/bin/python3.12`), numpy, PyYAML, pytest, ruff. Tools on HPCC: SignalP 6 GPU, hmmer 3.4, ncbi-blast 2.14.0+, TMHMM 2.0c. No new Python dependency.

**Spec:** `docs/superpowers/specs/2026-10-04-orchestrator-design.md` (revision 5; see "Changes to the spec" below). **Plan 1:** `docs/superpowers/plans/2026-10-04-cellsurface-sorting-hat-core.md` (must be implemented first). Reviews of this plan: `docs/superpowers/plans/2026-10-05-cellsurface-sorting-hat-modules-and-calibration-review-1.md`.

## Changes to the spec (owner decisions of 2026-10-05, after the reviews of this plan)

| Spec name (revision 5) | Name in this plan and in the code | Reason |
|---|---|---|
| `surface_glycoprotein[v]` | `signal_peptide_protein[v]` | R0 means "SignalP calls a signal peptide". Such proteins are secreted, wall-bound or GPI-anchored; not shown surface-exposed; glycosylation not assessed; plasma membrane mucins are missed. |
| `adhesion_repeat` | `tandem_repeat_protein` (evidence) | 20 of 27 repeat calls in *C. immitis* RS have no signal peptide (ubiquitin, calmodulin, ankyrin proteins). The adhesion call needs a signal peptide. |
| `adhesion_domain` | `wall_family_domain` (evidence; the `families` column names the family) | CFEM, Bys1, hydrophobin, Als and Pir families are linked to adhesion or wall function in at least one species; function is not shown here. |
| `cell_wall_adhesion_candidate[v]` | unchanged = (`tandem_repeat_protein` OR `wall_family_domain`) AND `signal_peptide_protein[v]`; the ungated form is removed | |
| `antigen_candidate` | `cocci_specificity_rank_top15` | It is a genus-specificity ranking, not epitope prediction; cut set after the anchors were seen. |
| `antigen_candidate_surface[v]` | `serodiagnostic_marker_candidate[v]` (headline) | Owner decision: the list serves serodiagnostic markers; a signal peptide is required. |
| `allergen_homolog_hit` | `iuis_allergen_similarity` (evidence; identity >= 35%, aligned length >= 80 aa, ONE local alignment; not the sliding-window Codex rule) | |
| `allergen_candidate` | `iuis_allergen_homolog` (evidence; identity >= 70% and coverage >= 80%, or an allergen-specific Pfam hit) | IgE cross-reactivity is a hypothesis at most. It is not a mechanism category and is removed from the `other_*` formulas. |
| `other_*` mechanism calls | `tandem_repeat_protein`, `wall_family_domain`, `cocci_specificity_rank_top15` | |
| Status of R0 | one entry **per species** (Phase C species test sets), not pooled | A pooled value describes neither species. |
| Af293 for calibration and the run-level check | UniProt UP000002530 first; the Fungi_5k Af293 file second, as another annotation | The Phase C truth was built on the UniProt gene models; only 4,743 of 9,161 Fungi_5k proteins are identical to them. |

Apply the same renames to the spec file (`sed` with the table above) in the pull request of this plan, and add the sentence "status entries are per species" to spec 3.2.

## Global Constraints

- The module table, run record and status source follow the "Module output contract" of Plan 1. Module tables are `.tsv.gz`; writes are atomic. A table with no usable result is refused (exit 2); some `error` rows make the run state `partial`; a module with no active family writes `unavailable` rows and run state `unavailable`.
- A module's identity (`version`, `params_hash`, `artefact_hash`) belongs to the **tool, database and parameters**, never to one input file. A measurement of a rule must stay valid for another proteome run with the same tool version.
- Status values: `estimated` needs at least 20 positives and 20 negatives, at least 20 independent clusters of each (when recorded), and a 95% interval half-width of at most 0.10 for **both** sensitivity and specificity. A module that was never tested on negatives is at most `smoke`. The status of a Phase C measurement is the weaker of the Phase C label and this rule. A truth set that helped to set the rule or its cutoffs, or that overlaps the module's reference set, caps the status at `smoke` (`--leakage` is required: `none`, `partial`, `tuned_on_truth`, `in_reference`, `unknown`).
- Intervals: the widest of the cluster-bootstrap percentile interval (2,000 resamples, fixed seed) and the Wilson interval on the number of independent clusters of that class. Specificity of Phase C comes with the per-stratum rates (N-int, N-sec, PM-TM); the pooled value depends on the mix of negatives.
- Specificity is never inferred from absence in a database. A negative needs an independent reason to be a negative; if there is none, the module gets no specificity and the report says so.
- A status entry is one species (resolved from `names.dmp`; a name with zero or several IDs stops the run). A protein's taxon may be a strain below the species. A status written for an older module identity is dropped when a new entry is merged.
- The annotation source is part of the claim: R0 was measured on UniProt gene models. A run on another annotation of the same species reports the R0 agreement on shared genes and is not `estimated` for the new gene models.
- Antigen, Cys-rich and expression lookups apply only to the taxon IDs in `--applicable-taxa` (default *C. immitis* RS, 246410, exact match; the tables belong to one annotation).
- SLURM scripts: `$SCRATCH` (`${SCRATCH:?}`), a `trap` that removes the scratch directory, copy-back with a temporary name then `mv`, never `BASH_SOURCE`, never `/tmp`. Run `hmmfetch -f` with model **names**, not accessions.
- A tool result file must belong to the FASTA that is run: the domain table must end with `# [ok]` and show `--cut_ga`; every domain-table target and BLAST query must be a FASTA ID.
- Use `/usr/bin/python3.12` for repository scripts. `pip` and `python` on the login node are Python 3.9.
- Report and docs text: Simplified Technical English. Say when a number was not measured. The report header carries the research-use statement (Task 0).

## What is calibrated, where, and what can be measured

"Sn" is sensitivity, "Sp" is specificity. *none* means that no truth exists today. Numbers marked (measured) were read from `_workdir/step1_compare/phasec/metrics.json` or counted from files on 2026-10-05.

| Module (call) | Calibration set and truth | Sn | Sp | Status now | What limits it |
|---|---|---|---|---|---|
| `step1_rule@R0` (`signal_peptide_protein`) | Phase C GO direct-evidence truth, **per species**, on UniProt gene models. *S. cerevisiae* 79 pos / 3,785 neg; *C. albicans* 153 / 459; *A. fumigatus* 19 / 45; *A. nidulans* 109 / 164; *C. neoformans* 7 / 32; *U. maydis* 9 / 28 | 0.848 [0.738, 0.943]; 0.477 [0.356, 0.591]; 0.947 [0.833, 1.000]; 0.688 [0.595, 0.779]; 0.857 [0.500, 1.000]; 1.000 [1.000, 1.000] (measured; the last two are widened with Wilson) | 0.966 [0.958, 0.972]; 0.943 [0.913, 0.968]; 1.000 (Wilson-widened); 0.988 [0.968, 1.000]; 0.938 [0.839, 1.000]; 0.893 [0.758, 1.000] (measured) | `smoke` for five species; `estimated` only for *A. nidulans* (the one species Phase C labels `estimate`) | Pooled sets hide a large difference (*S. cerevisiae* 0.85 against *C. albicans* 0.48). The pooled specificity is dominated by intracellular negatives; N-sec specificity is 0.91 (*S. cerevisiae*) and 0.93 (*C. albicans*). SignalP 6 was trained on proteins with experimental evidence; its overlap with the Phase C positives is not measured (`leakage = unknown` for R0). *Coccidioides* (Onygenales) is not covered: `unvalidated` there. |
| `tandem_repeat_protein` (`repeat02`, `repeat14`) | None with clade truth (synthetic series plus SOWgp only). Owner decision C1. | none | none | `unvalidated` | Calling known repeat adhesins is circular for settings tuned on SOWgp. 27 calls in RS, 20 without a signal peptide. |
| `wall_family_domain` (`pfam_adhesion`) | Per family: hits on the proteomes of 4 clades; non-member hits read by a person (Task 13). | per family | counts of reviewed non-member hits (no number called specificity) | `unvalidated`; all families start inactive | "Members known by function" are 0 to 3 per proteome; most non-member hits are uncharacterised family members. Domain presence is not function. Pth11-like receptors carry CFEM (`no_tm` removes them; helices that start inside the first 35 residues are not counted). |
| `allergen_homology` (`iuis_allergen_similarity`, `iuis_allergen_homolog`) | WHO/IUIS fungal sequences (111 usable of 116 with a sequence), leave-species-out against themselves (the 35%/80 aa rule and the 70%/80% rule, exactly as the engine applies them) | recovered share per rule (Task 13) | **none** | `unvalidated` | The reference set holds the 30 *A. fumigatus* allergens, so Af293 counts show self-recognition. Absence from IUIS means "never tested". WHO/IUIS lists no *Coccidioides* allergen; the only Onygenales entries are four *Trichophyton* proteases (contact route). |
| `antigen_lookup` (`cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate`) | Four *Coccidioides* anchors (PRA3, Ag2/PRA, SOWgp, PRA2; Ag2/PRA, PRA2, PRA3 are one protein family), CF antigen as cross-reactive control | 3 of 4 anchors in the top decile (pre-set test); the top-15% cut was set after the anchors were seen and is not supported by them | none independent | `smoke`, `leakage = tuned_on_truth` | 15% of the proteome is called; 143 of 1,371 (10.4%) of the calls have a signal peptide, hence the gated headline column. The ranking is similarity to IEDB antigens plus prevalence plus absence of orthologs; it is not epitope prediction and no antibody or T-cell measurement supports it. |
| `cys_rich` | No accuracy claim (finder README). | n/a | n/a | evidence only | |
| `expression` | RNA-seq, spherule against mycelium, 2 replicates. | n/a | n/a | evidence only | Host-phase evidence; says nothing about allergen exposure (conidia). |
| `tm` (TMHMM) | n/a here. | n/a | n/a | evidence only | TMHMM can read a signal peptide as a helix; `n_tm_mature` ignores helices that start at or before residue 35. |

Owner decisions C1 to C4 (answered 2026-10-05). They set what the calibration tasks use; the tasks that need new data are not code tasks and are not done yet.
- **C1** Truth for `tandem_repeat_protein`: curated 2a repeat adhesins (`data/curated/adhesins`) as positives against curated hard negatives and secreted non-repeat proteins, clustered by homology; SOWgp is marked `tuned` and excluded from the metrics. Family members for `wall_family_domain`: full-length domain plus the conserved Cys pattern (8 Cys for CFEM and hydrophobins) or a curated function; hits without the pattern are errors; uncharacterised pattern-positive hits are their own class.
- **C2** Allergen negatives: IEDB IgE assay records with a negative outcome for fungi, clustered with the positives. First step: count them. If there are too few, the module stays without a specificity.
- **C3** Independent truth for `serodiagnostic_marker_candidate`: a leave-species-out test (rebuild the ranking with one species' IEDB antigens held out: *A. fumigatus*, *H. capsulatum* or *C. albicans*). Feasibility is not checked (antigen counts per species; a comparable confounder set). Serology is future work.
- **C4** Sign-off of a Pfam family: the owner, per family, after reading the review table; the change is committed alone with the non-member hit table.

## Review Focus

1. A result file whose IDs do not match the FASTA (SignalP keeps the whole header; other tools cut at a space), an unfinished or empty domain table, or a BLAST table from another proteome must stop the module with exit code 2 and write nothing. (Task 9: `test_a_result_file_whose_ids_do_not_match_the_fasta_is_refused`, `test_pfam_command_refuses_a_domain_table_from_another_proteome`, `test_pfam_command_refuses_an_unfinished_domain_table`, `test_allergen_command_refuses_a_blast_table_from_another_proteome`, `test_a_lookup_that_matches_no_applicable_protein_is_refused`)
2. A measurement of rule R0 must still apply to a different proteome run with the same SignalP version. (Task 9: `test_the_signalp_module_identity_does_not_depend_on_the_proteome`)
3. A rule never tested on negatives, a truth set with too few independent clusters, or a truth set that overlaps the tuning data must not be called `estimated`. (Task 6: `test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate`, `test_one_cluster_per_class_cannot_give_an_estimate`, `test_a_tight_interval_with_few_recorded_clusters_is_not_an_estimate`, `test_all_positives_called_gives_a_width_not_a_zero_width_interval`; Task 9: `test_leakage_other_than_none_caps_an_estimate_at_smoke`)
4. A CFEM domain in a transmembrane receptor must not count as a family hit, and a signal peptide read as a helix must not remove a real wall protein. (Task 5: `test_tm_rows_count_helices_after_the_signal_peptide_window`; Task 2: `test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices`; Task 9: `test_pfam_command_applies_the_no_tm_condition_from_the_tm_table`)
5. While no Pfam family is active, the domain modules must be `unavailable`, not a table of zeros. (Task 9: `test_pfam_command_writes_unavailable_modules_when_no_family_is_active`)
6. The allergen recall must be computed against allergens of **other species**, with the rules the engine applies and the denominator taken from the FASTA; the allergen FASTA must not contain free text. (Task 4: `test_recall_by_rule_uses_the_rules_of_the_engine_and_the_fasta_denominator`, `test_lso_report_refuses_a_sequence_with_no_blast_line`, `test_allergen_fasta_is_clean_and_has_a_meta_table`)
7. A status measured on one module version must not be attached to the next. (Task 9: `test_a_changed_module_identity_drops_the_old_status_entries`)
8. A species estimate must not be replaced by a pooled value. (Task 7: `test_each_species_gets_its_own_entry_with_its_own_numbers`)
9. Proteins shorter than the detector minimum must not become `error`. (Task 3: `test_a_protein_shorter_than_the_detector_minimum_is_not_called_not_an_error`)
10. Two strains of one species: calls for identical sequences are compared by sequence hash, not ID. (Task 7: `test_only_identical_sequences_are_compared_by_sha256_not_by_id`)

---

### Task 0: Branch, and the rename patch for Plan 1 files

**Files:**
- Modify (by patch): `src/cellsurface_sorting_hat/categories.yaml`, `src/cellsurface_sorting_hat/outputs.py`, `src/cellsurface_sorting_hat/status.py`, `tests/cellsurface_sorting_hat/test_engine.py`, `test_cli.py`, `test_outputs.py`, `test_taxonomy_status.py`

**Interfaces:**
- Consumes: the Plan 1 implementation (merged).
- Produces: the call names of the table "Changes to the spec"; the `evidence:` list with the allergen, antigen and TM fields; the report header limits; `strata`, `n_clusters_pos` and `n_clusters_neg` as allowed keys of a `measure` (checked); per-stratum rates validated like other rates.

- [ ] **Step 1: Create the working branch**

```bash
git fetch origin
git checkout -b sorting-hat-modules origin/main      # Plan 1 (sorting-hat-core) must be merged first
git branch --show-current                            # expected: sorting-hat-modules
ls src/cellsurface_sorting_hat/engine.py             # expected: the file exists (Plan 1)
mkdir -p src/cellsurface_sorting_hat/modules src/cellsurface_sorting_hat/calibration
mkdir -p tests/cellsurface_sorting_hat/modules tests/cellsurface_sorting_hat/calibration
mkdir -p data/sorting_hat scripts/sorting_hat
```

If Plan 1 is not merged, branch from the Plan 1 branch tip instead and say so in the pull request.

- [ ] **Step 2: Save the patch and check that it applies**

Save the text below as `sorting_hat_plan2_rename.patch` in the repository root (do not commit it). It was made with `diff -u` from the Plan 1 files as they are written in Plan 1.

```diff
--- a/src/cellsurface_sorting_hat/categories.yaml
+++ b/src/cellsurface_sorting_hat/categories.yaml
@@ -29,43 +29,45 @@
   - antigen_lookup.prevalence
   - antigen_lookup.max_crossreact
   - allergen_homology.allergen_name
+  - allergen_homology.allergen_species
+  - allergen_homology.exposure
+  - allergen_homology.evidence
   - allergen_homology.identity
   - allergen_homology.coverage
   - allergen_homology.aligned_length
   - cys_rich.tier
   - expression.log2fc
+  - tm.n_tm
+  - tm.n_tm_mature
 
 calls:
-  - name: surface_glycoprotein
+  - name: signal_peptide_protein
     per_variant: true
     expr: {call: "{step1}"}
-  - name: adhesion_repeat
+  - name: tandem_repeat_protein
     expr:
       or: [{call: repeat02}, {call: repeat14}]
-  - name: adhesion_domain
+  - name: wall_family_domain
     expr: {flag: pfam_adhesion.hit}
-  - name: cell_wall_adhesion_ungated
-    expr:
-      or: [{ref: adhesion_repeat}, {ref: adhesion_domain}]
   - name: cell_wall_adhesion_candidate
     per_variant: true
     expr:
       and:
-        - or: [{ref: adhesion_repeat}, {ref: adhesion_domain}]
-        - {ref: surface_glycoprotein}
-  - name: antigen_candidate
+        - or: [{ref: tandem_repeat_protein}, {ref: wall_family_domain}]
+        - {ref: signal_peptide_protein}
+  - name: cocci_specificity_rank_top15
     expr:
       test: {module: antigen_lookup, field: percentile, op: "<=", value: "$antigen_percentile_max"}
-  - name: antigen_candidate_surface
+  - name: serodiagnostic_marker_candidate
     per_variant: true
     expr:
-      and: [{ref: antigen_candidate}, {ref: surface_glycoprotein}]
-  - name: allergen_homolog_hit
+      and: [{ref: cocci_specificity_rank_top15}, {ref: signal_peptide_protein}]
+  - name: iuis_allergen_similarity
     expr:
       and:
         - test: {module: allergen_homology, field: identity, op: ">=", value: "$allergen_hit_identity_min"}
         - test: {module: allergen_homology, field: aligned_length, op: ">=", value: "$allergen_hit_length_min"}
-  - name: allergen_candidate
+  - name: iuis_allergen_homolog
     expr:
       or:
         - and:
@@ -75,12 +77,12 @@
   - name: other_not_surface
     per_variant: true
     kind: other
-    surface: surface_glycoprotein
+    surface: signal_peptide_protein
     surface_is: not_called
-    mechanism: [adhesion_repeat, adhesion_domain, antigen_candidate, allergen_candidate]
+    mechanism: [tandem_repeat_protein, wall_family_domain, cocci_specificity_rank_top15]
   - name: other_surface_no_mechanism
     per_variant: true
     kind: other
-    surface: surface_glycoprotein
+    surface: signal_peptide_protein
     surface_is: called
-    mechanism: [adhesion_repeat, adhesion_domain, antigen_candidate, allergen_candidate]
+    mechanism: [tandem_repeat_protein, wall_family_domain, cocci_specificity_rank_top15]
--- a/src/cellsurface_sorting_hat/outputs.py
+++ b/src/cellsurface_sorting_hat/outputs.py
@@ -14,18 +14,28 @@
 MAX_LISTED_INVALID = 20
 
 KNOWN_LIMITS = [
-    "GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get only "
-    "`surface_glycoprotein`. There is no `cell_wall_protein` call in version 1.",
-    "`surface_glycoprotein` is defined by GO cell wall and extracellular region evidence. It is not "
-    'evidence of glycosylation. With the step 1 rule R0 the call means "SignalP calls a signal '
-    'peptide" and nothing more.',
-    "CFEM is filed under adhesion because class 2b-i is. Its confirmed fold is a hemophore. Binding "
-    "to a host receptor is not shown.",
-    "The repeat detectors have no clade truth set. Their calls are hypotheses.",
-    "The antigen call is the top 15% of a fixed Coccidioides ranking. It is a weak label, and the "
-    "ranking prints NOT CALIBRATED (3 of 4 anchors pass the top-decile test).",
+    "Research use only. This is not a regulatory allergenicity assessment and not a diagnostic "
+    "result. No row is supported by an IgE, antibody or T-cell measurement.",
+    "`signal_peptide_protein` means that SignalP calls a signal peptide (rule R0) and nothing more. "
+    "Such proteins are secreted, wall-bound or GPI-anchored. They are not shown to be exposed at the "
+    "cell surface, and glycosylation is not assessed. Plasma membrane mucins can be missed.",
+    "GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get no finer "
+    "label in version 1. There is no `cell_wall_protein` call.",
+    "`tandem_repeat_protein` and `wall_family_domain` are evidence. Repeat proteins include "
+    "intracellular ones (ubiquitin, calmodulin, ankyrin proteins). A domain of a family that is "
+    "linked to adhesion or wall function in at least one species (CFEM, Bys1, hydrophobin, Als) is "
+    "not shown to mediate adhesion here. The adhesion call needs a signal peptide.",
+    "`cocci_specificity_rank_top15` is the top 15% of a fixed Coccidioides immitis ranking "
+    "(similarity to IEDB antigens, prevalence, absence of orthologs in confounder fungi). It is not "
+    "epitope prediction. The cut was set after the four anchors were seen. The ranking prints NOT "
+    "CALIBRATED. `serodiagnostic_marker_candidate` adds a signal peptide. Peptide level only; glycan "
+    "epitopes are not assessed.",
+    "`iuis_allergen_similarity` and `iuis_allergen_homolog` are sequence similarity to allergens in "
+    "the WHO/IUIS fungal set (IgE binding in patients). They suggest possible IgE cross-reactivity "
+    "at most. A protein with no hit is not thereby a non-allergen. WHO/IUIS lists no Coccidioides "
+    "allergen.",
     "Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, "
-    "moonlighting proteins and biofilm are not categories.",
+    "non-protein adhesins, moonlighting proteins and biofilm are not categories.",
     "The taxon you give is recorded as given. It is not checked against the sequences.",
 ]
 
--- a/src/cellsurface_sorting_hat/status.py
+++ b/src/cellsurface_sorting_hat/status.py
@@ -19,6 +19,9 @@
     "truth_source",
     "n_pos",
     "n_neg",
+    "n_clusters_pos",
+    "n_clusters_neg",
+    "strata",
     "sensitivity",
     "specificity",
     "notes",
@@ -77,7 +80,12 @@
     for key in ("sensitivity", "specificity"):
         if key in measure:
             _check_rate(measure[key], f"{where}.{key}")
-    for key in ("n_pos", "n_neg"):
+    strata = measure.get("strata", {})
+    if not isinstance(strata, dict):
+        raise ValueError(f"{where}.strata: must be an object of rates")
+    for name, rate in strata.items():
+        _check_rate(rate, f"{where}.strata.{name}")
+    for key in ("n_pos", "n_neg", "n_clusters_pos", "n_clusters_neg"):
         if key in measure and not (isinstance(measure[key], int) and measure[key] >= 0):
             raise ValueError(f"{where}.{key}: must be a non-negative integer")
 
--- a/tests/cellsurface_sorting_hat/test_engine.py
+++ b/tests/cellsurface_sorting_hat/test_engine.py
@@ -51,15 +51,14 @@
 def test_packaged_config_loads_and_has_the_spec_calls():
     names = [c["name"] for c in load_config().calls]
     assert names == [
-        "surface_glycoprotein",
-        "adhesion_repeat",
-        "adhesion_domain",
-        "cell_wall_adhesion_ungated",
+        "signal_peptide_protein",
+        "tandem_repeat_protein",
+        "wall_family_domain",
         "cell_wall_adhesion_candidate",
-        "antigen_candidate",
-        "antigen_candidate_surface",
-        "allergen_homolog_hit",
-        "allergen_candidate",
+        "cocci_specificity_rank_top15",
+        "serodiagnostic_marker_candidate",
+        "iuis_allergen_similarity",
+        "iuis_allergen_homolog",
         "other_not_surface",
         "other_surface_no_mechanism",
     ]
@@ -67,40 +66,40 @@
 
 def test_only_available_step1_variants_get_per_variant_records():
     res = run(base_modules())
-    variants = {v for (_, call, v) in res if call == "surface_glycoprotein"}
+    variants = {v for (_, call, v) in res if call == "signal_peptide_protein"}
     assert variants == {"R0"}
 
 
 def test_a_true_or_input_wins_over_an_unknown_input():
     mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
     mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
-    r = run(mods)[("P", "adhesion_repeat", "")]
+    r = run(mods)[("P", "tandem_repeat_protein", "")]
     assert r.value == "called"
 
 
 def test_unknown_or_input_with_a_false_other_is_unknown():
     mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
-    assert run(mods)[("P", "adhesion_repeat", "")].value == "not_assessable"
+    assert run(mods)[("P", "tandem_repeat_protein", "")].value == "not_assessable"
 
 
 def test_false_surface_makes_gated_antigen_false_even_when_antigen_is_unknown():
     mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
     res = run(mods)
-    assert res[("P", "antigen_candidate", "")].value == "not_assessable"
-    assert res[("P", "antigen_candidate_surface", "R0")].value == "not_called"
+    assert res[("P", "cocci_specificity_rank_top15", "")].value == "not_assessable"
+    assert res[("P", "serodiagnostic_marker_candidate", "R0")].value == "not_called"
 
 
 def test_antigen_threshold_is_inclusive():
     mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15")}))
-    assert run(mods)[("P", "antigen_candidate", "")].value == "called"
+    assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "called"
     mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15.01")}))
-    assert run(mods)[("P", "antigen_candidate", "")].value == "not_called"
+    assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "not_called"
 
 
 def test_empty_or_nan_number_is_unknown():
     for bad in ("", "nan", "abc"):
         mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
-        assert run(mods)[("P", "antigen_candidate", "")].value == "not_assessable", bad
+        assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "not_assessable", bad
 
 
 @pytest.mark.parametrize(
@@ -120,28 +119,35 @@
         ),
         pfam_allergen=table("pfam_allergen", {"P": ok(hit=pfam)}),
     )
-    assert run(mods)[("P", "allergen_candidate", "")].value == expected
+    assert run(mods)[("P", "iuis_allergen_homolog", "")].value == expected
 
 
 def test_allergen_needs_no_surface_call():
     mods = base_modules(step1="not_called")
     mods["allergen_homology"] = table("allergen_homology", {"P": ok(identity="90", coverage="95")})
     res = run(mods)
-    assert res[("P", "surface_glycoprotein", "R0")].value == "not_called"
-    assert res[("P", "allergen_candidate", "")].value == "called"
+    assert res[("P", "signal_peptide_protein", "R0")].value == "not_called"
+    assert res[("P", "iuis_allergen_homolog", "")].value == "called"
 
 
 def test_other_not_surface_leaves_out_unknown_calls_and_names_them():
     mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
     r = run(mods)[("P", "other_not_surface", "R0")]
-    assert (r.value, r.other_basis) == ("called", "antigen_candidate")
+    assert (r.value, r.other_basis) == ("called", "cocci_specificity_rank_top15")
 
 
 def test_other_not_surface_is_false_when_a_mechanism_call_is_true():
+    mods = base_modules(repeat02=table("repeat02", {"P": ok(call="called")}))
+    assert run(mods)[("P", "other_not_surface", "R0")].value == "not_called"
+
+
+def test_an_allergen_flag_does_not_change_the_other_calls():
     mods = base_modules(
         allergen_homology=table("allergen_homology", {"P": ok(identity="90", coverage="95")})
     )
-    assert run(mods)[("P", "other_not_surface", "R0")].value == "not_called"
+    res = run(mods)
+    assert res[("P", "iuis_allergen_homolog", "")].value == "called"
+    assert res[("P", "other_not_surface", "R0")].value == "called"
 
 
 def test_other_surface_no_mechanism():
@@ -175,8 +181,8 @@
     mods = base_modules(step1="called")
     mods["repeat02"] = table("repeat02", {"P": ok(call="called")})  # repeat14 stays not_called
     res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
-    # adhesion_repeat is true because of repeat02 only: status smoke, not unvalidated
-    assert res[("P", "adhesion_repeat", "")].status == "smoke"
+    # tandem_repeat_protein is true because of repeat02 only: status smoke, not unvalidated
+    assert res[("P", "tandem_repeat_protein", "")].status == "smoke"
     # the gated call is true because of repeat02 and step 1: weakest is smoke
     r = res[("P", "cell_wall_adhesion_candidate", "R0")]
     assert (r.value, r.status) == ("called", "smoke")
@@ -207,10 +213,10 @@
         (lambda d: d["calls"].append(dict(d["calls"][0])), "duplicate call"),
         (lambda d: d["calls"][1].update(expr={"ref": "later_call"}), "unknown or later"),
         (lambda d: d["calls"][1].update(expr={"bogus": 1}), "bad node"),
-        (lambda d: d["calls"][5]["expr"]["test"].update(op="~"), "bad operator"),
-        (lambda d: d["calls"][5]["expr"]["test"].update(value="$nope"), "unknown threshold"),
+        (lambda d: d["calls"][4]["expr"]["test"].update(op="~"), "bad operator"),
+        (lambda d: d["calls"][4]["expr"]["test"].update(value="$nope"), "unknown threshold"),
         (
-            lambda d: d["calls"][1].update(expr={"ref": "surface_glycoprotein"}),
+            lambda d: d["calls"][1].update(expr={"ref": "signal_peptide_protein"}),
             "ungated call cannot use",
         ),
     ],
@@ -230,7 +236,7 @@
     assert (r.value, r.status, r.status_basis) == ("not_assessable", "unvalidated", "")
     # an OR with a false input and an unknown input is unknown and also has no deciding module
     mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
-    r = run(mods, status_of=lambda m, t: ("estimated", "t"))[("P", "adhesion_repeat", "")]
+    r = run(mods, status_of=lambda m, t: ("estimated", "t"))[("P", "tandem_repeat_protein", "")]
     assert (r.value, r.status) == ("not_assessable", "unvalidated")
 
 
@@ -243,13 +249,13 @@
         ("", "90", "not_assessable"),
     ],
 )
-def test_allergen_homolog_hit_is_the_35_percent_80_aa_evidence_tier(identity, length, expected):
+def test_iuis_allergen_similarity_is_the_35_percent_80_aa_evidence_tier(identity, length, expected):
     mods = base_modules(
         allergen_homology=table(
             "allergen_homology", {"P": ok(identity=identity, aligned_length=length, coverage="0")}
         )
     )
-    assert run(mods)[("P", "allergen_homolog_hit", "")].value == expected
+    assert run(mods)[("P", "iuis_allergen_similarity", "")].value == expected
 
 
 def test_evidence_rows_are_collected_only_for_ok_rows_with_a_value():
@@ -285,7 +291,7 @@
     "mutate,message",
     [
         (
-            lambda d: d["calls"][-2].update(surface="adhesion_repeat"),
+            lambda d: d["calls"][-2].update(surface="tandem_repeat_protein"),
             "needs an earlier per_variant surface",
         ),
         (
@@ -303,7 +309,7 @@
 @pytest.mark.parametrize("bad", ["inf", "-inf", "Infinity"])
 def test_a_non_finite_number_is_unknown(bad):
     mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
-    assert run(mods)[("P", "antigen_candidate", "")].value == "not_assessable"
+    assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "not_assessable"
 
 
 def test_a_false_and_takes_its_status_from_the_false_inputs_only():
@@ -311,6 +317,6 @@
     statuses = {R0: ("smoke", "t"), "antigen_lookup": ("estimated", "t")}
     mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")}))
     res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
-    r = res[("P", "antigen_candidate_surface", "R0")]
+    r = res[("P", "serodiagnostic_marker_candidate", "R0")]
     assert (r.value, r.status) == ("not_called", "smoke")
     assert r.status_basis == "step1_rule@R0:t"
--- a/tests/cellsurface_sorting_hat/test_cli.py
+++ b/tests/cellsurface_sorting_hat/test_cli.py
@@ -97,26 +97,27 @@
         return r["value"], r["status"]
 
     # SOW1: Coccidioides (taxon 41), not covered by the step 1 estimate (tested on taxon 40)
-    assert v("SOW1", "surface_glycoprotein", "R0") == ("called", "unvalidated")
-    assert v("SOW1", "adhesion_repeat") == ("called", "unvalidated")
+    assert v("SOW1", "signal_peptide_protein", "R0") == ("called", "unvalidated")
+    assert v("SOW1", "tandem_repeat_protein") == ("called", "unvalidated")
     assert v("SOW1", "cell_wall_adhesion_candidate", "R0") == ("called", "unvalidated")
-    assert v("SOW1", "antigen_candidate") == ("called", "unvalidated")
-    assert v("SOW1", "antigen_candidate_surface", "R0") == ("called", "unvalidated")
-    assert v("SOW1", "allergen_candidate")[0] == "not_called"
+    assert v("SOW1", "cocci_specificity_rank_top15") == ("called", "unvalidated")
+    assert v("SOW1", "serodiagnostic_marker_candidate", "R0") == ("called", "unvalidated")
+    assert v("SOW1", "iuis_allergen_homolog")[0] == "not_called"
     assert v("SOW1", "other_not_surface", "R0")[0] == "not_called"
     assert v("SOW1", "other_surface_no_mechanism", "R0")[0] == "not_called"
-    assert got[("SOW1", "surface_glycoprotein", "R0")]["status_basis"] == (
+    assert got[("SOW1", "signal_peptide_protein", "R0")]["status_basis"] == (
         "step1_rule@R0:taxon not tested"
     )
     # ENZ1: A. fumigatus (40), step 1 estimated; antigen not applicable; allergen homology called
-    assert v("ENZ1", "surface_glycoprotein", "R0") == ("not_called", "estimated")
-    assert v("ENZ1", "antigen_candidate")[0] == "not_assessable"
-    assert v("ENZ1", "antigen_candidate_surface", "R0") == ("not_called", "estimated")
-    assert v("ENZ1", "allergen_candidate")[0] == "called"
-    assert v("ENZ1", "other_not_surface", "R0")[0] == "not_called"
+    assert v("ENZ1", "signal_peptide_protein", "R0") == ("not_called", "estimated")
+    assert v("ENZ1", "cocci_specificity_rank_top15")[0] == "not_assessable"
+    assert v("ENZ1", "serodiagnostic_marker_candidate", "R0") == ("not_called", "estimated")
+    assert v("ENZ1", "iuis_allergen_homolog")[0] == "called"
+    r = got[("ENZ1", "other_not_surface", "R0")]  # allergen flags are not mechanism calls
+    assert (r["value"], r["other_basis"]) == ("called", "cocci_specificity_rank_top15")
     # STAR1: surface, no mechanism called, antigen left out of the basis
     r = got[("STAR1", "other_surface_no_mechanism", "R0")]
-    assert (r["value"], r["other_basis"]) == ("called", "antigen_candidate")
+    assert (r["value"], r["other_basis"]) == ("called", "cocci_specificity_rank_top15")
     # BAD1: invalid protein, every call unknown
     assert {r["value"] for k, r in got.items() if k[0] == "BAD1"} == {"not_assessable"}
 
@@ -156,7 +157,7 @@
         ]
     )
     assert code == 0
-    assert read_long(out)[("STAR1", "surface_glycoprotein", "R0")]["value"] == "called"
+    assert read_long(out)[("STAR1", "signal_peptide_protein", "R0")]["value"] == "called"
 
 
 def test_unknown_taxon_stops_the_run(tmp_path, write_module, nodes_dmp, capsys):
@@ -222,8 +223,8 @@
     )
     assert code == 0
     got = read_long(out)
-    assert got[("SOW1", "surface_glycoprotein", "R0")]["status"] == "unvalidated"  # taxon 41
-    assert got[("STAR1", "surface_glycoprotein", "R0")]["status"] == "estimated"  # --taxon 40
+    assert got[("SOW1", "signal_peptide_protein", "R0")]["status"] == "unvalidated"  # taxon 41
+    assert got[("STAR1", "signal_peptide_protein", "R0")]["status"] == "estimated"  # --taxon 40
 
 
 def test_partial_module_output_is_reported_and_missing_proteins_are_unknown(
@@ -250,9 +251,9 @@
     assert "| repeat14 | partial |" in (out / "report.md").read_text()
     got = read_long(out)
     # ENZ1: repeat02 false and repeat14 missing -> unknown, not false
-    assert got[("ENZ1", "adhesion_repeat", "")]["value"] == "not_assessable"
+    assert got[("ENZ1", "tandem_repeat_protein", "")]["value"] == "not_assessable"
     # SOW1: repeat02 true wins over anything
-    assert got[("SOW1", "adhesion_repeat", "")]["value"] == "called"
+    assert got[("SOW1", "tandem_repeat_protein", "")]["value"] == "called"
 
 
 def test_module_with_run_state_unavailable_is_treated_as_absent(tmp_path, write_module, nodes_dmp):
@@ -316,8 +317,8 @@
     assert [r["protein"] for r in rows] == ["SOW1", "ENZ1", "DUP1", "BAD1", "STAR1"]
     star = rows[-1]
     assert star["other_surface_no_mechanism[R0]"] == "called"
-    assert star["other_surface_no_mechanism[R0]_basis"] == "antigen_candidate"
-    assert star["surface_glycoprotein[R0]_status"] == "estimated"
+    assert star["other_surface_no_mechanism[R0]_basis"] == "cocci_specificity_rank_top15"
+    assert star["signal_peptide_protein[R0]_status"] == "estimated"
     run = json.loads((out / "run.json").read_text())
     assert run["n_proteins"] == 5 and run["n_invalid"] == 1
     assert run["taxa"] == {"40": 3, "41": 2}
@@ -327,9 +328,9 @@
     _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
     report = (out / "report.md").read_text()
     assert "## Known limits" in report
-    assert 'means "SignalP calls a signal peptide"' in report
-    assert "| surface_glycoprotein | R0 | 3 | 1 | 1 |" in report
-    assert "- antigen_candidate: " in report
+    assert "means that SignalP calls a signal peptide" in report
+    assert "| signal_peptide_protein | R0 | 3 | 1 | 1 |" in report
+    assert "- cocci_specificity_rank_top15: " in report
 
 
 def _run(tmp_path, nodes_dmp, fasta, taxon_map, wd, out="out"):
@@ -410,8 +411,8 @@
     assert code == 0
     assert "| repeat14 | bad_value | 1 |" in (out / "report.md").read_text()
     got = read_long(out)
-    assert got[("SOW1", "adhesion_repeat", "")]["value"] == "called"  # repeat02 is true
-    assert got[("STAR1", "adhesion_repeat", "")]["value"] == "not_called"
+    assert got[("SOW1", "tandem_repeat_protein", "")]["value"] == "called"  # repeat02 is true
+    assert got[("STAR1", "tandem_repeat_protein", "")]["value"] == "not_called"
 
 
 def test_a_module_whose_ids_match_no_protein_is_an_error(tmp_path, write_module, nodes_dmp):
@@ -465,7 +466,7 @@
         "Modules the rules need and that were not found: pfam_adhesion"
         in (out / "report.md").read_text()
     )
-    assert read_long(out)[("SOW1", "adhesion_domain", "")]["value"] == "not_assessable"
+    assert read_long(out)[("SOW1", "wall_family_domain", "")]["value"] == "not_assessable"
 
 
 @pytest.mark.parametrize(
@@ -576,7 +577,7 @@
     _module_bytes(wd, "repeat14", b"\xef\xbb\xbf" + raw)
     code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
     assert code == 0
-    assert read_long(out)[("SOW1", "surface_glycoprotein", "R0")]["value"] == "called"
+    assert read_long(out)[("SOW1", "signal_peptide_protein", "R0")]["value"] == "called"
 
 
 def test_a_repeated_id_in_the_taxon_map_names_both_lines(tmp_path, write_module, nodes_dmp, capsys):
--- a/tests/cellsurface_sorting_hat/test_outputs.py
+++ b/tests/cellsurface_sorting_hat/test_outputs.py
@@ -41,30 +41,39 @@
 
 
 def test_long_table_has_the_documented_columns(tmp_path):
-    write_long(tmp_path / "l.tsv.gz", [rec("A", "adhesion_repeat", "", "called", "smoke", "m:t")])
+    write_long(
+        tmp_path / "l.tsv.gz", [rec("A", "tandem_repeat_protein", "", "called", "smoke", "m:t")]
+    )
     rows = read(tmp_path / "l.tsv.gz")
     assert rows[0] == LONG_COLUMNS
-    assert rows[1] == ["A", "adhesion_repeat", "", "called", "smoke", "m:t", ""]
+    assert rows[1] == ["A", "tandem_repeat_protein", "", "called", "smoke", "m:t", ""]
 
 
 def test_wide_table_names_gated_columns_with_the_variant(tmp_path):
     records = [
-        rec("A", "surface_glycoprotein", "R0", "called", "estimated"),
-        rec("A", "other_not_surface", "R0", "called", "estimated", other="antigen_candidate"),
-        rec("B", "surface_glycoprotein", "R0", "not_called"),
+        rec("A", "signal_peptide_protein", "R0", "called", "estimated"),
+        rec(
+            "A",
+            "other_not_surface",
+            "R0",
+            "called",
+            "estimated",
+            other="cocci_specificity_rank_top15",
+        ),
+        rec("B", "signal_peptide_protein", "R0", "not_called"),
         rec("B", "other_not_surface", "R0", "not_assessable"),
     ]
     write_wide(tmp_path / "w.tsv.gz", records, ["A", "B"])
     rows = read(tmp_path / "w.tsv.gz")
     assert rows[0] == [
         "protein",
-        "surface_glycoprotein[R0]",
-        "surface_glycoprotein[R0]_status",
+        "signal_peptide_protein[R0]",
+        "signal_peptide_protein[R0]_status",
         "other_not_surface[R0]",
         "other_not_surface[R0]_status",
         "other_not_surface[R0]_basis",
     ]
-    assert rows[1][1:3] == ["called", "estimated"] and rows[1][-1] == "antigen_candidate"
+    assert rows[1][1:3] == ["called", "estimated"] and rows[1][-1] == "cocci_specificity_rank_top15"
 
 
 def test_report_warns_about_error_partial_and_inconsistent_modules():
@@ -74,7 +83,7 @@
             module_notes={"a": "no result table"},
             inconsistent={"c": 2},
         ),
-        [rec("A", "adhesion_repeat", "", "called")],
+        [rec("A", "tandem_repeat_protein", "", "called")],
     )
     assert "WARNING: Module a is in state error. no result table" in text
     assert "WARNING: Module b is in state partial." in text
@@ -90,7 +99,7 @@
 
 def test_report_counts_values_per_call_and_variant():
     records = [
-        rec("A", "surface_glycoprotein", "R0", v) for v in ("called", "called", "not_called")
+        rec("A", "signal_peptide_protein", "R0", v) for v in ("called", "called", "not_called")
     ]
     text = render_report(info(), records)
-    assert "| surface_glycoprotein | R0 | 2 | 1 | 0 |" in text
+    assert "| signal_peptide_protein | R0 | 2 | 1 | 0 |" in text
--- a/tests/cellsurface_sorting_hat/test_taxonomy_status.py
+++ b/tests/cellsurface_sorting_hat/test_taxonomy_status.py
@@ -189,3 +189,20 @@
     path.write_text("1\t|\t1\t|\tno rank\t|\n2\t|\tx\t|\tno rank\t|\n")
     with pytest.raises(TaxonError, match="nodes.dmp:2"):
         Lineage.from_nodes_dmp(path)
+
+
+def test_per_stratum_rates_and_cluster_counts_are_accepted_and_checked(tmp_path):
+    good = {
+        **MEASURE,
+        "n_clusters_pos": 40,
+        "n_clusters_neg": 300,
+        "strata": {"N-sec": {"value": 0.9, "lo": 0.8, "hi": 0.95}},
+    }
+    assert (
+        load_status_source(_status_file(tmp_path, good)).entries[0].measure["n_clusters_pos"] == 40
+    )
+    bad = {**MEASURE, "strata": {"N-sec": {"value": 2, "lo": 0, "hi": 3}}}
+    with pytest.raises(ValueError):
+        load_status_source(_status_file(tmp_path, bad))
+    with pytest.raises(ValueError):
+        load_status_source(_status_file(tmp_path, {**MEASURE, "n_clusters_pos": -1}))
```

```bash
patch -p1 --dry-run < sorting_hat_plan2_rename.patch    # expected: 7 files "checking file ...", no FAILED hunk
patch -p1 < sorting_hat_plan2_rename.patch
rm sorting_hat_plan2_rename.patch
```

- [ ] **Step 3: Run the Plan 1 tests with the new names**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q`
Expected: `145 passed` (only the Plan 1 tests exist at this point).

- [ ] **Step 4: Apply the renames to the spec and to the Plan 1 text**

```bash
sed -i -e 's/surface_glycoprotein/signal_peptide_protein/g' -e 's/adhesion_repeat/tandem_repeat_protein/g' \
  -e 's/adhesion_domain/wall_family_domain/g' -e 's/antigen_candidate_surface/serodiagnostic_marker_candidate/g' \
  -e 's/allergen_homolog_hit/iuis_allergen_similarity/g' -e 's/allergen_candidate/iuis_allergen_homolog/g' \
  docs/superpowers/specs/2026-10-04-orchestrator-design.md
sed -i 's/`antigen_candidate`/`cocci_specificity_rank_top15`/g; s/antigen_candidate/cocci_specificity_rank_top15/g' docs/superpowers/specs/2026-10-04-orchestrator-design.md
```

The spec's `surface_glycoprotein[v]` columns become `signal_peptide_protein[v]`. Plan 1's text keeps the old names; its code is patched in Step 2. Say this in the pull request.

- [ ] **Step 5: Commit**

```bash
git branch --show-current   # must print sorting-hat-modules
git add src tests docs/superpowers/specs/2026-10-04-orchestrator-design.md
git commit -m "feat(sorting-hat): rename calls to say what the tools measure; report limits; measure keys

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 1: Module writer and the SignalP wrapper

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/__init__.py`, `src/cellsurface_sorting_hat/modules/base.py`, `src/cellsurface_sorting_hat/modules/signalp.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_base_signalp.py`

**Interfaces:**
- Consumes: `write_atomic` (Plan 1 cache), `NA_INVALID`, `Protein` (Plan 1 fasta).
- Produces: `sha256_file`, `params_hash`, `artefact_hash`, `ModuleSpec(name, version, params, artefacts, tools, artefact_digest)`, `invalid_row(protein)`, `write_module(workdir, spec, columns, rows, run_state='ok', note='') -> run record dict`.
- Produces: `parse_signalp(path) -> {id: {prediction, sp_prob}}` (refuses a header that is not `Organism: Eukarya`), `signalp_rows(proteins, parsed)`, `COLUMNS`, `SignalPFormatError`. Only prediction `SP` is a call (rule R0).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_base_signalp.py` with exactly this content:

```python
import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.base import (
    ModuleSpec,
    artefact_hash,
    params_hash,
    sha256_file,
    write_module,
)
from cellsurface_sorting_hat.modules.signalp import (
    COLUMNS,
    SignalPFormatError,
    parse_signalp,
    signalp_rows,
)

SIGNALP = (
    "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 20260930100126\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
    "CIMG_00013-t26_1-p1 | transcript=CIMG_00013-t26_1 | gene=CIMG_00013 | organism=Coccidioides_immitis_RS"
    "\tOTHER\t1.000000\t0.000000\t\n"
    "CIMG_04613-t26_1-p1 | gene=CIMG_04613 | protein_length=324\tSP\t0.000100\t0.999800\tCS pos: 19-20. Pr: 0.9\n"
    "XP_3\tLIPO\t0.1\t0.2\tCS pos: 20-21. Pr: 0.8\n"
)


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


def test_parse_signalp_uses_the_first_token_of_the_header_as_id(tmp_path):
    path = tmp_path / "prediction_results.txt"
    path.write_text(SIGNALP)
    got = parse_signalp(path)
    assert set(got) == {"CIMG_00013-t26_1-p1", "CIMG_04613-t26_1-p1", "XP_3"}
    assert got["CIMG_04613-t26_1-p1"] == {"prediction": "SP", "sp_prob": 0.9998}


def test_only_sp_is_a_call_lipo_is_not(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(SIGNALP)
    rows = signalp_rows(
        [prot("CIMG_04613-t26_1-p1"), prot("XP_3"), prot("CIMG_00013-t26_1-p1")],
        parse_signalp(path),
    )
    assert [r["call"] for r in rows] == ["called", "not_called", "not_called"]
    assert rows[1]["prediction"] == "LIPO"


def test_missing_and_invalid_proteins_are_not_ok(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(SIGNALP)
    rows = signalp_rows([prot("ABSENT"), prot("BAD", state="na_invalid")], parse_signalp(path))
    assert [(r["id"], r["state"]) for r in rows] == [("ABSENT", "error"), ("BAD", "na_invalid")]


@pytest.mark.parametrize(
    "text,message",
    [
        ("", "no prediction rows"),
        ("A\tSP\t0.1\n", "at least 4"),
        ("A\tSP\t0.1\tx\n", "not a number"),
        ("A\tSP\t0.1\t0.9\nA\tSP\t0.1\t0.9\n", "duplicate ID"),
    ],
)
def test_bad_signalp_files_are_refused(tmp_path, text, message):
    path = tmp_path / "p.txt"
    path.write_text(text)
    with pytest.raises(SignalPFormatError, match=message):
        parse_signalp(path)


def test_hashes(tmp_path):
    a, b = tmp_path / "a.bin", tmp_path / "b.bin"
    a.write_bytes(b"1")
    b.write_bytes(b"2")
    assert sha256_file(a) != sha256_file(b)
    assert params_hash({"x": 1, "y": 2}) == params_hash({"y": 2, "x": 1})
    assert params_hash({"x": 1}) != params_hash({"x": 2})
    assert artefact_hash([a, b]) == artefact_hash([b, a])  # order does not matter
    assert artefact_hash([a]) != artefact_hash([b])


def test_write_module_writes_the_table_and_the_run_record(tmp_path):
    db = tmp_path / "db.hmm"
    db.write_bytes(b"hmm")
    spec = ModuleSpec("step1_rule@R0", "1", {"mode": "fast"}, (db,), {"signalp": "6.0h"})
    rows = [
        {"id": "A", "state": "ok", "call": "called", "sp_prob": "0.9", "prediction": "SP"},
        {"id": "B", "state": "na_invalid"},
    ]
    rec = write_module(tmp_path, spec, COLUMNS, rows)
    with gzip.open(tmp_path / "modules" / "step1_rule@R0.tsv.gz", "rt") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert table[0]["call"] == "called" and table[1]["call"] == ""
    on_disk = json.loads((tmp_path / "modules" / "step1_rule@R0.json").read_text())
    assert on_disk == rec
    assert rec["params_hash"] == params_hash({"mode": "fast"})
    assert rec["artefact_hash"] == artefact_hash([db]) and rec["n_rows"] == 2


def test_write_module_refuses_a_duplicate_id(tmp_path):
    with pytest.raises(ValueError, match="duplicate row"):
        write_module(
            tmp_path,
            ModuleSpec("m", "1"),
            [],
            [{"id": "A", "state": "ok"}, {"id": "A", "state": "ok"}],
        )
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_base_signalp.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.modules`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/__init__.py` with exactly this content:

```python
"""Module wrappers: turn tool outputs into the module tables that the engine reads."""
```

Create `src/cellsurface_sorting_hat/modules/base.py` with exactly this content:

```python
"""Shared helpers: hashes and the writer for a module table plus its run record."""

import csv
import gzip
import hashlib
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.fasta import NA_INVALID


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def params_hash(params):
    """Hash of a parameter dict (canonical JSON)."""
    return hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()


def artefact_hash(paths):
    """Hash of the databases or models a module uses: sorted ``name:sha256`` lines."""
    lines = sorted(f"{Path(p).name}:{sha256_file(p)}" for p in paths)
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    version: str
    params: dict = field(default_factory=dict)
    artefacts: tuple = ()
    tools: dict = field(default_factory=dict)
    artefact_digest: str = (
        ""  # use instead of hashing large files (for example a recorded Pfam sha256)
    )


def invalid_row(protein):
    return {"id": protein.id, "state": NA_INVALID}


def write_module(workdir, spec, columns, rows, run_state="ok", note=""):
    """Write ``<workdir>/modules/<name>.tsv.gz`` and ``<name>.json`` (see the module contract).

    ``columns`` are the fields after ``id`` and ``state``. A row that lacks a column gets an empty
    value. Every row needs ``id`` and ``state``.
    """
    folder = Path(workdir) / "modules"
    header = ["id", "state"] + list(columns)
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    seen = set()
    for r in rows:
        if r["id"] in seen:
            raise ValueError(f"{spec.name}: duplicate row for ID {r['id']!r}")
        seen.add(r["id"])
        writer.writerow([r.get(h, "") for h in header])
    write_atomic(folder / f"{spec.name}.tsv.gz", gzip.compress(buf.getvalue().encode(), mtime=0))
    record = {
        "module": spec.name,
        "version": spec.version,
        "params_hash": params_hash(spec.params),
        "artefact_hash": spec.artefact_digest or artefact_hash(spec.artefacts),
        "run_state": run_state,
        "params": spec.params,
        "tools": spec.tools,
        "n_rows": len(seen),
        "note": note,
    }
    write_atomic(
        folder / f"{spec.name}.json", (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    )
    return record
```

Create `src/cellsurface_sorting_hat/modules/signalp.py` with exactly this content:

```python
"""SignalP 6 results -> module ``step1_rule@R0`` (rule R0: SignalP calls a signal peptide, ``SP``)."""

from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

MODULE = "step1_rule@R0"


class SignalPFormatError(ValueError):
    """The SignalP file does not have the expected shape."""


def parse_signalp(path):
    """Return ``{protein ID: {"prediction": str, "sp_prob": float}}``.

    ``prediction_results.txt`` has two ``#`` header lines, then one TAB separated row per protein:
    ID (the FASTA header, the ID is its first token), prediction, OTHER probability, SP(Sec/SPI)
    probability, cleavage site. Prediction ``SP`` is the Sec/SPI signal peptide that rule R0 uses.
    """
    out = {}
    for n, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        if line.startswith("#"):
            if "Organism:" in line and "Eukarya" not in line:
                raise SignalPFormatError(f"{path}:{n}: not a Eukarya run: {line.strip()!r}")
            continue
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise SignalPFormatError(f"{path}:{n}: expected at least 4 TAB separated fields")
        pid = fields[0].split()[0]
        try:
            sp_prob = float(fields[3])
        except ValueError:
            raise SignalPFormatError(f"{path}:{n}: SP probability is not a number") from None
        if pid in out:
            raise SignalPFormatError(f"{path}:{n}: duplicate ID {pid!r}")
        out[pid] = {"prediction": fields[1], "sp_prob": sp_prob}
    if not out:
        raise SignalPFormatError(f"{path}: no prediction rows")
    return out


def signalp_rows(proteins, parsed):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            hit = parsed[p.id]
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if hit["prediction"] == "SP" else "not_called",
                    "sp_prob": f"{hit['sp_prob']:.6f}",
                    "prediction": hit["prediction"],
                }
            )
    return rows


COLUMNS = ["call", "sp_prob", "prediction"]
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_base_signalp.py -q`
Expected: `10 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules tests/cellsurface_sorting_hat/modules/test_base_signalp.py
git commit -m "feat(sorting-hat): module writer and SignalP R0 wrapper

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 2: Pfam family table, domain table parser and the specificity report

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/pfam.py`, `data/sorting_hat/family_table.tsv`
- Test: `tests/cellsurface_sorting_hat/modules/test_pfam.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `load_family_table(path) -> list[Family]`, `parse_domtblout(path)` (refuses a table without the `# [ok]` trailer or without `--cut_ga`), `check_hit_ids(hits, fasta_ids)`, `pfam_rows(proteins, hits, families, module, sp_calls=None, tm_counts=None)` (raises `NoActiveFamilyError` when no family of the module is active), `specificity_report(hit_ids, member_ids, universe_ids)`, `FamilyTableError`, `FAMILY_COLUMNS`, `MODULES`, `COLUMNS`.
- A family counts only if `active = yes` (with `active_by` and `active_date`). `second_condition` is empty, `signal_peptide` or `no_tm` (no helix that starts after residue 35).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_pfam.py` with exactly this content:

```python
import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.pfam import (
    FAMILY_COLUMNS,
    FamilyTableError,
    NoActiveFamilyError,
    check_hit_ids,
    load_family_table,
    parse_domtblout,
    pfam_rows,
    specificity_report,
)

OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"
TRAILER = "# [ok]\n"
DOMTBL = (
    OPTIONS
    + (
        "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description\n"
        "P1 - 289 CFEM PF05730.17 70 1e-20 60.1 8.9 1 1 1e-21 2e-20 59.0 8.9 1 70 20 90 20 91 0.9 -\n"
        "P2 - 400 Hydrophobin PF01185.24 60 1e-12 40.0 0.0 1 1 1e-13 3e-12 39.0 0.0 1 60 5 65 5 66 0.9 -\n"
        "P3 - 300 AltA1 PF16541.11 150 1e-30 90.0 0.0 1 1 1e-31 1e-30 89.0 0.0 1 150 10 160 10 161 0.9 -\n"
        "P4 - 500 Asp PF00026.29 300 1e-40 120.0 0.0 1 1 1e-41 1e-40 119.0 0.0 1 300 10 310 10 311 0.9 -\n"
    )
    + TRAILER
)


def write_table(path, rows):
    lines = ["\t".join(FAMILY_COLUMNS)]
    for r in rows:
        lines.append("\t".join(r.get(c, "") for c in FAMILY_COLUMNS))
    path.write_text("\n".join(lines) + "\n")


def fam(acc, module="pfam_adhesion", active="yes", second="", **extra):
    base = {
        "pfam_acc": acc,
        "name": acc,
        "class": "x",
        "module": module,
        "source_pmid": "1",
        "pfam_release": "38.2",
        "specificity_note": "n",
        "second_condition": second,
        "active": active,
        "active_by": "owner" if active == "yes" else "",
        "active_date": "2026-10-05" if active == "yes" else "",
    }
    base.update(extra)
    return base


def prot(pid):
    return Protein(pid, "MKT", "x" * 64, "ok", "", 0.0)


def test_domtblout_is_parsed_without_the_accession_version(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL)
    hits = parse_domtblout(path)
    assert [(h["target"], h["acc"]) for h in hits] == [
        ("P1", "PF05730"),
        ("P2", "PF01185"),
        ("P3", "PF16541"),
        ("P4", "PF00026"),
    ]
    assert hits[0]["ievalue"] == 2e-20


def test_only_active_families_of_the_module_count(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(
        path, [fam("PF05730"), fam("PF01185", active="no"), fam("PF16541", module="pfam_allergen")]
    )
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    prots = [prot(p) for p in ("P1", "P2", "P3", "P4")]
    adh = {r["id"]: r["hit"] for r in pfam_rows(prots, hits, fams, "pfam_adhesion")}
    all_ = {r["id"]: r["hit"] for r in pfam_rows(prots, hits, fams, "pfam_allergen")}
    assert adh == {"P1": "1", "P2": "0", "P3": "0", "P4": "0"}  # PF01185 is inactive
    assert all_ == {"P1": "0", "P2": "0", "P3": "1", "P4": "0"}


def test_second_condition_needs_a_signal_peptide_call(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF00026", second="signal_peptide")])
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    p4 = [prot("P4")]
    assert pfam_rows(p4, hits, fams, "pfam_adhesion")[0]["hit"] == "0"  # no SP information
    assert pfam_rows(p4, hits, fams, "pfam_adhesion", {"P4": "not_called"})[0]["hit"] == "0"
    assert pfam_rows(p4, hits, fams, "pfam_adhesion", {"P4": "called"})[0]["hit"] == "1"


@pytest.mark.parametrize(
    "mutate,message",
    [
        (lambda r: r.update(pfam_acc="PF05730.17"), "no version"),
        (lambda r: r.update(module="other"), "module must be"),
        (lambda r: r.update(active="maybe"), "active must be"),
        (lambda r: r.update(second_condition="tm"), "second_condition"),
        (lambda r: r.update(active="yes", active_by=""), "needs active_by"),
    ],
)
def test_bad_family_rows_are_refused(tmp_path, mutate, message):
    row = fam("PF05730")
    mutate(row)
    path = tmp_path / "f.tsv"
    write_table(path, [row])
    with pytest.raises(FamilyTableError, match=message):
        load_family_table(path)


def test_duplicate_families_are_refused(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730"), fam("PF05730")])
    with pytest.raises(FamilyTableError, match="duplicate"):
        load_family_table(path)


def test_specificity_report_counts_and_lists_the_non_member_hits():
    rep = specificity_report(
        hit_ids={"a", "b", "x"}, member_ids={"a", "b", "c"}, universe_ids=set("abcxyz")
    )
    assert (rep["tp"], rep["fp"], rep["fn"], rep["tn"]) == (2, 1, 1, 2)
    assert rep["nonmember_hits"] == ["x"] and rep["missed_members"] == ["c"]
    assert rep["sensitivity"] == pytest.approx(2 / 3) and rep["specificity"] == pytest.approx(2 / 3)


def test_specificity_report_needs_members_inside_the_universe():
    with pytest.raises(ValueError):
        specificity_report({"a"}, {"zz"}, {"a"})


def test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730", second="no_tm")])
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    p1 = [prot("P1")]
    assert pfam_rows(p1, hits, fams, "pfam_adhesion")[0]["hit"] == "0"  # no TM information
    assert (
        pfam_rows(p1, hits, fams, "pfam_adhesion", tm_counts={"P1": 7})[0]["hit"] == "0"
    )  # a receptor
    assert pfam_rows(p1, hits, fams, "pfam_adhesion", tm_counts={"P1": 0})[0]["hit"] == "1"


def test_a_domain_table_without_the_ok_trailer_or_cut_ga_is_refused(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL.replace(TRAILER, ""))
    with pytest.raises(ValueError, match="no '# \\[ok\\]' trailer"):
        parse_domtblout(path)
    path.write_text(DOMTBL.replace("--cut_ga", "-E 1e-5"))
    with pytest.raises(ValueError, match="not run with --cut_ga"):
        parse_domtblout(path)


def test_domain_table_targets_that_are_not_in_the_fasta_are_refused(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL)
    hits = parse_domtblout(path)
    with pytest.raises(ValueError, match="not in the FASTA"):
        check_hit_ids(hits, ["P1", "P2"])  # P3 and P4 are missing
    check_hit_ids(hits, ["P1", "P2", "P3", "P4", "P5"])  # no error


def test_a_module_with_no_active_family_raises_instead_of_writing_zeros(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730", active="no")])
    with pytest.raises(NoActiveFamilyError):
        pfam_rows([prot("P1")], [], load_family_table(path), "pfam_adhesion")
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_pfam.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.modules.pfam`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/pfam.py` with exactly this content:

```python
"""Pfam family table and ``hmmsearch --domtblout`` results -> modules ``pfam_adhesion``, ``pfam_allergen``."""

import csv
from dataclasses import dataclass
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

FAMILY_COLUMNS = [
    "pfam_acc",
    "name",
    "class",
    "module",
    "source_pmid",
    "pfam_release",
    "specificity_note",
    "second_condition",
    "active",
    "active_by",
    "active_date",
]
MODULES = ("pfam_adhesion", "pfam_allergen")
SECOND_CONDITIONS = ("", "signal_peptide", "no_tm")
COLUMNS = ["hit", "families"]


class FamilyTableError(ValueError):
    """The family table is not valid."""


class NoActiveFamilyError(ValueError):
    """No family of the module has passed its specificity test and sign-off."""


@dataclass(frozen=True)
class Family:
    pfam_acc: str
    name: str
    cls: str
    module: str
    second_condition: str
    active: bool


def load_family_table(path):
    rows, seen = [], set()
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(FAMILY_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise FamilyTableError(f"{path}: missing column(s) {sorted(missing)}")
        for n, r in enumerate(reader, 2):
            where = f"{path}:{n}"
            acc = r["pfam_acc"].strip()
            if not acc.startswith("PF") or "." in acc:
                raise FamilyTableError(f"{where}: pfam_acc must look like PF05730 (no version)")
            if acc in seen:
                raise FamilyTableError(f"{where}: duplicate {acc}")
            seen.add(acc)
            if r["module"] not in MODULES:
                raise FamilyTableError(f"{where}: module must be one of {MODULES}")
            if r["second_condition"] not in SECOND_CONDITIONS:
                raise FamilyTableError(
                    f"{where}: second_condition must be one of {SECOND_CONDITIONS}"
                )
            if r["active"] not in ("yes", "no"):
                raise FamilyTableError(f"{where}: active must be yes or no")
            if r["active"] == "yes" and not (r["active_by"].strip() and r["active_date"].strip()):
                raise FamilyTableError(f"{where}: an active family needs active_by and active_date")
            rows.append(
                Family(
                    acc,
                    r["name"],
                    r["class"],
                    r["module"],
                    r["second_condition"],
                    r["active"] == "yes",
                )
            )
    return rows


def parse_domtblout(path):
    """Return a list of ``{"target", "acc", "ievalue", "score"}`` (accession without version).

    hmmsearch writes ``# [ok]`` as the last line of a finished run and ``--cut_ga`` in its option
    settings. A file without the trailer is truncated; a file without ``--cut_ga`` used another
    cutoff. Both are refused.
    """
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    if not any(line.strip() == "# [ok]" for line in lines[-3:]):
        raise ValueError(f"{path}: no '# [ok]' trailer; the hmmsearch run is not finished")
    if not any("--cut_ga" in line for line in lines if line.startswith("#")):
        raise ValueError(f"{path}: hmmsearch was not run with --cut_ga")
    hits = []
    for n, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 22:
            raise ValueError(f"{path}:{n}: expected at least 22 fields, found {len(f)}")
        hits.append(
            {
                "target": f[0],
                "acc": f[4].split(".")[0],
                "ievalue": float(f[12]),
                "score": float(f[13]),
            }
        )
    return hits


def check_hit_ids(hits, fasta_ids):
    """A domain table whose targets are not FASTA IDs belongs to another proteome."""
    unknown = sorted({h["target"] for h in hits} - set(fasta_ids))
    if unknown:
        raise ValueError(
            f"{len(unknown)} domain-table target(s) are not in the FASTA, for example {unknown[0]!r}"
        )


def pfam_rows(proteins, hits, families, module, sp_calls=None, tm_counts=None):
    """Rows for ``module``. A protein is a hit if it has a domain of an active family of this module.

    ``second_condition = signal_peptide`` counts a hit only when ``sp_calls[id]`` is ``called``.
    ``second_condition = no_tm`` counts a hit only when ``tm_counts[id]`` is 0 (for example a CFEM
    domain in a receptor with transmembrane helices is not a cell wall CFEM protein). If the needed
    table is None, such a family never counts.
    """
    wanted = {f.pfam_acc: f for f in families if f.module == module and f.active}
    if not wanted:
        raise NoActiveFamilyError(f"{module}: no family of this module is active")
    found = {}
    for h in hits:
        fam = wanted.get(h["acc"])
        if fam is None:
            continue
        if (
            fam.second_condition == "signal_peptide"
            and (sp_calls or {}).get(h["target"]) != "called"
        ):
            continue
        if fam.second_condition == "no_tm" and (tm_counts or {}).get(h["target"]) != 0:
            continue
        found.setdefault(h["target"], set()).add(fam.pfam_acc)
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        accs = sorted(found.get(p.id, ()))
        rows.append(
            {"id": p.id, "state": "ok", "hit": "1" if accs else "0", "families": ",".join(accs)}
        )
    return rows


def specificity_report(hit_ids, member_ids, universe_ids):
    """Compare the proteins hit by a family with the known members (a specificity test).

    ``universe_ids`` are all proteins searched. Returns counts and the non-member hits, which a
    person must read before the family is made active.
    """
    hit, members, universe = set(hit_ids), set(member_ids), set(universe_ids)
    if not members <= universe:
        raise ValueError("member_ids must be a subset of universe_ids")
    tp, fp, fn = hit & members, hit - members, members - hit
    tn = universe - hit - members
    return {
        "tp": len(tp),
        "fp": len(fp),
        "fn": len(fn),
        "tn": len(tn),
        "sensitivity": len(tp) / len(members) if members else None,
        "specificity": len(tn) / (len(tn) + len(fp)) if (len(tn) + len(fp)) else None,
        "nonmember_hits": sorted(fp),
        "missed_members": sorted(fn),
    }
```

The family table is TAB separated (the first line is the header). All fifteen families start inactive. A family is made active only after its specificity review and a sign-off (Task 13). The model names (column 2) are checked against the Pfam database by the job script (Task 10) and were checked against Pfam 38.2 on 2026-10-05.

Create `data/sorting_hat/family_table.tsv` with exactly this content:

```text
pfam_acc	name	class	module	source_pmid	pfam_release	specificity_note	second_condition	active	active_by	active_date
PF05730	CFEM	2b-i CFEM-domain surface protein	pfam_adhesion	docs/reports/2026-09-29-class2b-structure.md; PMID 28513415 names CFEM proteins	38.2	hemophore fold, not shown to bind a host receptor; Pth11-like GPCRs also carry a CFEM domain, so no_tm keeps only proteins without TMHMM helices; specificity test not run	no_tm	no		
PF04681	Bys1	2b-iii Bys1 invasin (CalA)	pfam_adhesion	docs/reports/2026-09-29-class2b-structure.md	38.2	Bys1 is not a thaumatin Pfam; calB and calC paralogs are the control; specificity test not run		no		
PF01185	Hydrophobin	2c hydrophobin	pfam_adhesion	PMID 28513415 names RodA and RodB	38.2	specificity test not run		no		
PF06766	Hydrophobin_2	2c hydrophobin	pfam_adhesion	PMID 28513415 names RodA and RodB	38.2	specificity test not run		no		
PF28987	DewD	2c hydrophobin (DewD)	pfam_adhesion	analysis/cys_candidates/01_known_family_hmm.sh	38.2	specificity test not run		no		
PF22354	Eas	2c hydrophobin (Eas)	pfam_adhesion	Pfam 38.2 name Eas, Hydrophobin	38.2	not in the earlier HMM script; specificity test not run		no		
PF11766	Candida_ALS_N	2a Als adhesin N-terminal domain	pfam_adhesion	PMID 28513415 names the Als family	38.2	Candida-specific; specificity test not run		no		
PF05792	Candida_ALS	2a Als adhesin repeat region	pfam_adhesion	PMID 28513415 names the Als family	38.2	Candida-specific; specificity test not run		no		
PF16541	AltA1	Alt a 1 allergen family	pfam_allergen	WHO/IUIS Alt a 1; docs/reports/2026-10-04-fungal-allergen-scoping.md	38.2	allergen-specific; specificity test not run		no		
PF25312	Allergen_Asp_f_4	Asp f 4 allergen family	pfam_allergen	WHO/IUIS Asp f 4; docs/reports/2026-10-04-fungal-allergen-scoping.md	38.2	allergen-specific; specificity test not run		no		
PF10182	Flo11	2a yeast flocculin domain	pfam_adhesion	Pfam 38.2 Flo11; yeast adhesins (FLO11)	38.2	S. cerevisiae S288C carries a Flo8 defect and does not flocculate; specificity test not run		no		
PF10528	GLEYA	2a yeast adhesin lectin domain	pfam_adhesion	Pfam 38.2 GLEYA (Epa/Flo adhesins)	38.2	specificity test not run		no		
PF11765	Hyphal_reg_CWP	2a Hwp1-like wall protein	pfam_adhesion	Pfam 38.2; PMID 28513415 names Hwp1	38.2	Candida-specific; specificity test not run		no		
PF00399	PIR	covalently bound yeast wall protein repeat (Pir)	pfam_adhesion	Pfam 38.2 PIR; Pir proteins are structural wall proteins, not adhesins	38.2	structural wall family; counted as a wall family only; specificity test not run		no		
PF07691	PA14	2a PA14 adhesin domain (Flo5, Epa)	pfam_adhesion	Pfam 38.2 PA14	38.2	broad family with bacterial members and beta-glucosidases; needs a signal peptide	signal_peptide	no
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_pfam.py -q`
Expected: `15 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/pfam.py data/sorting_hat/family_table.tsv tests/cellsurface_sorting_hat/modules/test_pfam.py
git commit -m "feat(sorting-hat): Pfam family table and domain table parser

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 3: Repeat detector wrappers

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/repeats.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_repeats.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `parse_repeat_table(path)`, `repeat_rows(proteins, parsed, min_coverage=0.25, min_copies=2.5, min_len=80)`, `run_detector(...)`, `COLUMNS`, `DETECTOR_MIN_LEN`.
- The call is `period > 0 and coverage >= min_coverage and copies >= min_copies` (the rule of `03_repeat_surface_candidates.py`, without its secretion gate: the gate is applied by `cell_wall_adhesion_candidate`). Both detector scripts skip proteins shorter than 80 aa; such a protein is `ok` with `not_called`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_repeats.py` with exactly this content:

```python
import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.repeats import parse_repeat_table, repeat_rows, run_detector


def prot(pid, state="ok", length=100):
    return Protein(pid, "M" * length, "x" * 64, state, "", 0.0)


REPEATS = (
    "strain\tprotein\tlength\trep_period\trep_score\trep_start\trep_end\trep_n_copies\trep_coverage\n"
    "S\tA\t300\t47\t0.8\t80\t270\t4.0\t0.62\n"
    "S\tB\t300\t0\t0.1\t0\t0\t0\t0.0\n"
    "S\tC\t300\t12\t0.5\t10\t60\t4.1\t0.17\n"
    "S\tD\t300\t12\t0.5\t10\t60\t2.4\t0.40\n"
)


def test_repeat_call_needs_period_coverage_and_copies(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows([prot(x) for x in "ABCD"] + [prot("E")], parse_repeat_table(path))
    assert [r.get("call") for r in rows] == [
        "called",
        "not_called",
        "not_called",
        "not_called",
        None,
    ]
    assert rows[-1]["state"] == "error"  # not in the detector output


def test_repeat_thresholds_are_parameters(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows([prot("C")], parse_repeat_table(path), min_coverage=0.10)
    assert rows[0]["call"] == "called"


def test_repeat_table_needs_its_columns_and_unique_proteins(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text("protein\trep_period\nA\t3\n")
    with pytest.raises(ValueError, match="missing column"):
        parse_repeat_table(path)
    path.write_text("protein\trep_period\trep_n_copies\trep_coverage\nA\t3\t3\t0.5\nA\t3\t3\t0.5\n")
    with pytest.raises(ValueError, match="duplicate protein"):
        parse_repeat_table(path)


def test_run_detector_builds_the_command(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr("subprocess.run", lambda cmd, check: calls.append((cmd, check)))
    cmd = run_detector(
        tmp_path, "repeat14", "in.faa", "out.tsv", python="py", extra=["--mode", "exact"]
    )
    assert calls == [(cmd, True)]
    assert cmd[0] == "py" and cmd[1].endswith("analysis/cocci_repeats/14_repeat_detect_general.py")
    assert cmd[2:] == ["in.faa", "--out", "out.tsv", "--mode", "exact"]


def test_a_protein_shorter_than_the_detector_minimum_is_not_called_not_an_error(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows(
        [prot("SHORT", length=50), prot("LONG_MISSING", length=200)], parse_repeat_table(path)
    )
    assert (rows[0]["state"], rows[0]["call"]) == ("ok", "not_called")
    assert rows[1]["state"] == "error"  # long enough to be profiled, but absent from the table
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_repeats.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.repeats`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/repeats.py` with exactly this content:

```python
"""Repeat detector tables -> modules ``repeat02`` and ``repeat14`` (call = repeat per 03_repeat_surface_candidates.py)."""

import csv
import subprocess
import sys
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

DEFAULT_MIN_COVERAGE = 0.25
DEFAULT_MIN_COPIES = 2.5
DETECTOR_MIN_LEN = 80  # the --min-len default of both detector scripts
SCRIPTS = {"repeat02": "02_repeat_profile.py", "repeat14": "14_repeat_detect_general.py"}


def parse_repeat_table(path):
    """Return ``{protein ID: {"period": int, "copies": float, "coverage": float}}``."""
    out = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"protein", "rep_period", "rep_n_copies", "rep_coverage"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
        for n, r in enumerate(reader, 2):
            if r["protein"] in out:
                raise ValueError(f"{path}:{n}: duplicate protein {r['protein']!r}")
            out[r["protein"]] = {
                "period": int(float(r["rep_period"] or 0)),
                "copies": float(r["rep_n_copies"] or 0),
                "coverage": float(r["rep_coverage"] or 0),
            }
    return out


def repeat_rows(
    proteins,
    parsed,
    min_coverage=DEFAULT_MIN_COVERAGE,
    min_copies=DEFAULT_MIN_COPIES,
    min_len=DETECTOR_MIN_LEN,
):
    """The detectors skip proteins shorter than ``min_len`` (default 80). Such a protein cannot hold
    a repeat array that the detectors can see, so it is ``not_called`` with ``period`` 0."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed and len(p.sequence) < min_len:
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "not_called",
                    "period": 0,
                    "copies": "0.00",
                    "coverage": "0.000",
                }
            )
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            r = parsed[p.id]
            called = r["period"] > 0 and r["coverage"] >= min_coverage and r["copies"] >= min_copies
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if called else "not_called",
                    "period": r["period"],
                    "copies": f"{r['copies']:.2f}",
                    "coverage": f"{r['coverage']:.3f}",
                }
            )
    return rows


COLUMNS = ["call", "period", "copies", "coverage"]


def run_detector(repo_root, module, fasta, out_tsv, python=sys.executable, extra=()):
    """Run an existing detector script from ``analysis/cocci_repeats`` on one FASTA."""
    script = Path(repo_root) / "analysis" / "cocci_repeats" / SCRIPTS[module]
    cmd = [python, str(script), str(fasta), "--out", str(out_tsv), *extra]
    subprocess.run(cmd, check=True)
    return cmd
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_repeats.py -q`
Expected: `5 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/repeats.py tests/cellsurface_sorting_hat/modules/test_repeats.py
git commit -m "feat(sorting-hat): repeat detector wrappers

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 4: Allergen homology and the leave-species-out recall

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/allergen.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_allergen.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `build_allergen_fasta(isoallergen_tsv, out_fasta, allergen_tsv=None) -> (n_written, skipped)` (removes whitespace, refuses free text and fragments under 30 residues, writes `<out_fasta>.meta.tsv` with species, order, exposure and the IUIS evidence text), `read_meta`, `parse_blast`, `check_blast_ids`, `allergen_rows(proteins, best, meta=None)` (fields `identity`, `aligned_length`, `coverage`, `allergen_name`, `allergen_species`, `exposure`, `evidence`), `species_code`, `parse_blast_hits`, `best_hit_other_species`, `recall_by_rule`, `lso_report(blast_path, fasta_ids, rules=None)`, `BLAST_FIELDS`, `COLUMNS`.
- The module reports one local alignment (the best bit score). The engine's `iuis_allergen_similarity` (identity >= 35% and aligned length >= 80) is not the sliding 80-residue window of the Codex rule; a window count measured inside the BLAST alignments was 181 against 105 Af293 proteins (bioinformatics review, lower bound).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_allergen.py` with exactly this content:

```python
import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.allergen import (
    allergen_name,
    allergen_rows,
    best_hit_other_species,
    build_allergen_fasta,
    check_blast_ids,
    lso_report,
    parse_blast,
    parse_blast_hits,
    read_meta,
    recall_by_rule,
    species_code,
)


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


SEQ1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"  # 33 residues
SEQ2 = "MNLLPQWERTYIPASDFGHKLCVNMQRSTWYAAA"  # 34 residues
IUIS = (
    "AllergenID\tIsoName\tName\tSequence\n"
    f"11\tAsp f 1.0101\tAsp f 1\t{SEQ1[:20]} {SEQ1[20:]}\n"
    f"12\tAsp f 2.0101\tAsp f 2\t{SEQ2}\n"
    "13\tNoSeq\tNo seq\t\n"
    "14\tEpi p 1.0101\tEpi p 1\tN-TERMINAL PEPTIDE: ADGIVAVELDTY >INTERNAL PEPTIDE: RGSFXK\n"
    "15\tPep 1.0101\tPep 1\tADGIVAVELDTY\n"
)
ALLERGENS = (
    "AllergenID\tName\tSpecies\tTaxOrder\tExposure\tAllergenicity\n"
    "11\tAsp f 1\tAspergillus fumigatus\tEurotiales\tAirway\tIgE binding in 75% of 40 sera\n"
    "12\tAsp f 2\tAspergillus fumigatus\tEurotiales\tAirway\t\n"
)


def test_allergen_fasta_is_clean_and_has_a_meta_table(tmp_path):
    src, al, out = tmp_path / "i.tsv", tmp_path / "a.tsv", tmp_path / "a.faa"
    src.write_text(IUIS)
    al.write_text(ALLERGENS)
    n, skipped = build_allergen_fasta(src, out, al)
    assert n == 2
    assert (
        out.read_text() == f">Asp_f_1.0101|11\n{SEQ1}\n>Asp_f_2.0101|12\n{SEQ2}\n"
    )  # spaces removed
    assert skipped == [
        ("Epi_p_1.0101|14", "not a protein sequence"),
        ("Pep_1.0101|15", "peptide fragment of 12 residues"),
    ]
    meta = read_meta(str(out) + ".meta.tsv")
    assert (
        meta["Asp_f_1.0101|11"]["exposure"] == "Airway"
        and "75% of 40 sera" in meta["Asp_f_1.0101|11"]["evidence"]
    )
    assert (
        meta["Asp_f_2.0101|12"]["evidence"] == ""
        and meta["Asp_f_2.0101|12"]["species"] == "Aspergillus fumigatus"
    )


def test_duplicate_allergen_ids_are_refused(tmp_path):
    src = tmp_path / "i.tsv"
    src.write_text(
        f"AllergenID\tIsoName\tName\tSequence\n1\tA 1\tA 1\t{SEQ1}\n1\tA 1\tA 1\t{SEQ2}\n"
    )
    with pytest.raises(ValueError, match="duplicate allergen ID"):
        build_allergen_fasta(src, tmp_path / "a.faa")


BLAST = (
    "P1\tAsp_f_1.0101|11\t99.0\t100\t120\t125\t200\t1e-50\n"
    "P1\tAsp_f_2.0101|12\t45.0\t90\t120\t300\t60\t1e-5\n"
    "P2\tAsp_f_2.0101|12\t38.0\t85\t200\t300\t50\t1e-3\n"
)


def test_blast_best_hit_is_by_bit_score_and_fields_are_derived(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text(BLAST)
    best = parse_blast(path)
    meta = {
        "Asp_f_1.0101|11": {
            "species": "Aspergillus fumigatus",
            "exposure": "Airway",
            "evidence": "75% of 40",
        }
    }
    rows = {
        r["id"]: r
        for r in allergen_rows(
            [prot("P1"), prot("P2"), prot("P3"), prot("BAD", "na_invalid")], best, meta
        )
    }
    assert rows["P1"]["identity"] == "99.00" and rows["P1"]["allergen_name"] == "Asp_f_1.0101"
    assert rows["P1"]["coverage"] == "80.0" and rows["P1"]["exposure"] == "Airway"
    assert (rows["P2"]["identity"], rows["P2"]["aligned_length"], rows["P2"]["coverage"]) == (
        "38.00",
        "85",
        "28.3",
    )
    assert (rows["P3"]["identity"], rows["P3"]["aligned_length"], rows["P3"]["allergen_name"]) == (
        "0",
        "0",
        "",
    )
    assert rows["BAD"]["state"] == "na_invalid"


def test_blast_line_with_the_wrong_field_count_is_refused(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text("P1\tS\t99.0\n")
    with pytest.raises(ValueError, match="expected 8 fields"):
        parse_blast(path)


def test_a_blast_table_from_another_proteome_is_refused():
    best = {"X9": {"subject": "S", "identity": 90.0, "length": 100, "qlen": 100, "slen": 100}}
    with pytest.raises(ValueError, match="not in the FASTA"):
        check_blast_ids(best, ["P1", "P2"])
    check_blast_ids({"P1": best["X9"]}, ["P1", "P2"])  # no error


def test_names_and_species_codes():
    assert allergen_name("Asp_f_1.0101|11") == "Asp_f_1.0101"
    assert species_code("Asp_f_1.0101|11") == "Asp_f"
    assert species_code("Cand_a_3.0101|5") == "Cand_a"


def h(q, s, ident, length, qlen=100, slen=100, bits=100.0):
    return {
        "query": q,
        "subject": s,
        "identity": ident,
        "length": length,
        "qlen": qlen,
        "slen": slen,
        "bitscore": bits,
    }


IDS = ["Asp_f_1.0101|1", "Asp_f_2.0101|2", "Asp_n_1.0101|3", "Alt_a_1.0101|4"]


def test_other_species_hits_exclude_the_same_species():
    hits = [
        h(IDS[0], IDS[1], 99, 100, bits=300),
        h(IDS[0], IDS[2], 60, 90, bits=100),
        h(IDS[0], IDS[3], 40, 85, bits=50),
    ]
    best = best_hit_other_species(hits, IDS)
    assert best[IDS[0]]["subject"] == IDS[2]  # the 99% hit is another allergen of the same species


def test_recall_by_rule_uses_the_rules_of_the_engine_and_the_fasta_denominator():
    best = {IDS[0]: h(IDS[0], IDS[2], 60, 90), IDS[1]: h(IDS[1], IDS[3], 72, 85)}
    rules = [("similarity", 35.0, 80, 0.0), ("homolog", 70.0, 0, 80.0)]
    got = {r["rule"]: (r["recovered"], r["n"]) for r in recall_by_rule(IDS, best, rules)}
    assert got == {
        "similarity": (2, 4),
        "homolog": (1, 4),
    }  # IDS[2], IDS[3] have no hit: counted as missed


def test_lso_report_refuses_a_sequence_with_no_blast_line(tmp_path):
    path = tmp_path / "b.tsv"
    lines = [f"{i}\t{i}\t100.0\t100\t100\t100\t500\t0" for i in IDS[:3]]  # IDS[3] absent
    path.write_text("\n".join(lines) + "\n")
    assert len(parse_blast_hits(path)) == 3
    with pytest.raises(ValueError, match="no BLAST line"):
        lso_report(path, IDS)
    lines.append(f"{IDS[3]}\t{IDS[3]}\t100.0\t100\t100\t100\t500\t0")
    path.write_text("\n".join(lines) + "\n")
    rep = lso_report(path, IDS)
    assert (rep["n_sequences"], rep["n_species"]) == (4, 3)
    assert [r["recovered"] for r in rep["recall"]] == [0, 0]
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_allergen.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.allergen`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/allergen.py` with exactly this content:

```python
"""WHO/IUIS fungal allergen sequences and BLASTP results -> module ``allergen_homology``.

Module fields: ``identity`` (percent identity of the best local alignment), ``aligned_length``
(alignment length, gaps included), ``coverage`` (alignment length as a percent of the allergen
sequence, capped at 100), ``allergen_name``, ``allergen_species``, ``exposure`` and ``evidence``
(the IUIS evidence text of the matched allergen, or empty). A protein with no hit has ``0``, ``0``,
``0`` and empty text (state ``ok``). The engine makes two calls from these fields:
``iuis_allergen_similarity`` (identity >= 35% and aligned length >= 80, ONE local alignment; this is
not the sliding 80-residue window of the Codex rule) and ``iuis_allergen_homolog`` (identity >= 70%
and coverage >= 80%).
"""

import csv
import re
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

BLAST_FIELDS = "qseqid sseqid pident length qlen slen bitscore evalue"
COLUMNS = [
    "identity",
    "aligned_length",
    "coverage",
    "allergen_name",
    "allergen_species",
    "exposure",
    "evidence",
]
META_COLUMNS = ["id", "name", "species", "tax_order", "exposure", "evidence", "length"]
RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYXBZUJO")
MIN_USABLE_LENGTH = (
    30  # shorter entries are peptide fragments; they cannot reach an 80 aa alignment
)
_ID_BAD = re.compile(r"[^A-Za-z0-9_.\-]+")


def build_allergen_fasta(isoallergen_tsv, out_fasta, allergen_tsv=None):
    """Write the IUIS fungal sequences as FASTA and ``<out_fasta>.meta.tsv``.

    ID = ``<IsoAllergenName>|<AllergenID>``. All whitespace is removed from a sequence. An entry
    whose sequence has characters that are not residues (free text such as "N-terminal peptide: ...")
    or is shorter than 30 residues is not written and is listed in the returned ``skipped``. The
    meta table carries species, order, exposure and the IUIS evidence text of each allergen (from
    ``allergen_tsv``, the ``iuis_fungal_allergens.tsv`` table, when given).

    Returns ``(n_written, skipped)`` where ``skipped`` is a list of ``(id, reason)``.
    """
    allergen = {}
    if allergen_tsv:
        with open(allergen_tsv, encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                allergen[str(r["AllergenID"])] = r
    n, skipped, seen = 0, [], set()
    meta_path = Path(str(out_fasta) + ".meta.tsv")
    with (
        open(isoallergen_tsv, encoding="utf-8-sig", newline="") as fh,
        open(out_fasta, "w") as out,
        open(meta_path, "w") as meta,
    ):
        meta.write("\t".join(META_COLUMNS) + "\n")
        for r in csv.DictReader(fh, delimiter="\t"):
            raw = (r.get("Sequence") or "").strip()
            name = _ID_BAD.sub("_", (r.get("IsoName") or r.get("Name") or "").strip())
            pid = f"{name}|{r['AllergenID']}"
            seq = "".join(raw.split()).upper()
            if not seq or seq == "NAN":
                continue
            if not set(seq) <= RESIDUES:
                skipped.append((pid, "not a protein sequence"))
                continue
            if len(seq) < MIN_USABLE_LENGTH:
                skipped.append((pid, f"peptide fragment of {len(seq)} residues"))
                continue
            if pid in seen:
                raise ValueError(f"duplicate allergen ID {pid}")
            seen.add(pid)
            out.write(f">{pid}\n{seq}\n")
            a = allergen.get(str(r["AllergenID"]), {})
            evidence = " ".join(str(a.get("Allergenicity") or "").split())[:300]
            meta.write(
                "\t".join(
                    [
                        pid,
                        a.get("Name", ""),
                        a.get("Species", ""),
                        a.get("TaxOrder", ""),
                        a.get("Exposure", ""),
                        evidence,
                        str(len(seq)),
                    ]
                ).replace("nan", "")
                + "\n"
            )
            n += 1
    return n, skipped


def read_meta(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def parse_blast(path):
    """Return the best local alignment per query: ``{query: {...}}`` (highest bit score)."""
    best = {}
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) != len(BLAST_FIELDS.split()):
            raise ValueError(
                f"{path}:{n}: expected {len(BLAST_FIELDS.split())} fields, found {len(f)}"
            )
        q, s, pident, length, qlen, slen, bits, ev = f
        hit = {
            "subject": s,
            "identity": float(pident),
            "length": int(length),
            "qlen": int(qlen),
            "slen": int(slen),
            "bitscore": float(bits),
            "evalue": float(ev),
        }
        cur = best.get(q)
        if cur is None or (hit["bitscore"], hit["identity"]) > (cur["bitscore"], cur["identity"]):
            best[q] = hit
    return best


def allergen_name(subject):
    """``Asp_f_1.0101|11`` -> ``Asp_f_1.0101``."""
    return subject.split("|", 1)[0]


def species_code(subject):
    """``Asp_f_1.0101|11`` -> ``Asp_f`` (genus and species letters of the allergen name)."""
    parts = allergen_name(subject).split("_")
    return "_".join(parts[:2]) if len(parts) >= 3 else parts[0]


def allergen_rows(proteins, best, meta=None):
    meta = meta or {}
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        h = best.get(p.id)
        if h is None:
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "identity": "0",
                    "aligned_length": "0",
                    "coverage": "0",
                    "allergen_name": "",
                    "allergen_species": "",
                    "exposure": "",
                    "evidence": "",
                }
            )
            continue
        m = meta.get(h["subject"], {})
        cov = min(100.0, 100.0 * h["length"] / h["slen"]) if h["slen"] else 0.0
        rows.append(
            {
                "id": p.id,
                "state": "ok",
                "identity": f"{h['identity']:.2f}",
                "aligned_length": str(h["length"]),
                "coverage": f"{cov:.1f}",
                "allergen_name": allergen_name(h["subject"]),
                "allergen_species": m.get("species", ""),
                "exposure": m.get("exposure", ""),
                "evidence": m.get("evidence", ""),
            }
        )
    return rows


def check_blast_ids(best, fasta_ids):
    """A BLAST table whose queries are not FASTA IDs belongs to another proteome."""
    unknown = sorted(set(best) - set(fasta_ids))
    if unknown:
        raise ValueError(
            f"{len(unknown)} BLAST query ID(s) are not in the FASTA, for example {unknown[0]!r}"
        )


def parse_blast_hits(path):
    """All rows of a BLAST table as dicts."""
    hits = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) != len(BLAST_FIELDS.split()):
            raise ValueError(
                f"{path}:{n}: expected {len(BLAST_FIELDS.split())} fields, found {len(f)}"
            )
        hits.append(
            {
                "query": f[0],
                "subject": f[1],
                "identity": float(f[2]),
                "length": int(f[3]),
                "qlen": int(f[4]),
                "slen": int(f[5]),
                "bitscore": float(f[6]),
                "evalue": float(f[7]),
            }
        )
    return hits


def best_hit_other_species(hits, ids):
    """Per query, the best hit (bit score) to an allergen of another species (see ``species_code``)."""
    best = {}
    species = {i: species_code(i) for i in ids}
    for h in hits:
        q, s = h["query"], h["subject"]
        if q == s or q not in species or s not in species or species[q] == species[s]:
            continue
        cur = best.get(q)
        if cur is None or (h["bitscore"], h["identity"]) > (cur["bitscore"], cur["identity"]):
            best[q] = h
    return best


def recall_by_rule(ids, best_other, rules):
    """Share of allergens whose best hit in ANOTHER species meets each rule.

    ``rules`` is a list of ``(name, min_identity, min_aligned_length, min_coverage)``. These are the
    rules the engine applies, so each number belongs to one call. The denominator is ``ids`` (all
    sequences of the FASTA); a sequence with no hit counts as not recovered. This is sensitivity only.
    """
    out = []
    for name, min_id, min_len, min_cov in rules:
        k = 0
        for i in ids:
            h = best_other.get(i)
            if (
                h
                and h["identity"] >= min_id
                and h["length"] >= min_len
                and 100.0 * h["length"] / h["slen"] >= min_cov
            ):
                k += 1
        out.append({"rule": name, "recovered": k, "n": len(ids)})
    return out


def lso_report(blast_path, fasta_ids, rules=None):
    """Leave-species-out recall of the allergen set against itself (all-against-all BLAST table).

    Every sequence of ``fasta_ids`` must appear as a query (a sequence with no self-hit means that
    BLAST masked or dropped it); otherwise the function refuses.
    """
    rules = rules or [
        ("iuis_allergen_similarity", 35.0, 80, 0.0),
        ("iuis_allergen_homolog", 70.0, 0, 80.0),
    ]
    hits = parse_blast_hits(blast_path)
    queries = {h["query"] for h in hits}
    missing = sorted(set(fasta_ids) - queries)
    if missing:
        raise ValueError(
            f"{len(missing)} sequence(s) have no BLAST line, for example {missing[0]!r}"
        )
    best = best_hit_other_species(hits, fasta_ids)
    species = {species_code(i) for i in fasta_ids}
    return {
        "n_sequences": len(fasta_ids),
        "n_species": len(species),
        "recall": recall_by_rule(sorted(fasta_ids), best, rules),
    }
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_allergen.py -q`
Expected: `9 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/allergen.py tests/cellsurface_sorting_hat/modules/test_allergen.py
git commit -m "feat(sorting-hat): allergen homology module and leave-species-out recall

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 5: Lookups: antigen ranking, Cys-rich tiers, expression, TM helices

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/lookups.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_lookups.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `read_table`, `load_protein_map`, `ranking_by_gene`, `gene_of_ranking_id`, `antigen_rows` (adds `idmap_method = gene_best_transcript`), `cys_rows`, `expression_rows`, `tm_rows` (fields `n_tm`, `n_tm_mature`, `topology`), and the `*_COLUMNS` lists.
- States: `not_applicable` (taxon not in `applicable_taxa`), `not_in_reference` (no map entry or table row), `ok`. The RefSeq protein maps to a gene and the gene to its best-ranked transcript: 126 RS isoforms inherit the values of another transcript (spec 3.2 asks for an exact sequence match first; not implemented, listed as a gap).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_lookups.py` with exactly this content:

```python
import gzip

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.lookups import (
    antigen_rows,
    cys_rows,
    expression_rows,
    gene_of_ranking_id,
    load_protein_map,
    ranking_by_gene,
    read_table,
    tm_rows,
)

RANKING = (
    "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
    "CIMG_04613-t26_1-p1\t651\t7.12\t2.5\t0.0\t0.9201\t0.0\n"
    "CIMG_04613-t26_2-p1\t900\t9.9\t2.0\t0.0\t0.9201\t0.0\n"
    "CIMG_09560-t26_1-p1\t986\t10.79\t0.1\t1.0\t0.9877\t69.4\n"
)
PMAP = "protein_id\tgene_id\tproduct\tlength\nXP_1\tCIMG_04613\tp\t324\nXP_2\tCIMG_09560\tp\t100\nXP_9\tCIMG_99999\tp\t10\n"


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


def test_gene_of_ranking_id():
    assert gene_of_ranking_id("CIMG_04613-t26_1-p1") == "CIMG_04613"


def test_ranking_keeps_the_best_row_per_gene_and_counts_genes_with_several(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(RANKING)
    by_gene, several = ranking_by_gene(path)
    assert by_gene["CIMG_04613"]["rank"] == "651" and several == 1


def test_antigen_states_and_fields(tmp_path):
    (tmp_path / "r.tsv").write_text(RANKING)
    (tmp_path / "m.tsv").write_text(PMAP)
    by_gene, _ = ranking_by_gene(tmp_path / "r.tsv")
    pmap = load_protein_map(tmp_path / "m.tsv")
    taxa = {
        "XP_1": 246410,
        "XP_2": 246410,
        "XP_9": 246410,
        "XP_X": 246410,
        "AF1": 746128,
        "BAD": 246410,
    }
    rows = antigen_rows(
        [prot(p, "na_invalid" if p == "BAD" else "ok") for p in taxa], taxa, pmap, by_gene, {246410}
    )
    got = {r["id"]: r for r in rows}
    assert (
        got["XP_1"]["state"] == "ok"
        and got["XP_1"]["percentile"] == "7.12"
        and got["XP_1"]["max_crossreact"] == "0.0"
    )
    assert got["XP_2"]["percentile"] == "10.79"
    assert got["XP_9"]["state"] == "not_in_reference"  # gene not in the ranking
    assert got["XP_X"]["state"] == "not_in_reference"  # protein not in the map
    assert got["AF1"]["state"] == "not_applicable"  # another taxon
    assert got["BAD"]["state"] == "na_invalid"


def test_read_table_counts_repeated_keys(tmp_path):
    path = tmp_path / "t.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("k\tv\na\t1\na\t2\nb\t3\n")
    rows, repeated = read_table(path, "k")
    assert rows["a"]["v"] == "1" and repeated == 1
    with pytest.raises(ValueError, match="missing column"):
        read_table(path, "nope")


def test_cys_expression_and_tm_rows():
    taxa = {"A": 246410, "B": 246410, "C": 746128}
    cys = cys_rows(
        [prot(p) for p in "ABC"],
        taxa,
        {"A": {"tier": "cys_rich_sp_unassigned", "cys_frac": "0.1"}},
        {246410},
    )
    assert [(r["id"], r["state"]) for r in cys] == [
        ("A", "ok"),
        ("B", "not_in_reference"),
        ("C", "not_applicable"),
    ]
    assert cys[0]["tier"] == "cys_rich_sp_unassigned"
    expr = expression_rows(
        [prot("A"), prot("C")],
        taxa,
        {"A": "CIMG_1"},
        {"CIMG_1": {"log2fc_48h": "9.09", "padj_48h": "1e-5", "log2fc_8d": "8.0"}},
        {246410},
    )
    assert expr[0]["log2fc"] == "9.09" and expr[1]["state"] == "not_applicable"
    tm = tm_rows([prot("A"), prot("B")], {"A": {"pred_hel": "7", "topology": "o10-32i"}})
    assert tm[0]["n_tm"] == "7" and tm[1]["state"] == "error"


def test_tm_rows_count_helices_after_the_signal_peptide_window():
    from cellsurface_sorting_hat.modules.lookups import tm_rows

    tmhmm = {
        "SP_ONLY": {"pred_hel": "1", "topology": "o10-32i"},
        "RECEPTOR": {
            "pred_hel": "7",
            "topology": "i40-62o70-92i100-122o130-152i160-182o190-212i220-242o",
        },
        "MIXED": {"pred_hel": "2", "topology": "o7-25i40-62o"},
        "NONE": {"pred_hel": "0", "topology": "o"},
    }
    got = {r["id"]: (r["n_tm"], r["n_tm_mature"]) for r in tm_rows([prot(k) for k in tmhmm], tmhmm)}
    assert got == {"SP_ONLY": ("1", 0), "RECEPTOR": ("7", 7), "MIXED": ("2", 1), "NONE": ("0", 0)}
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_lookups.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.lookups`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/lookups.py` with exactly this content:

```python
"""Lookups into precomputed tables: antigen ranking, Cys-rich tiers, spherule expression, TM helices.

These modules need an ID map. For *C. immitis* RS the RefSeq protein ID (``XP_...``) maps to a
gene ID (``CIMG_...``) through ``protein_map.tsv``. The antigen ranking and the spherule table are
keyed on that gene ID. A protein with no map entry or no table row is ``not_in_reference``. A protein
from a taxon where the table does not apply is ``not_applicable``.
"""

import csv
import gzip
import re
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

ANTIGEN_COLUMNS = [
    "percentile",
    "antigenicity",
    "specificity",
    "prevalence",
    "max_crossreact",
    "rank",
    "idmap_method",
]
CYS_COLUMNS = ["tier", "cys_frac"]
EXPRESSION_COLUMNS = ["log2fc", "padj", "log2fc_8d"]
TM_COLUMNS = ["n_tm", "n_tm_mature", "topology"]
SIGNAL_PEPTIDE_WINDOW = (
    35  # a helix that starts at or before this residue may be the signal peptide
)


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def read_table(path, key, delimiter="\t"):
    """Read a TSV into ``{row[key]: row}``. A repeated key keeps the first row and is counted."""
    rows, repeated = {}, 0
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if key not in (reader.fieldnames or []):
            raise ValueError(f"{path}: missing column {key!r}")
        for r in reader:
            if r[key] in rows:
                repeated += 1
            else:
                rows[r[key]] = r
    return rows, repeated


def load_protein_map(path):
    rows, _ = read_table(path, "protein_id")
    return {pid: r["gene_id"] for pid, r in rows.items()}


def gene_of_ranking_id(protein):
    """``CIMG_04613-t26_1-p1`` -> ``CIMG_04613``."""
    return protein.split("-", 1)[0]


def ranking_by_gene(ranking_tsv):
    """Best (lowest rank) ranking row per gene. Returns ``(rows, n_genes_with_several_rows)``."""
    rows, several = {}, set()
    with _open_text(ranking_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = gene_of_ranking_id(r["protein"])
            cur = rows.get(gene)
            if cur is not None:
                several.add(gene)
            if cur is None or int(r["rank"]) < int(cur["rank"]):
                rows[gene] = r
    return rows, len(several)


def _applicable(taxa, pid, applicable_taxa):
    return taxa[pid] in applicable_taxa


def antigen_rows(proteins, taxa, protein_map, by_gene, applicable_taxa):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        else:
            gene = protein_map.get(p.id) or (
                gene_of_ranking_id(p.id) if p.id.startswith("CIMG_") else None
            )
            r = by_gene.get(gene)
            if r is None:
                rows.append({"id": p.id, "state": "not_in_reference"})
            else:
                rows.append(
                    {
                        "id": p.id,
                        "state": "ok",
                        "percentile": r["percentile"],
                        "antigenicity": r["antigenicity"],
                        "specificity": r["specificity"],
                        "prevalence": r["prevalence"],
                        "max_crossreact": r["max_fungal_crossreact_pid"],
                        "rank": r["rank"],
                        "idmap_method": "gene_best_transcript",
                    }
                )
    return rows


def cys_rows(proteins, taxa, cys_table, applicable_taxa):
    """``cys_table``: ``{protein_id: row}`` from ``candidates.tsv.gz`` of ``analysis/cys_candidates``."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        elif p.id not in cys_table:
            rows.append({"id": p.id, "state": "not_in_reference"})
        else:
            r = cys_table[p.id]
            rows.append({"id": p.id, "state": "ok", "tier": r["tier"], "cys_frac": r["cys_frac"]})
    return rows


def expression_rows(proteins, taxa, protein_map, spherule, applicable_taxa):
    """``spherule``: ``{gene_id: row}`` from ``spherule_surface_table.tsv.gz``."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        else:
            r = spherule.get(protein_map.get(p.id, ""))
            if r is None:
                rows.append({"id": p.id, "state": "not_in_reference"})
            else:
                rows.append(
                    {
                        "id": p.id,
                        "state": "ok",
                        "log2fc": r["log2fc_48h"],
                        "padj": r["padj_48h"],
                        "log2fc_8d": r["log2fc_8d"],
                    }
                )
    return rows


def tm_rows(proteins, tmhmm):
    """``tmhmm``: ``{protein_id: row}`` from the TMHMM table (columns ``pred_hel``, ``topology``)."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in tmhmm:
            rows.append({"id": p.id, "state": "error"})
        else:
            r = tmhmm[p.id]
            starts = [int(a) for a, _ in re.findall(r"(\d+)-(\d+)", r["topology"])]
            mature = sum(1 for a in starts if a > SIGNAL_PEPTIDE_WINDOW)
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "n_tm": r["pred_hel"],
                    "n_tm_mature": mature,
                    "topology": r["topology"],
                }
            )
    return rows
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_lookups.py -q`
Expected: `6 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/lookups.py tests/cellsurface_sorting_hat/modules/test_lookups.py
git commit -m "feat(sorting-hat): antigen, Cys-rich, expression and TM lookups

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 6: Intervals and measures

**Files:**
- Create: `src/cellsurface_sorting_hat/calibration/__init__.py`, `src/cellsurface_sorting_hat/calibration/intervals.py`, `src/cellsurface_sorting_hat/calibration/measure.py`
- Test: `tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py`

**Interfaces:**
- Consumes: `write_atomic` (Plan 1), `load_status_source`, status constants (Plan 1, patched in Task 0).
- Produces: `wilson(k, n)` (refuses k > n), `cluster_bootstrap(y, call, clusters, n_boot=2000, seed=1)` (widest of bootstrap and Wilson on the number of clusters; returns the cluster counts), `build_measure(...)` (records `n_clusters_pos`, `n_clusters_neg`), `status_from_measure(measure)`, `make_entry(taxa, measure, source='', cap=None)`, `write_status_source(workdir, module, entries)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py` with exactly this content:

```python
import json

import numpy as np
import pytest

from cellsurface_sorting_hat.calibration.intervals import cluster_bootstrap, wilson
from cellsurface_sorting_hat.calibration.measure import (
    build_measure,
    make_entry,
    status_from_measure,
    write_status_source,
)
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source


def test_wilson_matches_known_values():
    v, lo, hi = wilson(8, 10)
    assert (
        v == 0.8 and lo == pytest.approx(0.4902, abs=1e-3) and hi == pytest.approx(0.9433, abs=1e-3)
    )
    assert wilson(0, 10)[1] == 0.0 and wilson(10, 10)[2] == pytest.approx(1.0)
    assert wilson(0, 0) is None


def test_bootstrap_point_estimates_and_ordering():
    y = np.array([1] * 30 + [0] * 70)
    call = np.array([True] * 24 + [False] * 6 + [True] * 7 + [False] * 63)
    clusters = [f"c{i}" for i in range(100)]
    r = cluster_bootstrap(y, call, clusters, n_boot=500, seed=3)
    assert r["sensitivity"]["value"] == pytest.approx(0.8) and r["specificity"][
        "value"
    ] == pytest.approx(0.9)
    for k in (r["sensitivity"], r["specificity"]):
        assert 0 <= k["lo"] <= k["value"] <= k["hi"] <= 1
    assert r == cluster_bootstrap(y, call, clusters, n_boot=500, seed=3)  # same seed, same result


def test_clusters_widen_the_interval_compared_with_single_proteins():
    y = np.array([1] * 40 + [0] * 40)
    call = np.array(([True] * 20 + [False] * 20) * 2)
    single = cluster_bootstrap(y, call, [str(i) for i in range(80)], n_boot=800, seed=1)[
        "sensitivity"
    ]
    # 8 clusters of 5 positives that are all called or all missed: far less independent information
    pos_clusters = [f"p{i // 5}" for i in range(40)]
    call2 = np.array([(i // 5) % 2 == 0 for i in range(40)] + [False] * 40)
    grouped = cluster_bootstrap(
        y, call2, pos_clusters + [f"n{i}" for i in range(40)], n_boot=800, seed=1
    )["sensitivity"]
    assert (grouped["hi"] - grouped["lo"]) > (single["hi"] - single["lo"])


def test_bootstrap_without_negatives_has_no_specificity():
    r = cluster_bootstrap([1, 1, 1], [True, False, True], ["a", "b", "c"], n_boot=100)
    assert r["specificity"] is None and r["sensitivity"] is not None


@pytest.mark.parametrize(
    "args", [([1, 0], [True], ["a", "b"]), ([2, 0], [True, False], ["a", "b"])]
)
def test_bootstrap_refuses_bad_input(args):
    with pytest.raises(ValueError):
        cluster_bootstrap(*args)


def test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate():
    def m(n_pos, sens, n_neg=None, spec=None):
        out = {"n_pos": n_pos, "sensitivity": sens}
        if spec:
            out.update(n_neg=n_neg, specificity=spec)
        return out

    tight = {"value": 0.6, "lo": 0.55, "hi": 0.65}
    wide = {"value": 0.8, "lo": 0.5, "hi": 0.95}
    sp_ok = {"value": 0.96, "lo": 0.95, "hi": 0.97}
    sp_wide = {"value": 0.9, "lo": 0.6, "hi": 0.99}
    assert status_from_measure(m(232, tight, 4244, sp_ok)) == "estimated"
    assert status_from_measure(m(25, wide, 100, sp_ok)) == "smoke"  # sensitivity half-width 0.225
    assert status_from_measure(m(16, tight, 100, sp_ok)) == "smoke"  # fewer than 20 positives
    assert status_from_measure(m(232, tight, 15, sp_ok)) == "smoke"  # fewer than 20 negatives
    assert status_from_measure(m(232, tight, 4244, sp_wide)) == "smoke"  # specificity too wide
    assert status_from_measure(m(232, tight)) == "smoke"  # never tested on negatives
    # negatives were counted but no specificity was measured: still not an estimate
    assert status_from_measure({"n_pos": 232, "n_neg": 500, "sensitivity": tight}) == "smoke"
    assert status_from_measure({"n_pos": 0}) == "unvalidated"


def test_a_leakage_cap_limits_an_estimate_to_smoke():
    measure = {
        "n_pos": 232,
        "n_neg": 4244,
        "sensitivity": {"value": 0.6, "lo": 0.55, "hi": 0.65},
        "specificity": {"value": 0.96, "lo": 0.95, "hi": 0.97},
    }
    assert make_entry([1], measure)["status"] == "estimated"
    assert make_entry([1], measure, cap="smoke")["status"] == "smoke"


def test_build_measure_counts_and_notes():
    m = build_measure(
        "set",
        "truth.tsv",
        [1, 1, 0, 0],
        [True, False, False, False],
        ["a", "b", "c", "d"],
        notes="n",
        n_boot=50,
    )
    assert (m["n_pos"], m["n_neg"], m["calibration_set"], m["notes"]) == (2, 2, "set", "n")


def test_write_status_source_uses_the_module_identity_and_validates(tmp_path):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    m = build_measure("set", "t", [1, 0], [True, False], ["a", "b"], n_boot=20)
    path = write_status_source(tmp_path, "pfam_adhesion", [make_entry([40], m)])
    rec = load_status_source(path)
    run = json.loads((tmp_path / "modules" / "pfam_adhesion.json").read_text())
    assert rec.identity.params_hash == run["params_hash"] and rec.entries[0].taxa == (40,)
    with pytest.raises(ValueError, match="appears in two entries"):
        write_status_source(
            tmp_path, "pfam_adhesion", [make_entry([40], m), make_entry([40, 41], m)]
        )


def test_all_positives_called_gives_a_width_not_a_zero_width_interval():
    y = [1] * 20 + [0] * 20
    r = cluster_bootstrap(y, [True] * 20 + [False] * 20, [f"c{i}" for i in range(40)], n_boot=300)
    assert (
        r["sensitivity"]["value"] == 1.0 and r["sensitivity"]["lo"] < 0.9
    )  # Wilson on 20 clusters
    assert r["specificity"]["value"] == 1.0 and r["specificity"]["lo"] < 0.9
    assert (r["n_clusters_pos"], r["n_clusters_neg"]) == (20, 20)


def test_one_cluster_per_class_cannot_give_an_estimate():
    y = [1] * 30 + [0] * 30
    call = [True] * 24 + [False] * 6 + [False] * 28 + [True] * 2
    clusters = ["pos"] * 30 + ["neg"] * 30
    m = build_measure("s", "t", y, call, clusters, n_boot=200)
    assert (m["n_clusters_pos"], m["n_clusters_neg"]) == (1, 1)
    assert status_from_measure(m) == "smoke"  # many proteins, almost no independent clusters
    assert (m["sensitivity"]["hi"] - m["sensitivity"]["lo"]) / 2 > 0.10


def test_wilson_refuses_k_above_n():
    with pytest.raises(ValueError):
        wilson(5, 3)


def test_a_tight_interval_with_few_recorded_clusters_is_not_an_estimate():
    tight = {"value": 0.6, "lo": 0.55, "hi": 0.65}
    base = {
        "n_pos": 300,
        "n_neg": 3000,
        "sensitivity": tight,
        "specificity": {"value": 0.96, "lo": 0.95, "hi": 0.97},
    }
    assert status_from_measure({**base, "n_clusters_pos": 40, "n_clusters_neg": 300}) == "estimated"
    assert status_from_measure({**base, "n_clusters_pos": 5, "n_clusters_neg": 300}) == "smoke"
    assert status_from_measure({**base, "n_clusters_pos": 40, "n_clusters_neg": 3}) == "smoke"
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.calibration`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/calibration/__init__.py` with exactly this content:

```python
"""Calibration: sensitivity and specificity of each module, and the status sources that record them."""
```

Create `src/cellsurface_sorting_hat/calibration/intervals.py` with exactly this content:

```python
"""Confidence intervals for sensitivity and specificity.

``wilson`` is for simple counts. ``cluster_bootstrap`` resamples homology clusters (not single
proteins), so related proteins do not make the interval too narrow. Models are not refitted; the
interval describes sampling of the calibration set only.
"""

import math

import numpy as np


def wilson(k, n, z=1.959964):
    """Wilson score interval for k successes in n trials: ``(value, lo, hi)``; None when n is 0."""
    if n <= 0:
        return None
    if not 0 <= k <= n:
        raise ValueError(f"k must be between 0 and n (k={k}, n={n})")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)


def _rates(y, call):
    pos, neg = y == 1, y == 0
    sens = (call & pos).sum() / pos.sum() if pos.sum() else np.nan
    spec = (~call & neg).sum() / neg.sum() if neg.sum() else np.nan
    return sens, spec


def cluster_bootstrap(y, call, clusters, n_boot=2000, seed=1):
    """Sensitivity and specificity with a 95% interval that respects clusters.

    ``y`` is 1 for a known positive and 0 for a known negative; ``call`` is True when the module
    called the protein. The interval is the widest of two: the percentile interval of a cluster
    bootstrap, and the Wilson interval on the number of independent clusters of that class (the
    effective sample size). The second guards against two failures of the bootstrap: zero width when
    every resample gives the same value (all positives called), and too narrow an interval when there
    are few clusters. Returns ``{"sensitivity", "specificity", "n_clusters_pos", "n_clusters_neg"}``;
    a rate is None when its class is empty. Each rate is ``{"value", "lo", "hi"}``, ``lo <= value <= hi``.
    """
    y = np.asarray(y, dtype=int)
    call = np.asarray(call, dtype=bool)
    clusters = np.asarray(clusters, dtype=str)
    if not (len(y) == len(call) == len(clusters)):
        raise ValueError("y, call and clusters must have the same length")
    if not set(np.unique(y)) <= {0, 1}:
        raise ValueError("y must hold 0 and 1 only")
    uniq, inverse = np.unique(clusters, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(len(uniq))]
    point = _rates(y, call)
    rng = np.random.default_rng(seed)
    boot = np.full((n_boot, 2), np.nan)
    for b in range(n_boot):
        draw = rng.integers(0, len(uniq), size=len(uniq))
        idx = np.concatenate([members[i] for i in draw])
        boot[b] = _rates(y[idx], call[idx])
    n_clusters = {1: len(set(clusters[y == 1])), 0: len(set(clusters[y == 0]))}
    out = {"n_clusters_pos": n_clusters[1], "n_clusters_neg": n_clusters[0]}
    for j, (name, label) in enumerate((("sensitivity", 1), ("specificity", 0))):
        value = point[j]
        if np.isnan(value):
            out[name] = None
            continue
        col = boot[:, j][~np.isnan(boot[:, j])]
        lo, hi = np.percentile(col, [2.5, 97.5]) if len(col) else (value, value)
        nc = n_clusters[label]
        _, wlo, whi = wilson(int(round(value * nc)), nc)
        out[name] = {
            "value": float(value),
            "lo": float(min(lo, wlo, value)),
            "hi": float(max(hi, whi, value)),
        }
    return out
```

Create `src/cellsurface_sorting_hat/calibration/measure.py` with exactly this content:

```python
"""Build ``measure`` objects and write status sources (the files that record calibration)."""

import json
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.calibration.intervals import cluster_bootstrap
from cellsurface_sorting_hat.status import ESTIMATED, SMOKE, UNVALIDATED, load_status_source

MIN_POSITIVES = 20
MIN_NEGATIVES = 20
MIN_CLUSTERS = 20
MAX_HALF_WIDTH = 0.10


def build_measure(calibration_set, truth_source, y, call, clusters, notes="", n_boot=2000, seed=1):
    rates = cluster_bootstrap(y, call, clusters, n_boot=n_boot, seed=seed)
    y = list(y)
    measure = {
        "calibration_set": calibration_set,
        "truth_source": truth_source,
        "n_pos": sum(1 for v in y if v == 1),
        "n_neg": sum(1 for v in y if v == 0),
        "n_clusters_pos": rates["n_clusters_pos"],
        "n_clusters_neg": rates["n_clusters_neg"],
    }
    for key in ("sensitivity", "specificity"):
        if rates[key] is not None:
            measure[key] = rates[key]
    if notes:
        measure["notes"] = notes
    return measure


def status_from_measure(
    measure,
    min_positives=MIN_POSITIVES,
    min_negatives=MIN_NEGATIVES,
    max_half_width=MAX_HALF_WIDTH,
    min_clusters=MIN_CLUSTERS,
):
    """``estimated`` needs at least 20 positives and 20 negatives (and, when recorded, at least 20
    independent clusters of each) and a 95% interval half-width of at most 0.10 for both
    sensitivity and specificity. Phase C uses the same floor for recall. A measure with any
    positive but without a specificity is at most ``smoke``: a rule that was never tested on
    negatives is not an estimate. No sensitivity or no positive gives ``unvalidated``."""
    sens, spec = measure.get("sensitivity"), measure.get("specificity")
    if not sens or measure.get("n_pos", 0) < 1:
        return UNVALIDATED
    narrow = (sens["hi"] - sens["lo"]) / 2 <= max_half_width
    if spec:
        narrow = narrow and (spec["hi"] - spec["lo"]) / 2 <= max_half_width
    enough = measure["n_pos"] >= min_positives and measure.get("n_neg", 0) >= min_negatives
    # independent clusters, when the measure records them (a Phase C measure does not)
    for key in ("n_clusters_pos", "n_clusters_neg"):
        if key in measure and measure[key] < min_clusters:
            enough = False
    if spec and enough and narrow:
        return ESTIMATED
    return SMOKE


def make_entry(taxa, measure, source="", cap=None):
    """A status entry. ``cap="smoke"`` limits the status, for example when the truth set overlaps
    the data that set the rule."""
    status = status_from_measure(measure)
    if cap == SMOKE and status == ESTIMATED:
        status = SMOKE
    return {"taxa": [int(t) for t in taxa], "status": status, "source": source, "measure": measure}


def write_status_source(workdir, module, entries):
    """Write ``<workdir>/status/<module>.json`` using the identity in ``modules/<module>.json``.

    The file is validated by reading it back with ``load_status_source``. Entries must not share a
    tested taxon (the engine would pick the first of equal depth).
    """
    record = json.loads((Path(workdir) / "modules" / f"{module}.json").read_text())
    seen = {}
    for e in entries:
        for t in e["taxa"]:
            if t in seen:
                raise ValueError(
                    f"taxon {t} appears in two entries ({seen[t]} and {e['measure']['calibration_set']})"
                )
            seen[t] = e["measure"]["calibration_set"]
    data = {k: record[k] for k in ("module", "version", "params_hash", "artefact_hash")}
    data["entries"] = entries
    path = Path(workdir) / "status" / f"{module}.json"
    write_atomic(path, (json.dumps(data, indent=2, sort_keys=True) + "\n").encode())
    load_status_source(path)  # raises if the file is not valid
    return path
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py -q`
Expected: `14 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/calibration tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py
git commit -m "feat(sorting-hat): cluster-aware intervals and status measures

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 7: Phase C converter and the strain stability comparison

**Files:**
- Create: `src/cellsurface_sorting_hat/calibration/phasec.py`, `src/cellsurface_sorting_hat/calibration/stability.py`, `data/sorting_hat/phasec_set_species.tsv`
- Test: `tests/cellsurface_sorting_hat/calibration/test_phasec.py`, `tests/cellsurface_sorting_hat/calibration/test_stability.py`

**Interfaces:**
- Consumes: `make_entry` (Task 6), `wilson` (Task 6); `TaxonError` (Plan 1).
- Produces: `read_names`, `species_taxid`, `entries_from_phasec(metrics_json, set_taxa, candidate='R0', variant='V-go', truth='direct')` (one entry per species test set; status is the weaker of the Phase C label and `status_from_measure`; per-stratum specificity in `strata`; zero-width intervals at 0 or 1 are widened with Wilson), `SETS`; `compare_runs(out_a, out_b)`.
- Specificity is one minus the false-positive rate with the interval ends swapped (valid: the percentile interval is equivariant under `1 - x`).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/calibration/test_phasec.py` with exactly this content:

```python
import json

import pytest

from cellsurface_sorting_hat.calibration.measure import write_status_source
from cellsurface_sorting_hat.calibration.phasec import (
    SETS,
    entries_from_phasec,
    read_names,
    species_taxid,
)
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source
from cellsurface_sorting_hat.taxonomy import TaxonError


def cell(r, f):
    return {
        "recall": {"value": r[0], "lo": r[1], "hi": r[2]},
        "fpr": {"value": f[0], "lo": f[1], "hi": f[2]} if f else {"value": None},
    }


def make_set(label, pos, neg, c, strata=None):
    metrics = {"all": {"V-go": {"R0": c}}}
    for name, f in (strata or {}).items():
        metrics[name] = {"V-go": {"R0": {"recall": {"value": None}, "fpr": f}}}
    return {
        "label": label,
        "n_direct_positives": pos,
        "truth": {"direct": {"n": {"all": {"pos": pos, "neg": neg}}, "metrics": metrics}},
    }


def _metrics(tmp_path, scer_label="estimate"):
    data = {
        "test_sets": {
            "S1:Scer_SGD": make_set(
                scer_label,
                79,
                2000,
                cell((0.848, 0.77, 0.92), (0.03, 0.02, 0.04)),
                {
                    "N-sec": {"value": 0.084, "lo": 0.06, "hi": 0.11},
                    "N-int": {"value": 0.002, "lo": 0.001, "hi": 0.004},
                },
            ),
            "S1:Calb_CGD": make_set(
                "smoke test", 153, 2244, cell((0.477, 0.36, 0.59), (0.04, 0.03, 0.05))
            ),
            "S3-Eurotiomycetes:Afum_ASPFU": make_set(
                "smoke test", 19, 44, cell((0.947, 0.78, 1.0), (0.02, 0.0, 0.08))
            ),
        }
    }
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps(data))
    return path


SET_TAXA = {"S1:Scer_SGD": [4932], "S1:Calb_CGD": [5476], "S3-Eurotiomycetes:Afum_ASPFU": [746128]}


def test_each_species_gets_its_own_entry_with_its_own_numbers(tmp_path):
    entries = entries_from_phasec(_metrics(tmp_path), SET_TAXA)
    by = {e["measure"]["calibration_set"]: e for e in entries}
    assert [e["taxa"] for e in entries] == [[4932], [5476], [746128]]
    assert by["S1:Scer_SGD"]["measure"]["sensitivity"] == {"value": 0.848, "lo": 0.77, "hi": 0.92}
    assert by["S1:Calb_CGD"]["measure"]["sensitivity"]["value"] == 0.477  # not a pooled value
    sp = by["S1:Scer_SGD"]["measure"]["specificity"]
    assert (
        sp["value"] == pytest.approx(0.97)
        and sp["lo"] == pytest.approx(0.96)
        and sp["hi"] == pytest.approx(0.98)
    )
    assert (by["S1:Scer_SGD"]["measure"]["n_pos"], by["S1:Scer_SGD"]["measure"]["n_neg"]) == (
        79,
        2000,
    )


def test_per_stratum_specificity_is_kept(tmp_path):
    scer = entries_from_phasec(_metrics(tmp_path), SET_TAXA)[0]["measure"]
    assert scer["strata"]["N-sec"]["value"] == pytest.approx(0.916)
    assert scer["strata"]["N-int"]["lo"] == pytest.approx(0.996)
    assert "strata" not in entries_from_phasec(_metrics(tmp_path), SET_TAXA)[1]["measure"]


def test_the_status_is_the_weaker_of_the_phase_c_label_and_the_rule(tmp_path):
    entries = {
        e["measure"]["calibration_set"]: e
        for e in entries_from_phasec(_metrics(tmp_path), SET_TAXA)
    }
    assert entries["S1:Scer_SGD"]["status"] == "estimated"  # label estimate and the rule agrees
    assert (
        entries["S3-Eurotiomycetes:Afum_ASPFU"]["status"] == "smoke"
    )  # label smoke test, 19 positives
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 400, 4000, cell((0.6, 0.3, 0.9), (0.05, 0.04, 0.06))
    )
    path = tmp_path / "wide.json"
    path.write_text(json.dumps(data))
    wide = {e["measure"]["calibration_set"]: e for e in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]
    assert wide["status"] == "smoke"  # labelled estimate, but the sensitivity interval is wide


def test_a_set_with_no_specificity_cell_is_at_most_smoke_even_if_labelled_estimate(tmp_path):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 400, 4000, cell((0.5, 0.48, 0.52), None)
    )
    path = tmp_path / "m2.json"
    path.write_text(json.dumps(data))
    e = {x["measure"]["calibration_set"]: x for x in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]
    assert "specificity" not in e["measure"] and e["status"] == "smoke"


def test_a_tight_estimate_with_both_rates_is_estimated(tmp_path):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 400, 4000, cell((0.5, 0.48, 0.52), (0.05, 0.04, 0.06))
    )
    path = tmp_path / "m3.json"
    path.write_text(json.dumps(data))
    e = {x["measure"]["calibration_set"]: x for x in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]
    assert e["status"] == "estimated"


def test_entries_make_a_valid_status_source(tmp_path):
    write_module(tmp_path, ModuleSpec("step1_rule@R0", "1"), [], [{"id": "A", "state": "ok"}])
    path = write_status_source(
        tmp_path, "step1_rule@R0", entries_from_phasec(_metrics(tmp_path), SET_TAXA)
    )
    assert len(load_status_source(path).entries) == 3


def test_default_sets_are_single_species_and_share_no_source():
    assert all(len(group) == 1 for group in SETS.values())
    sources = [s for group in SETS.values() for s in group]
    assert len(sources) == len(set(sources)) == 6


def test_species_taxid_reads_names_dmp_and_refuses_missing_or_ambiguous_names(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(
        "4932\t|\tSaccharomyces cerevisiae\t|\t\t|\tscientific name\t|\n"
        "4932\t|\tbaker's yeast\t|\t\t|\tcommon name\t|\n"
        "111\t|\tTwin\t|\t\t|\tscientific name\t|\n"
        "222\t|\tTwin\t|\t\t|\tscientific name\t|\n"
    )
    table = read_names(names)
    assert species_taxid(table, "Saccharomyces cerevisiae") == 4932
    for bad in ("Twin", "Nothing here"):
        with pytest.raises(TaxonError):
            species_taxid(table, bad)


def test_a_zero_width_interval_at_a_rate_of_one_is_widened_with_wilson(tmp_path):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 9, 28, cell((1.0, 1.0, 1.0), (0.0, 0.0, 0.0))
    )
    path = tmp_path / "edge.json"
    path.write_text(json.dumps(data))
    m = {e["measure"]["calibration_set"]: e for e in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]["measure"]
    assert (
        m["sensitivity"]["lo"] == pytest.approx(0.701, abs=1e-3) and m["sensitivity"]["hi"] == 1.0
    )  # 9 of 9
    assert m["specificity"]["lo"] == pytest.approx(0.879, abs=1e-3)  # 28 of 28
```

Create `tests/cellsurface_sorting_hat/calibration/test_stability.py` with exactly this content:

```python
import csv
import gzip

from cellsurface_sorting_hat.calibration.stability import compare_runs


def write_run(path, proteins, calls):
    path.mkdir()
    with gzip.open(path / "proteins.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
        )
        for pid, sha in proteins:
            w.writerow([pid, sha, 1, "ok", "", 0, 0])
    with gzip.open(path / "calls.long.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for pid, call, variant, value in calls:
            w.writerow([pid, call, variant, value, "unvalidated", "", ""])


def test_only_identical_sequences_are_compared_by_sha256_not_by_id(tmp_path):
    write_run(
        tmp_path / "a",
        [("a1", "s1"), ("a2", "s2"), ("a3", "s3")],
        [("a1", "x", "", "called"), ("a2", "x", "", "called"), ("a3", "x", "", "called")],
    )
    write_run(
        tmp_path / "b",
        [("b9", "s1"), ("b8", "s2"), ("b7", "s4")],
        [("b9", "x", "", "called"), ("b8", "x", "", "not_called"), ("b7", "x", "", "called")],
    )
    r = compare_runs(tmp_path / "a", tmp_path / "b")
    assert (r["n_a"], r["n_b"], r["n_identical_sequences"]) == (3, 3, 2)
    assert dict(r["per_call"][("x", "")]) == {"agree": 1, "differ": 1}
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_phasec.py tests/cellsurface_sorting_hat/calibration/test_stability.py -q`
Expected: collection error, `ModuleNotFoundError: ... calibration.phasec`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/calibration/phasec.py` with exactly this content:

```python
"""Turn the Phase C ``metrics.json`` into status entries for the parameter-free rule R0 (one per species).

R0 ("SignalP calls a signal peptide") has no fitted parameter, so its sensitivity and false-positive
rate on a test set are measurements of the rule itself. Only the pooled test sets in ``SETS`` are
used, so no tested taxon appears in two entries. Specificity is ``1 - FPR`` over all negative
classes of the set, with the interval ends swapped.
"""

import json
from pathlib import Path

from cellsurface_sorting_hat.calibration.intervals import wilson
from cellsurface_sorting_hat.calibration.measure import make_entry
from cellsurface_sorting_hat.taxonomy import TaxonError

# Phase C species test set -> the source ID (row of species.tsv). One entry per species: a pooled
# estimate describes neither species (S. cerevisiae and C. albicans, A. fumigatus and A. nidulans
# differ strongly).
SETS = {
    "S1:Scer_SGD": ("Scer_SGD",),
    "S1:Calb_CGD": ("Calb_CGD",),
    "S3-Eurotiomycetes:Afum_ASPFU": ("Afum_ASPFU",),
    "S3-Eurotiomycetes:Anid_EMENI": ("Anid_EMENI",),
    "S3-Basidiomycota:Cneo_H99_GOA": ("Cneo_H99_GOA",),
    "S3-Basidiomycota:Umay_MYCMD": ("Umay_MYCMD",),
}
STRATA = ("N-int", "N-sec", "PM-TM")


def read_names(names_dmp):
    """``names.dmp`` -> ``{scientific name: [taxon IDs]}``."""
    out = {}
    for line in Path(names_dmp).read_text(encoding="utf-8-sig").splitlines():
        f = [x.strip() for x in line.split("|")]
        if len(f) >= 4 and f[3] == "scientific name":
            out.setdefault(f[1], []).append(int(f[0]))
    return out


def species_taxid(names, scientific_name):
    """The one taxon ID for a scientific name; refuses a missing or ambiguous name."""
    ids = names.get(scientific_name, [])
    if len(ids) != 1:
        raise TaxonError(f"{scientific_name!r} has {len(ids)} taxon IDs in names.dmp; expected 1")
    return ids[0]


def rate_from_phasec(cell):
    if not cell or cell.get("value") is None:
        return None
    return {"value": cell["value"], "lo": cell["lo"], "hi": cell["hi"]}


def entries_from_phasec(metrics_json, set_taxa, candidate="R0", variant="V-go", truth="direct"):
    """Return a list of entries (see ``make_entry``), one per test set in ``set_taxa``.

    ``set_taxa`` maps a Phase C test set key to its species-level taxon IDs. The status is the
    weaker of the Phase C label (``estimate`` or ``smoke test``) and ``status_from_measure``. A
    Phase C measure records no cluster counts, so the cluster floor is not applied; Phase C's own
    bootstrap resampled clusters. The rates of the negative classes (N-int, N-sec, PM-TM) are kept
    in ``strata``, because the pooled specificity depends on the mix of negatives.
    """
    m = json.loads(Path(metrics_json).read_text())
    entries = []
    for key, taxa in set_taxa.items():
        ts = m["test_sets"][key]
        metrics = ts["truth"][truth]["metrics"]
        cell = metrics["all"][variant][candidate]
        sens = rate_from_phasec(cell["recall"])
        fpr = rate_from_phasec(cell["fpr"])
        measure = {
            "calibration_set": key,
            "truth_source": f"{Path(metrics_json).name}:test_sets/{key}/truth/{truth}",
            "n_pos": int(ts["truth"][truth]["n"]["all"]["pos"]),
            "n_neg": int(ts["truth"][truth]["n"]["all"]["neg"]),
            "notes": f"Phase C {ts['label']}; rule {candidate}, variant {variant}; GO direct evidence",
        }
        if sens:
            measure["sensitivity"] = widen_at_boundary(sens, measure["n_pos"])
        if fpr:
            measure["specificity"] = widen_at_boundary(_spec_from_fpr(fpr), measure["n_neg"])
        strata = {}
        for name in STRATA:
            f = rate_from_phasec(
                metrics.get(name, {}).get(variant, {}).get(candidate, {}).get("fpr")
            )
            if f:
                strata[name] = _spec_from_fpr(f)
        if strata:
            measure["strata"] = strata
        entry = make_entry(taxa, measure, source=measure["truth_source"])
        if ts["label"] != "estimate" and entry["status"] == "estimated":
            entry["status"] = "smoke"
        entries.append(entry)
    return entries


def _spec_from_fpr(fpr):
    return {"value": 1 - fpr["value"], "lo": 1 - fpr["hi"], "hi": 1 - fpr["lo"]}


def widen_at_boundary(rate, n):
    """A percentile bootstrap gives a zero-width interval at a rate of 0 or 1 (every resample
    agrees). Combine it with the Wilson interval on ``n`` counted proteins, which is wider."""
    if not rate or n <= 0 or rate["value"] not in (0.0, 1.0):
        return rate
    k = int(round(rate["value"] * n))
    _, lo, hi = wilson(k, n)
    return {"value": rate["value"], "lo": min(rate["lo"], lo), "hi": max(rate["hi"], hi)}
```

Create `src/cellsurface_sorting_hat/calibration/stability.py` with exactly this content:

```python
"""Stability of calls between two runs (for example two strains of one species).

Proteins are matched by sequence sha256 (``proteins.tsv.gz``), not by ID. Only proteins whose
sequence is identical in both runs are compared; the others are counted. A difference in calls for
an identical sequence comes from the taxon, the status or a module, not from the sequence.
"""

import csv
import gzip
from collections import Counter, defaultdict


def _read(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _calls(out_dir):
    calls = defaultdict(dict)
    for r in _read(f"{out_dir}/calls.long.tsv.gz"):
        calls[r["protein"]][(r["call"], r["variant"])] = r["value"]
    return calls


def compare_runs(out_a, out_b):
    """Return ``{"n_a", "n_b", "n_identical_sequences", "per_call": {(call, variant): Counter}}``."""
    pa = {r["id"]: r["sha256"] for r in _read(f"{out_a}/proteins.tsv.gz")}
    pb = {r["id"]: r["sha256"] for r in _read(f"{out_b}/proteins.tsv.gz")}
    by_sha_b = defaultdict(list)
    for pid, sha in pb.items():
        by_sha_b[sha].append(pid)
    ca, cb = _calls(out_a), _calls(out_b)
    per_call, n_identical, unique = defaultdict(Counter), 0, set()
    for pid, sha in pa.items():
        if sha not in by_sha_b:
            continue
        n_identical += 1
        unique.add(sha)
        other = sorted(by_sha_b[sha])[0]  # identical sequences have identical calls inside one run
        for key, value in ca[pid].items():
            per_call[key]["agree" if cb[other].get(key) == value else "differ"] += 1
    return {
        "n_a": len(pa),
        "n_b": len(pb),
        "n_identical_sequences": n_identical,  # proteins of run A, counted by ID
        "n_unique_identical_sequences": len(unique),
        "per_call": dict(per_call),
    }
```

Create `data/sorting_hat/phasec_set_species.tsv` with exactly this content:

```text
set_key	scientific_name
S1:Scer_SGD	Saccharomyces cerevisiae
S1:Calb_CGD	Candida albicans
S3-Eurotiomycetes:Afum_ASPFU	Aspergillus fumigatus
S3-Eurotiomycetes:Anid_EMENI	Aspergillus nidulans
S3-Basidiomycota:Cneo_H99_GOA	Cryptococcus neoformans
S3-Basidiomycota:Umay_MYCMD	Ustilago maydis
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_phasec.py tests/cellsurface_sorting_hat/calibration/test_stability.py -q`
Expected: `10 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/calibration/phasec.py src/cellsurface_sorting_hat/calibration/stability.py data/sorting_hat/phasec_set_species.tsv tests/cellsurface_sorting_hat/calibration/test_phasec.py tests/cellsurface_sorting_hat/calibration/test_stability.py
git commit -m "feat(sorting-hat): per-species Phase C converter and strain stability comparison

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 8: Proteome provenance

**Files:**
- Create: `src/cellsurface_sorting_hat/proteomes.py`
- Test: `tests/cellsurface_sorting_hat/test_proteomes.py`

**Interfaces:**
- Consumes: `read_fasta`, `OK` (Plan 1); `sha256_file` (Task 1); `write_atomic` (Plan 1).
- Produces: `write_provenance(path, name, source, fasta, retrieved, expected_min, expected_max) -> dict`, `ProteomeError`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_proteomes.py` with exactly this content:

```python
import json

import pytest

from cellsurface_sorting_hat.proteomes import ProteomeError, write_provenance


def test_provenance_checks_the_count_and_records_the_digest(tmp_path):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n>B\nMKT\n>C\nMKS\n")
    rec = write_provenance(
        tmp_path / "prov.json", "toy", "https://example.org/x", fasta, "2026-10-05", 2, 5
    )
    assert (rec["n_proteins"], rec["n_unique_sequences"], rec["n_invalid"]) == (3, 2, 0)
    assert json.loads((tmp_path / "prov.json").read_text())["sha256"] == rec["sha256"]
    with pytest.raises(ProteomeError, match="expected 10 to 20"):
        write_provenance(tmp_path / "p2.json", "toy", "x", fasta, "2026-10-05", 10, 20)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_proteomes.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.proteomes`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/proteomes.py` with exactly this content:

```python
"""Provenance record for a proteome FASTA that a run uses."""

import json
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.fasta import OK, read_fasta
from cellsurface_sorting_hat.modules.base import sha256_file


class ProteomeError(ValueError):
    """The proteome does not look like the one that was expected."""


def write_provenance(path, name, source, fasta, retrieved, expected_min, expected_max):
    """Check the FASTA and write ``path`` (JSON). ``source`` is a URL or a file path; ``retrieved``
    is the date as text. The protein count must lie in ``[expected_min, expected_max]``."""
    proteins = read_fasta(fasta)
    n = len(proteins)
    if not expected_min <= n <= expected_max:
        raise ProteomeError(f"{name}: {n} proteins, expected {expected_min} to {expected_max}")
    record = {
        "name": name,
        "source": source,
        "retrieved": retrieved,
        "fasta": str(Path(fasta).resolve()),
        "sha256": sha256_file(fasta),
        "n_proteins": n,
        "n_invalid": sum(p.state != OK for p in proteins),
        "n_unique_sequences": len({p.sha256 for p in proteins}),
        "first_ids": [p.id for p in proteins[:3]],
    }
    write_atomic(path, (json.dumps(record, indent=2, sort_keys=True) + "\n").encode())
    return record
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_proteomes.py -q`
Expected: `1 passed`

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --files $(git diff --name-only --cached; git ls-files -o --exclude-standard src tests) 2>/dev/null || (ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat)
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/proteomes.py tests/cellsurface_sorting_hat/test_proteomes.py
git commit -m "feat(sorting-hat): proteome provenance record

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 9: Module command and calibration command

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/cli.py`, `src/cellsurface_sorting_hat/calibration/panel.py`, `src/cellsurface_sorting_hat/calibration/cli.py`
- Modify: `pyproject.toml`, `tests/surface_glyco/test_package.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_module_cli.py`, `tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py`

**Interfaces:**
- Consumes: everything from Tasks 1 to 8; `read_taxon_map`, `assign_taxa`, `RunError`, `InputError` (Plan 1 cli).
- Produces the command `cellsurface_sorting_hat_module {signalp,pfam,repeat02,repeat14,allergen,antigen,cys,expression,tm}` (exit 0 or 2; refuses a table without a usable result; `partial` and `unavailable` run states) and the command `cellsurface_sorting_hat_calibrate {phasec,truth,allergen-lso,pfam-specificity,panel}`. `truth` requires `--leakage {none,partial,tuned_on_truth,in_reference,unknown}`; anything but `none` caps the status at `smoke`; it records the sensitivity bound with not-assessable positives counted as missed. `panel` is report-only (verdicts `agree`, `disagree`, `not_in_run`, `known_miss`, `excluded`, `excluded_leakage`).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_module_cli.py` with exactly this content:

```python
"""The module command end to end on tiny inputs."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.modules.cli import main
from cellsurface_sorting_hat.modules.pfam import FAMILY_COLUMNS

FASTA = ">XP_1 a\nMKTAYIAKQRQ\n>XP_2 b\nMNLLPQWERT\n>BAD\nMK*T\n"

OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"


def write_domtbl(path, body):
    path.write_text(OPTIONS + body + "# [ok]\n")


def read(workdir, name):
    with gzip.open(workdir / "modules" / f"{name}.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "p.faa"
    path.write_text(FASTA)
    return path


def test_signalp_command_writes_the_r0_module(tmp_path, fasta, capsys):
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# h\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nXP_1 x\tSP\t0.1\t0.9\t\nXP_2 y\tOTHER\t0.9\t0.1\t\n"
    )
    wd = tmp_path / "wd"
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(wd),
            "--results",
            str(res),
            "--signalp-version",
            "6.0h-gpu",
        ]
    )
    assert code == 0 and json.loads(capsys.readouterr().out)["module"] == "step1_rule@R0"
    rows = read(wd, "step1_rule@R0")
    assert (rows["XP_1"]["call"], rows["XP_2"]["call"], rows["BAD"]["state"]) == (
        "called",
        "not_called",
        "na_invalid",
    )
    run = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert run["tools"] == {"signalp": "6.0h-gpu"} and run["params"]["rule"] == "R0"


def test_pfam_command_writes_both_modules_and_uses_the_recorded_pfam_digest(tmp_path, fasta):
    table = tmp_path / "f.tsv"
    row = dict.fromkeys(FAMILY_COLUMNS, "")
    row.update(
        pfam_acc="PF05730",
        name="CFEM",
        module="pfam_adhesion",
        **{"class": "2b-i"},
        active="yes",
        active_by="owner",
        active_date="2026-10-05",
    )
    row2 = dict(row, pfam_acc="PF16541", name="AltA1", module="pfam_allergen")
    table.write_text(
        "\n".join(
            "\t".join(r[c] for c in FAMILY_COLUMNS)
            for r in [dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True)), row, row2]
        )
        + "\n"
    )
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    wd = tmp_path / "wd"
    assert (
        main(
            [
                "pfam",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--domtbl",
                str(dom),
                "--family-table",
                str(table),
                "--hmmer-version",
                "3.4",
                "--pfam-release",
                "38.2",
                "--pfam-sha256",
                "abc",
            ]
        )
        == 0
    )
    assert read(wd, "pfam_adhesion")["XP_1"]["hit"] == "1"
    assert read(wd, "pfam_allergen")["XP_1"]["hit"] == "0"
    rec = json.loads((wd / "modules" / "pfam_adhesion.json").read_text())
    assert rec["artefact_hash"].startswith("abc:") and rec["params"]["families"] == ["PF05730"]


def test_repeat_command(tmp_path, fasta):
    tab = tmp_path / "r.tsv"
    tab.write_text(
        "protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t5\t3.0\t0.5\nXP_2\t0\t0\t0\n"
    )
    wd = tmp_path / "wd"
    assert main(["repeat02", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "repeat02")["XP_1"]["call"] == "called"


def test_allergen_command(tmp_path, fasta):
    blast = tmp_path / "b.tsv"
    blast.write_text("XP_2\tAsp_f_1.0101|11\t82.0\t90\t120\t100\t150\t1e-20\n")
    ref = tmp_path / "a.faa"
    ref.write_text(">Asp_f_1.0101|11\nMKT\n")
    wd = tmp_path / "wd"
    assert (
        main(
            [
                "allergen",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--blast",
                str(blast),
                "--allergen-fasta",
                str(ref),
                "--blast-version",
                "2.16.0+",
            ]
        )
        == 0
    )
    r = read(wd, "allergen_homology")["XP_2"]
    assert (r["identity"], r["coverage"], r["allergen_name"]) == ("82.00", "90.0", "Asp_f_1.0101")


def test_antigen_command_needs_a_taxon_and_marks_other_taxa_not_applicable(tmp_path, fasta, capsys):
    rank = tmp_path / "r.tsv"
    rank.write_text(
        "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
        "CIMG_1-t26_1-p1\t5\t0.05\t2\t1\t1\t0\n"
    )
    pmap = tmp_path / "m.tsv"
    pmap.write_text("protein_id\tgene_id\tproduct\tlength\nXP_1\tCIMG_1\tp\t11\n")
    wd = tmp_path / "wd"
    base = [
        "antigen",
        "--fasta",
        str(fasta),
        "--workdir",
        str(wd),
        "--ranking",
        str(rank),
        "--protein-map",
        str(pmap),
    ]
    assert main(base) == 2 and "give --taxon or --taxon-map" in capsys.readouterr().err
    assert main(base + ["--taxon", "746128"]) == 0
    assert read(wd, "antigen_lookup")["XP_1"]["state"] == "not_applicable"
    assert main(base + ["--taxon", "246410"]) == 0
    assert read(wd, "antigen_lookup")["XP_1"]["percentile"] == "0.05"


def test_tm_command(tmp_path, fasta):
    tab = tmp_path / "t.tsv"
    tab.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t0\t0\t0\to\nXP_2\t10\t20\t5\t1\ti5-27o\n"
    )
    wd = tmp_path / "wd"
    assert main(["tm", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "tm")["XP_2"]["n_tm"] == "1"


def test_input_errors_exit_2(tmp_path, fasta, capsys):
    assert (
        main(
            [
                "signalp",
                "--fasta",
                str(fasta),
                "--workdir",
                str(tmp_path),
                "--results",
                str(tmp_path / "none"),
                "--signalp-version",
                "6",
            ]
        )
        == 2
    )
    assert "cellsurface_sorting_hat_module: error" in capsys.readouterr().err


def test_a_result_file_whose_ids_do_not_match_the_fasta_is_refused(tmp_path, fasta, capsys):
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# h\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nsp|Q1|X\tSP\t0.1\t0.9\t\n"
    )
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--results",
            str(res),
            "--signalp-version",
            "6.0h-gpu",
        ]
    )
    assert code == 2 and "no protein has a result" in capsys.readouterr().err
    assert not (tmp_path / "wd" / "modules").exists()  # nothing was written


def test_the_signalp_module_identity_does_not_depend_on_the_proteome(tmp_path, fasta):
    """A measurement of rule R0 belongs to the tool version, so it must stay valid for another proteome."""
    other = tmp_path / "o.faa"
    other.write_text(">Z1\nMKTAYI\n")
    ids = {}
    for name, f, pid in (("a", fasta, "XP_1"), ("b", other, "Z1")):
        res = tmp_path / f"{name}.txt"
        res.write_text(
            f"# h\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n{pid}\tSP\t0.1\t0.9\t\n"
        )
        wd = tmp_path / name
        assert (
            main(
                [
                    "signalp",
                    "--fasta",
                    str(f),
                    "--workdir",
                    str(wd),
                    "--results",
                    str(res),
                    "--signalp-version",
                    "6.0h-gpu",
                ]
            )
            == 0
        )
        ids[name] = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert (
        ids["a"]["artefact_hash"] == ids["b"]["artefact_hash"]
        and ids["a"]["params_hash"] == ids["b"]["params_hash"]
    )
    changed = tmp_path / "c"
    assert (
        main(
            [
                "signalp",
                "--fasta",
                str(fasta),
                "--workdir",
                str(changed),
                "--results",
                str(tmp_path / "a.txt"),
                "--signalp-version",
                "6.0i-gpu",
            ]
        )
        == 0
    )
    assert (
        json.loads((changed / "modules" / "step1_rule@R0.json").read_text())["artefact_hash"]
        != ids["a"]["artefact_hash"]
    )


def test_pfam_command_applies_the_no_tm_condition_from_the_tm_table(tmp_path, fasta):
    table = tmp_path / "f.tsv"
    row = dict.fromkeys(FAMILY_COLUMNS, "")
    row.update(
        pfam_acc="PF05730",
        name="CFEM",
        module="pfam_adhesion",
        **{"class": "2b-i"},
        second_condition="no_tm",
        active="yes",
        active_by="owner",
        active_date="2026-10-05",
    )
    header = dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True))
    table.write_text(
        "\n".join("\t".join(r[c] for c in FAMILY_COLUMNS) for r in [header, row]) + "\n"
    )
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom,
        "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
        "XP_2 - 10 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n",
    )
    tm = tmp_path / "t.tsv"
    tm.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t25\t20\t1\ti5-27o\nXP_2\t10\t150\t20\t7\ti5-27o40-62i70-92o\n"
    )
    wd = tmp_path / "wd"
    assert main(["tm", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tm)]) == 0
    assert (
        main(
            [
                "pfam",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--domtbl",
                str(dom),
                "--family-table",
                str(table),
                "--hmmer-version",
                "3.4",
                "--pfam-release",
                "38.2",
                "--pfam-sha256",
                "abc",
                "--tm-module",
                "tm",
            ]
        )
        == 0
    )
    rows = read(wd, "pfam_adhesion")
    assert (rows["XP_1"]["hit"], rows["XP_2"]["hit"]) == (
        "1",
        "0",
    )  # XP_1: only the signal peptide read as a helix; XP_2: helices after residue 35


def _family_table(tmp_path, active):
    table = tmp_path / "f.tsv"
    base = dict.fromkeys(FAMILY_COLUMNS, "")
    base.update(
        pfam_acc="PF05730", name="CFEM", module="pfam_adhesion", **{"class": "2b-i"}, active="no"
    )
    if active:
        base.update(active="yes", active_by="owner", active_date="2026-10-05")
    header = dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True))
    table.write_text(
        "\n".join("\t".join(r[c] for c in FAMILY_COLUMNS) for r in [header, base]) + "\n"
    )
    return table


def _pfam_args(tmp_path, fasta, table, dom, wd):
    return [
        "pfam",
        "--fasta",
        str(fasta),
        "--workdir",
        str(wd),
        "--domtbl",
        str(dom),
        "--family-table",
        str(table),
        "--hmmer-version",
        "3.4",
        "--pfam-release",
        "38.2",
        "--pfam-sha256",
        "abc",
    ]


def test_pfam_command_writes_unavailable_modules_when_no_family_is_active(tmp_path, fasta):
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    wd = tmp_path / "wd"
    assert main(_pfam_args(tmp_path, fasta, _family_table(tmp_path, active=False), dom, wd)) == 0
    rec = json.loads((wd / "modules" / "pfam_adhesion.json").read_text())
    assert rec["run_state"] == "unavailable" and "no active family" in rec["note"]
    assert {r["state"] for r in read(wd, "pfam_adhesion").values() if r["id"] != "BAD"} == {
        "unavailable"
    }


def test_pfam_command_refuses_a_domain_table_from_another_proteome(tmp_path, fasta, capsys):
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "OTHER1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    code = main(
        _pfam_args(tmp_path, fasta, _family_table(tmp_path, active=True), dom, tmp_path / "wd")
    )
    assert code == 2 and "not in the FASTA" in capsys.readouterr().err


def test_pfam_command_refuses_an_unfinished_domain_table(tmp_path, fasta, capsys):
    dom = tmp_path / "d.domtbl"
    dom.write_text(
        OPTIONS
        + "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    code = main(
        _pfam_args(tmp_path, fasta, _family_table(tmp_path, active=True), dom, tmp_path / "wd")
    )
    assert code == 2 and "[ok]" in capsys.readouterr().err


def test_allergen_command_refuses_a_blast_table_from_another_proteome(tmp_path, fasta, capsys):
    blast = tmp_path / "b.tsv"
    blast.write_text("OTHER\tAsp_f_1.0101|11\t82.0\t90\t120\t100\t150\t1e-20\n")
    ref = tmp_path / "a.faa"
    ref.write_text(">Asp_f_1.0101|11\nMKT\n")
    code = main(
        [
            "allergen",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--blast",
            str(blast),
            "--allergen-fasta",
            str(ref),
            "--blast-version",
            "2.14.0+",
        ]
    )
    assert code == 2 and "not in the FASTA" in capsys.readouterr().err


def test_partial_results_set_the_run_state_partial(tmp_path, fasta):
    tab = tmp_path / "r.tsv"
    tab.write_text("protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t5\t3.0\t0.5\n")
    wd = tmp_path / "wd"
    long_fasta = tmp_path / "long.faa"
    long_fasta.write_text(">XP_1\n" + "M" * 100 + "\n>XP_2\n" + "M" * 100 + "\n")
    assert (
        main(["repeat02", "--fasta", str(long_fasta), "--workdir", str(wd), "--table", str(tab)])
        == 0
    )
    rec = json.loads((wd / "modules" / "repeat02.json").read_text())
    assert rec["run_state"] == "partial" and "1 protein(s) have no result" in rec["note"]


def test_a_lookup_that_matches_no_applicable_protein_is_refused(tmp_path, fasta, capsys):
    cand = tmp_path / "c.tsv"
    cand.write_text("protein_id\ttier\tcys_frac\nCIMG_9\tcys_rich_sp_unassigned\t0.1\n")
    code = main(
        [
            "cys",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--taxon",
            "246410",
            "--candidates",
            str(cand),
        ]
    )
    assert code == 2 and "no protein has a result" in capsys.readouterr().err


def test_signalp_file_from_a_non_eukarya_run_is_refused(tmp_path, fasta, capsys):
    res = tmp_path / "p.txt"
    res.write_text(
        "# SignalP-6.0\tOrganism: Other\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nXP_1\tSP\t0.1\t0.9\t\n"
    )
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--results",
            str(res),
            "--signalp-version",
            "6",
        ]
    )
    assert code == 2 and "not a Eukarya run" in capsys.readouterr().err
```

Create `tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py` with exactly this content:

```python
import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.calibration.cli import main
from cellsurface_sorting_hat.calibration.panel import panel_check
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source

OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"


def write_domtbl(path, body):
    path.write_text(OPTIONS + body + "# [ok]\n")


def write_calls(path, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for p, c, v, val in rows:
            w.writerow([p, c, v, val, "unvalidated", "", ""])


def write_truth(path, rows):
    path.write_text("id\tlabel\tcluster\n" + "".join(f"{i}\t{y}\t{c}\n" for i, y, c in rows))


def test_truth_command_writes_an_entry_with_sensitivity_and_specificity(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("allergen_homology", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(40)]
    calls = [
        (i, "iuis_allergen_homolog", "", "called" if k < 12 or 20 <= k < 24 else "not_called")
        for k, i in enumerate(ids)
    ]
    calls.append(("U1", "iuis_allergen_homolog", "", "not_assessable"))
    write_calls(tmp_path / "c.tsv.gz", calls)
    write_truth(
        tmp_path / "t.tsv",
        [(i, 1 if k < 20 else 0, f"c{k}") for k, i in enumerate(ids)]
        + [("U1", 1, "cu"), ("GONE", 0, "cg")],
    )
    code = main(
        [
            "truth",
            "--workdir",
            str(tmp_path),
            "--module",
            "allergen_homology",
            "--calls-long",
            str(tmp_path / "c.tsv.gz"),
            "--call",
            "iuis_allergen_homolog",
            "--truth",
            str(tmp_path / "t.tsv"),
            "--calibration-set",
            "toy",
            "--leakage",
            "none",
            "--taxa",
            "746128",
            "--n-boot",
            "200",
        ]
    )
    assert code == 0
    entry = load_status_source(tmp_path / "status" / "allergen_homology.json").entries[0]
    m = entry.measure
    assert (m["n_pos"], m["n_neg"]) == (20, 20)
    assert m["sensitivity"]["value"] == pytest.approx(12 / 20) and m["specificity"][
        "value"
    ] == pytest.approx(16 / 20)
    assert "truth rows without a call: 1" in m["notes"] and "not assessable: 1" in m["notes"]
    assert "sensitivity if not assessable positives count as missed: 0.571" in m["notes"]
    assert entry.status == "smoke"  # 20 positives but a wide interval


def test_a_second_set_is_added_and_the_same_set_is_replaced(tmp_path):
    write_module(tmp_path, ModuleSpec("repeat02", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(10)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "tandem_repeat_protein", "", "called" if k < 5 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 5 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(name, taxon):
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "repeat02",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "tandem_repeat_protein",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--leakage",
                "none",
                "--taxa",
                str(taxon),
                "--n-boot",
                "50",
            ]
        )

    assert run("setA", 4932) == 0 and run("setB", 5476) == 0 and run("setA", 4932) == 0
    sets = [
        e.measure["calibration_set"]
        for e in load_status_source(tmp_path / "status" / "repeat02.json").entries
    ]
    assert sorted(sets) == ["setA", "setB"]


def test_truth_with_no_matching_protein_is_an_error(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("repeat02", "1"), [], [{"id": "A", "state": "ok"}])
    write_calls(tmp_path / "c.tsv.gz", [("X", "tandem_repeat_protein", "", "called")])
    write_truth(tmp_path / "t.tsv", [("Y", 1, "c")])
    assert (
        main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "repeat02",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "tandem_repeat_protein",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                "s",
                "--leakage",
                "none",
                "--taxa",
                "1",
            ]
        )
        == 2
    )
    assert "no truth protein has a call" in capsys.readouterr().err


def test_phasec_command_resolves_species_names_to_taxa(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(
        "".join(
            f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n"
            for t, n in [
                (4932, "Saccharomyces cerevisiae"),
                (5476, "Candida albicans"),
                (746128, "Aspergillus fumigatus"),
                (162425, "Aspergillus nidulans"),
                (5207, "Cryptococcus neoformans"),
                (5270, "Ustilago maydis"),
            ]
        )
    )
    sp = tmp_path / "sets.tsv"
    sp.write_text(
        "set_key\tscientific_name\nS1:Scer_SGD\tSaccharomyces cerevisiae\nS1:Calb_CGD\tCandida albicans\n"
        "S3-Eurotiomycetes:Afum_ASPFU\tAspergillus fumigatus\nS3-Eurotiomycetes:Anid_EMENI\tAspergillus nidulans\n"
        "S3-Basidiomycota:Cneo_H99_GOA\tCryptococcus neoformans\nS3-Basidiomycota:Umay_MYCMD\tUstilago maydis\n"
    )

    def cell(r, f):
        return {
            "recall": {"value": r, "lo": r - 0.05, "hi": r + 0.05},
            "fpr": {"value": f, "lo": f / 2, "hi": f * 2},
        }

    def ts(label, pos, c):
        return {
            "label": label,
            "n_direct_positives": pos,
            "truth": {
                "direct": {
                    "n": {"all": {"pos": pos, "neg": 100}},
                    "metrics": {"all": {"V-go": {"R0": c}}},
                }
            },
        }

    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "test_sets": {
                    "S1:Scer_SGD": ts("estimate", 232, cell(0.6, 0.04)),
                    "S1:Calb_CGD": ts("smoke test", 153, cell(0.5, 0.04)),
                    "S3-Eurotiomycetes:Afum_ASPFU": ts("smoke test", 19, cell(0.9, 0.01)),
                    "S3-Eurotiomycetes:Anid_EMENI": ts("estimate", 109, cell(0.7, 0.01)),
                    "S3-Basidiomycota:Cneo_H99_GOA": ts("smoke test", 7, cell(0.9, 0.08)),
                    "S3-Basidiomycota:Umay_MYCMD": ts("smoke test", 9, cell(0.9, 0.08)),
                }
            }
        )
    )
    write_module(
        tmp_path / "wd", ModuleSpec("step1_rule@R0", "1"), [], [{"id": "A", "state": "ok"}]
    )
    assert (
        main(
            [
                "phasec",
                "--workdir",
                str(tmp_path / "wd"),
                "--metrics",
                str(metrics),
                "--set-species",
                str(sp),
                "--names-dmp",
                str(names),
            ]
        )
        == 0
    )
    entries = load_status_source(tmp_path / "wd" / "status" / "step1_rule@R0.json").entries
    assert {e.taxa for e in entries} == {(4932,), (5476,), (746128,), (162425,), (5207,), (5270,)}
    assert {e.measure["calibration_set"]: e.status for e in entries}[
        "S3-Basidiomycota:Umay_MYCMD"
    ] == "smoke"


def test_a_changed_module_identity_drops_the_old_status_entries(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("repeat02", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(10)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "tandem_repeat_protein", "", "called" if k < 5 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 5 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(name, taxon):
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "repeat02",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "tandem_repeat_protein",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--taxa",
                str(taxon),
                "--leakage",
                "none",
                "--n-boot",
                "50",
            ]
        )

    assert run("setA", 4932) == 0
    write_module(
        tmp_path, ModuleSpec("repeat02", "2", {"a": 2}), [], [{"id": "A", "state": "ok"}]
    )  # new version
    assert run("setB", 5476) == 0
    assert "module identity changed" in capsys.readouterr().err
    entries = load_status_source(tmp_path / "status" / "repeat02.json").entries
    assert [e.measure["calibration_set"] for e in entries] == [
        "setB"
    ]  # setA was measured on version 1


def test_panel_check_counts_agreement_and_leaves_known_misses_out(tmp_path):
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            ("A", "tandem_repeat_protein", "", "called"),
            ("B", "tandem_repeat_protein", "", "not_called"),
            ("C", "tandem_repeat_protein", "", "called"),
        ],
    )
    panel = tmp_path / "p.tsv"
    panel.write_text(
        "protein\tcall\tvariant\texpected\tsource\nA\ttandem_repeat_protein\t\tcalled\tx\nB\ttandem_repeat_protein\t\tcalled\tx\n"
        "C\ttandem_repeat_protein\t\tknown_miss\tx\nZ\ttandem_repeat_protein\t\tcalled\tx\n"
    )
    rows, summary = panel_check(tmp_path / "c.tsv.gz", panel)
    assert summary == {"agree": 1, "disagree": 1, "known_miss": 1, "not_in_run": 1}
    assert [r["observed"] for r in rows] == ["called", "not_called", "called", "missing"]


def test_allergen_lso_command_prints_recall_per_rule(tmp_path, capsys):
    fasta = tmp_path / "a.faa"
    fasta.write_text(">Asp_f_1.0101|1\nM\n>Asp_n_1.0101|2\nM\n>Alt_a_1.0101|3\nM\n")
    blast = tmp_path / "b.tsv"
    blast.write_text(
        "Asp_f_1.0101|1\tAsp_f_1.0101|1\t100\t100\t100\t100\t500\t0\n"
        "Asp_f_1.0101|1\tAsp_n_1.0101|2\t60\t90\t100\t100\t100\t1e-20\n"
        "Asp_n_1.0101|2\tAsp_n_1.0101|2\t100\t100\t100\t100\t500\t0\n"
        "Asp_n_1.0101|2\tAsp_f_1.0101|1\t60\t90\t100\t100\t100\t1e-20\n"
        "Alt_a_1.0101|3\tAlt_a_1.0101|3\t100\t100\t100\t100\t500\t0\n"
    )
    assert main(["allergen-lso", "--blast", str(blast), "--allergen-fasta", str(fasta)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "sequences\t3\tspecies\t3"
    assert out[1] == "iuis_allergen_similarity\t2/3" and out[2] == "iuis_allergen_homolog\t0/3"


def test_leakage_other_than_none_caps_an_estimate_at_smoke(tmp_path):
    write_module(tmp_path, ModuleSpec("antigen_lookup", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(60)]
    # 30 positives, all called; 30 negatives, none called: a tight interval that would be an estimate
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "cocci_specificity_rank_top15", "", "called" if k < 30 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 30 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(leakage, name):
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "antigen_lookup",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "cocci_specificity_rank_top15",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--taxa",
                "246410" if name == "a" else "5476",
                "--leakage",
                leakage,
                "--n-boot",
                "100",
            ]
        )

    assert run("none", "a") == 0 and run("tuned_on_truth", "b") == 0
    status = {
        e.measure["calibration_set"]: (e.status, e.measure["notes"])
        for e in load_status_source(tmp_path / "status" / "antigen_lookup.json").entries
    }
    assert status["a"][0] == "estimated" and status["b"][0] == "smoke"
    assert "leakage: tuned_on_truth" in status["b"][1]


def test_pfam_specificity_command_lists_non_member_hits_with_their_helices(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n>B\nMKT\n>C\nMKT\n>D\nMKT\n")
    dom = tmp_path / "d.domtbl"
    row = "{} - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    write_domtbl(dom, row.format("A") + row.format("B") + row.format("D"))
    members = tmp_path / "m.tsv"
    members.write_text("pfam_acc\tprotein_id\nPF05730\tA\nPF05730\tC\nPF01185\tD\n")
    tm = tmp_path / "t.tsv"
    tm.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nB\t11\t150\t20\t7\to5-27i\n"
    )
    assert (
        main(
            [
                "pfam-specificity",
                "--family",
                "PF05730",
                "--domtbl",
                str(dom),
                "--members",
                str(members),
                "--universe-fasta",
                str(fasta),
                "--tm-table",
                str(tm),
            ]
        )
        == 0
    )
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "family\tPF05730\thits\t3\tmembers\t2"
    assert "tp\t1" in out and "fp\t2" in out and "fn\t1" in out
    assert (
        "nonmember_hit\tB\tn_tm=7" in out
        and "nonmember_hit\tD\tn_tm=NA" in out
        and "missed_member\tC" in out
    )


def test_pfam_specificity_refuses_members_outside_the_proteome(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n")
    dom = tmp_path / "d.domtbl"
    write_domtbl(dom, "")
    members = tmp_path / "m.tsv"
    members.write_text("pfam_acc\tprotein_id\nPF05730\tZ\n")
    assert (
        main(
            [
                "pfam-specificity",
                "--family",
                "PF05730",
                "--domtbl",
                str(dom),
                "--members",
                str(members),
                "--universe-fasta",
                str(fasta),
            ]
        )
        == 2
    )
    assert "not in the proteome" in capsys.readouterr().err


def test_panel_marks_proteins_seen_in_the_reference_as_excluded(tmp_path):
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            ("A", "iuis_allergen_homolog", "", "called"),
            ("B", "iuis_allergen_homolog", "", "called"),
        ],
    )
    panel = tmp_path / "p.tsv"
    panel.write_text(
        "protein\tcall\tvariant\texpected\tsource\ttuning\nA\tiuis_allergen_homolog\t\tcalled\tx\tin_reference\n"
        "B\tiuis_allergen_homolog\t\tcalled\tx\t\n"
    )
    _, summary = panel_check(tmp_path / "c.tsv.gz", panel)
    assert summary == {"excluded_leakage": 1, "agree": 1}
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_module_cli.py tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.cli`.

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/cli.py` with exactly this content:

```python
"""Command ``cellsurface_sorting_hat_module``: turn one tool's output into a module table."""

import argparse
import json
import sys
from pathlib import Path

from cellsurface_sorting_hat.cli import InputError, RunError, assign_taxa, read_taxon_map
from cellsurface_sorting_hat.fasta import FastaError, read_fasta
from cellsurface_sorting_hat.modules import allergen, lookups, pfam, repeats, signalp
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module

APPLICABLE_RS = {
    246410
}  # C. immitis RS (taxon_id in analysis/step1_compare/phasec/tc_taxon_clades.tsv)


def _common(p):
    p.add_argument("--fasta", required=True)
    p.add_argument("--workdir", required=True)


def _taxa_args(p):
    p.add_argument("--taxon", type=int)
    p.add_argument("--taxon-map")
    p.add_argument("--applicable-taxa", type=int, nargs="+", default=sorted(APPLICABLE_RS))


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_module")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("signalp", help="SignalP 6 prediction_results.txt -> step1_rule@R0")
    _common(p)
    p.add_argument("--results", required=True)
    p.add_argument("--signalp-version", required=True, help="for example 6.0h-gpu (recorded)")
    p.add_argument("--signalp-mode", default="fast")

    p = sub.add_parser("pfam", help="hmmsearch --domtblout -> pfam_adhesion and pfam_allergen")
    _common(p)
    p.add_argument("--domtbl", required=True)
    p.add_argument("--family-table", required=True)
    p.add_argument("--pfam-release", required=True)
    p.add_argument(
        "--pfam-sha256", required=True, help="sha256 of Pfam-A.hmm, from provenance.json"
    )
    p.add_argument("--hmmer-version", required=True, help="for example 3.4 (recorded)")
    p.add_argument(
        "--sp-module", help="step1_rule@R0 table in --workdir, for second_condition=signal_peptide"
    )
    p.add_argument("--tm-module", help="tm table in --workdir, for second_condition=no_tm")

    for name in ("repeat02", "repeat14"):
        p = sub.add_parser(name, help="repeat detector table -> " + name)
        _common(p)
        p.add_argument("--table", required=True)
        p.add_argument("--min-coverage", type=float, default=repeats.DEFAULT_MIN_COVERAGE)
        p.add_argument("--min-copies", type=float, default=repeats.DEFAULT_MIN_COPIES)
        p.add_argument(
            "--script", help="detector script; its sha256 is part of the module identity"
        )

    p = sub.add_parser("allergen", help="BLASTP results against the IUIS fungal allergens")
    _common(p)
    p.add_argument("--blast", required=True, help="outfmt 6 with: " + allergen.BLAST_FIELDS)
    p.add_argument("--allergen-fasta", required=True)
    p.add_argument("--blast-version", required=True)
    p.add_argument("--evalue", default="1")
    p.add_argument("--seg", default="no", help="BLAST low-complexity masking used by the job")

    p = sub.add_parser("antigen", help="antigen ranking lookup (Coccidioides)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--ranking", required=True)
    p.add_argument("--protein-map", required=True)

    p = sub.add_parser("cys", help="Cys-rich tiers (analysis/cys_candidates candidates.tsv.gz)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--candidates", required=True)

    p = sub.add_parser("expression", help="spherule table lookup (Coccidioides)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--table", required=True)
    p.add_argument("--protein-map", required=True)

    p = sub.add_parser("tm", help="TMHMM table -> tm")
    _common(p)
    p.add_argument("--table", required=True)
    p.add_argument("--tmhmm-version", default="2.0c")
    return ap


def _taxa(args, proteins):
    tmap = read_taxon_map(args.taxon_map) if args.taxon_map else {}
    if args.taxon is None and not tmap:
        raise RunError("give --taxon or --taxon-map")
    return assign_taxa(proteins, args.taxon, tmap)


def _write(workdir, spec, columns, rows):
    """Write the module. Refuse a table without one usable result: every protein ``error`` (an ID
    mismatch) or no applicable protein ``ok`` (a lookup that matched nothing). Some ``error`` rows
    make the run state ``partial``."""
    counted = [r for r in rows if r["state"] not in ("na_invalid", "not_applicable")]
    if counted and not any(r["state"] == "ok" for r in counted):
        states = sorted({r["state"] for r in counted})
        raise RunError(
            f"{spec.name}: no protein has a result (states: {', '.join(states)}); "
            "check that the IDs agree"
        )
    errors = sum(1 for r in rows if r["state"] == "error")
    state = "partial" if errors else "ok"
    note = f"{errors} protein(s) have no result" if errors else ""
    return write_module(workdir, spec, columns, rows, run_state=state, note=note)


def _family_digest(families):
    """Hash of what decides a call (accession, module, condition, active), not of the free text."""
    import hashlib

    rows = sorted((f.pfam_acc, f.module, f.second_condition, f.active) for f in families)
    return hashlib.sha256(repr(rows).encode()).hexdigest()


def run(args):
    proteins = read_fasta(args.fasta)
    ids = [p.id for p in proteins]
    w = args.workdir
    if args.cmd == "signalp":
        spec = ModuleSpec(
            "step1_rule@R0",
            "1",
            {"rule": "R0", "mode": args.signalp_mode, "organism": "eukarya"},
            (),
            {"signalp": args.signalp_version},
            artefact_digest=_tool_digest("signalp", args.signalp_version),
        )
        rows = signalp.signalp_rows(proteins, signalp.parse_signalp(args.results))
        return _write(w, spec, signalp.COLUMNS, rows)
    if args.cmd == "pfam":
        families = pfam.load_family_table(args.family_table)
        hits = pfam.parse_domtblout(args.domtbl)
        pfam.check_hit_ids(hits, ids)
        sp_calls = _module_column(w, args.sp_module, "call") if args.sp_module else None
        tm_counts = None
        if args.tm_module:
            raw = _module_column(w, args.tm_module, "n_tm_mature")
            tm_counts = {k: int(v) for k, v in raw.items() if v.isdigit()}
        out = None
        for module in pfam.MODULES:
            active = sorted(f.pfam_acc for f in families if f.module == module and f.active)
            params = {"pfam_release": args.pfam_release, "cut": "ga", "families": active}
            spec = ModuleSpec(
                module,
                "1",
                params,
                (),
                {"hmmer": args.hmmer_version},
                artefact_digest=args.pfam_sha256 + ":" + _family_digest(families),
            )
            if not active:  # no family has passed its specificity test: no domain test was made
                rows = [{"id": p.id, "state": "unavailable"} for p in proteins]
                out = write_module(
                    w, spec, pfam.COLUMNS, rows, run_state="unavailable", note="no active family"
                )
                continue
            rows = pfam.pfam_rows(proteins, hits, families, module, sp_calls, tm_counts)
            out = _write(w, spec, pfam.COLUMNS, rows)
        return out
    if args.cmd in ("repeat02", "repeat14"):
        params = {
            "min_coverage": args.min_coverage,
            "min_copies": args.min_copies,
            "min_len": repeats.DETECTOR_MIN_LEN,
            "script_sha256": _file_digest(args.script) if args.script else "",
        }
        spec = ModuleSpec(args.cmd, "1", params)
        rows = repeats.repeat_rows(
            proteins, repeats.parse_repeat_table(args.table), args.min_coverage, args.min_copies
        )
        return _write(w, spec, repeats.COLUMNS, rows)
    if args.cmd == "allergen":
        best = allergen.parse_blast(args.blast)
        allergen.check_blast_ids(best, ids)
        meta_path = Path(str(args.allergen_fasta) + ".meta.tsv")
        meta = allergen.read_meta(meta_path) if meta_path.exists() else {}
        spec = ModuleSpec(
            "allergen_homology",
            "1",
            {"evalue": args.evalue, "program": "blastp", "seg": args.seg},
            (args.allergen_fasta,),
            {"blast": args.blast_version},
        )
        rows = allergen.allergen_rows(proteins, best, meta)
        return _write(w, spec, allergen.COLUMNS, rows)
    if args.cmd == "tm":
        table, _ = lookups.read_table(args.table, "protein_id")
        spec = ModuleSpec(
            "tm",
            "1",
            {"signal_peptide_window": lookups.SIGNAL_PEPTIDE_WINDOW},
            (),
            {"tmhmm": args.tmhmm_version},
            artefact_digest=_tool_digest("tmhmm", args.tmhmm_version),
        )
        return _write(w, spec, lookups.TM_COLUMNS, lookups.tm_rows(proteins, table))
    taxa = _taxa(args, proteins)
    applicable = set(args.applicable_taxa)
    if args.cmd == "antigen":
        by_gene, several = lookups.ranking_by_gene(args.ranking)
        pmap = lookups.load_protein_map(args.protein_map)
        params = {"applicable_taxa": sorted(applicable), "genes_with_several_ranking_rows": several}
        spec = ModuleSpec("antigen_lookup", "1", params, (args.ranking, args.protein_map))
        rows = lookups.antigen_rows(proteins, taxa, pmap, by_gene, applicable)
        return _write(w, spec, lookups.ANTIGEN_COLUMNS, rows)
    if args.cmd == "cys":
        table, _ = lookups.read_table(args.candidates, "protein_id")
        spec = ModuleSpec(
            "cys_rich", "1", {"applicable_taxa": sorted(applicable)}, (args.candidates,)
        )
        rows = lookups.cys_rows(proteins, taxa, table, applicable)
        return _write(w, spec, lookups.CYS_COLUMNS, rows)
    if args.cmd == "expression":
        table, _ = lookups.read_table(args.table, "gene_id")
        pmap = lookups.load_protein_map(args.protein_map)
        spec = ModuleSpec(
            "expression",
            "1",
            {"applicable_taxa": sorted(applicable)},
            (args.table, args.protein_map),
        )
        rows = lookups.expression_rows(proteins, taxa, pmap, table, applicable)
        return _write(w, spec, lookups.EXPRESSION_COLUMNS, rows)
    raise AssertionError(args.cmd)


def _tool_digest(tool, version):
    """Identity of a tool: the measurement of a rule belongs to the tool version, not to one input."""
    import hashlib

    return hashlib.sha256(f"{tool}:{version}".encode()).hexdigest()


def _module_column(workdir, module, column):
    import csv
    import gzip

    with gzip.open(Path(workdir) / "modules" / f"{module}.tsv.gz", "rt", newline="") as fh:
        return {r["id"]: r.get(column, "") for r in csv.DictReader(fh, delimiter="\t")}


def _file_digest(path):
    from cellsurface_sorting_hat.modules.base import sha256_file

    return sha256_file(path)


def main(argv=None):
    try:
        record = run(build_parser().parse_args(argv))
    except (RunError, InputError, FastaError, ValueError, KeyError, OSError) as err:
        print(f"cellsurface_sorting_hat_module: error: {err}", file=sys.stderr)
        return 2
    print(json.dumps({k: record[k] for k in ("module", "run_state", "n_rows")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

Create `src/cellsurface_sorting_hat/calibration/panel.py` with exactly this content:

```python
"""Report-only check of calls against a panel of proteins with expected values."""

import csv
import gzip
from collections import Counter

EXPECTED = ("called", "not_called", "known_miss", "excluded")
LEAKY = ("in_reference", "tuned", "partial")


def read_calls(path):
    with gzip.open(path, "rt", newline="") as fh:
        return {
            (r["protein"], r["call"], r["variant"]): r["value"]
            for r in csv.DictReader(fh, delimiter="\t")
        }


def panel_check(calls_long, panel_tsv):
    """Compare ``calls_long`` with ``panel_tsv`` (columns protein, call, variant, expected, source; optional column tuning:
    ``in_reference``, ``tuned`` or ``partial`` for a protein that the module, its reference set or its
    cutoffs have seen).

    ``known_miss`` and ``excluded`` rows are listed and not counted as agreement or disagreement.
    Returns ``(rows, summary)``. Nothing here is an acceptance gate.
    """
    calls = read_calls(calls_long)
    rows, summary = [], Counter()
    with open(panel_tsv, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["expected"] not in EXPECTED:
                raise ValueError(f"{r['protein']}: expected must be one of {EXPECTED}")
            observed = calls.get((r["protein"], r["call"], r["variant"]), "missing")
            tuning = (r.get("tuning") or "").strip()
            if r["expected"] in ("known_miss", "excluded"):
                verdict = r["expected"]
            elif tuning in LEAKY:
                verdict = "excluded_leakage"  # seen in tuning or part of the reference set
            elif observed == "missing":
                verdict = "not_in_run"  # the protein is not in the proteome that was run
            else:
                verdict = "agree" if observed == r["expected"] else "disagree"
            summary[verdict] += 1
            rows.append({**r, "observed": observed, "verdict": verdict})
    return rows, dict(summary)
```

Create `src/cellsurface_sorting_hat/calibration/cli.py` with exactly this content:

```python
"""Command ``cellsurface_sorting_hat_calibrate``: write status sources and check panels."""

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

from cellsurface_sorting_hat.calibration import phasec
from cellsurface_sorting_hat.calibration.measure import (
    build_measure,
    make_entry,
    write_status_source,
)
from cellsurface_sorting_hat.calibration.panel import panel_check
from cellsurface_sorting_hat.fasta import read_fasta
from cellsurface_sorting_hat.modules import allergen, pfam
from cellsurface_sorting_hat.taxonomy import TaxonError


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_calibrate")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser(
        "phasec", help="status source for step1_rule@R0 from the Phase C metrics.json"
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--set-species", required=True, help="TSV: set_key, scientific_name")
    p.add_argument("--names-dmp", required=True)

    p = sub.add_parser(
        "truth", help="sensitivity and specificity of one call against a truth table"
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--module", required=True, help="module that receives the status entry")
    p.add_argument(
        "--calls-long", required=True, help="calls.long.tsv.gz of a run on the truth proteins"
    )
    p.add_argument("--call", required=True)
    p.add_argument("--variant", default="")
    p.add_argument("--truth", required=True, help="TSV: id, label (1 or 0), cluster")
    p.add_argument("--calibration-set", required=True)
    p.add_argument("--taxa", type=int, nargs="+", required=True, help="tested taxa (species level)")
    p.add_argument("--notes", default="")
    p.add_argument(
        "--leakage",
        required=True,
        choices=["none", "partial", "tuned_on_truth", "in_reference", "unknown"],
        help="did the truth proteins help to set the rule or its cutoffs? anything but 'none' caps the status at smoke",
    )
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser(
        "allergen-lso", help="leave-species-out recall of the allergen set (sensitivity only)"
    )
    p.add_argument("--blast", required=True, help="allergens against allergens, outfmt 6")
    p.add_argument(
        "--allergen-fasta", required=True, help="the searched FASTA; fixes the denominator"
    )

    p = sub.add_parser(
        "pfam-specificity", help="specificity test of one Pfam family on one proteome"
    )
    p.add_argument("--family", required=True, help="for example PF05730")
    p.add_argument("--domtbl", required=True, help="hmmsearch --domtblout of the proteome")
    p.add_argument(
        "--members",
        required=True,
        help="TSV: pfam_acc, protein_id (members known from curation, not from Pfam)",
    )
    p.add_argument("--universe-fasta", required=True, help="the proteome that was searched")
    p.add_argument(
        "--tm-table",
        help="TMHMM table (protein_id, pred_hel) to show helices beside each non-member hit",
    )

    p = sub.add_parser("panel", help="report-only check of calls against a panel")
    p.add_argument("--calls-long", required=True)
    p.add_argument("--panel", required=True)
    return ap


def _set_taxa(path, names_dmp):
    names = phasec.read_names(names_dmp)
    by_set = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_set.setdefault(r["set_key"], []).append(
                phasec.species_taxid(names, r["scientific_name"])
            )
    return by_set


def _merge_entries(path, new_entries, workdir=None, module=None):
    """Entries already in the status source stay, except those of the same calibration set.

    Old entries are dropped (with a message) when the module identity in the workdir differs from
    the identity in the old file: a measurement of an older version must not carry the new one.
    """
    if not Path(path).exists():
        return list(new_entries)
    old_file = json.loads(Path(path).read_text())
    old = old_file["entries"]
    if workdir is not None:
        now = json.loads((Path(workdir) / "modules" / f"{module}.json").read_text())
        keys = ("version", "params_hash", "artefact_hash")
        if any(old_file.get(k) != now.get(k) for k in keys):
            print(
                f"{path}: module identity changed; {len(old)} old entr(ies) dropped",
                file=sys.stderr,
            )
            return list(new_entries)
    names = {e["measure"]["calibration_set"] for e in new_entries}
    return [e for e in old if e["measure"]["calibration_set"] not in names] + list(new_entries)


def run(args):
    if args.cmd == "phasec":
        entries = phasec.entries_from_phasec(
            args.metrics, _set_taxa(args.set_species, args.names_dmp)
        )
        path = Path(args.workdir) / "status" / "step1_rule@R0.json"
        return write_status_source(
            args.workdir,
            "step1_rule@R0",
            _merge_entries(path, entries, args.workdir, "step1_rule@R0"),
        )
    if args.cmd == "truth":
        with gzip.open(args.calls_long, "rt", newline="") as fh:
            calls = {
                r["protein"]: r["value"]
                for r in csv.DictReader(fh, delimiter="\t")
                if r["call"] == args.call and r["variant"] == args.variant
            }
        y, called, clusters, unmatched, unknown, unknown_pos = [], [], [], 0, 0, 0
        with open(args.truth, encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                value = calls.get(r["id"])
                if value is None:
                    unmatched += 1
                elif value == "not_assessable":
                    unknown += 1
                    unknown_pos += int(r["label"]) == 1
                else:
                    y.append(int(r["label"]))
                    called.append(value == "called")
                    clusters.append(r["cluster"])
        if not y:
            raise ValueError("no truth protein has a call in --calls-long")
        pos_unknown = unknown_pos
        n_pos_called = sum(1 for lab, c in zip(y, called, strict=True) if lab == 1 and c)
        n_pos_all = sum(1 for lab in y if lab == 1) + pos_unknown
        bound = n_pos_called / n_pos_all if n_pos_all else float("nan")
        notes = (
            f"{args.notes} truth rows without a call: {unmatched}; not assessable: {unknown}; "
            f"sensitivity if not assessable positives count as missed: {bound:.3f}"
        ).strip()
        measure = build_measure(
            args.calibration_set,
            str(args.truth),
            y,
            called,
            clusters,
            notes,
            args.n_boot,
            args.seed,
        )
        path = Path(args.workdir) / "status" / f"{args.module}.json"
        measure["notes"] = f"{measure['notes']} leakage: {args.leakage}".strip()
        cap = None if args.leakage == "none" else "smoke"
        entry = make_entry(args.taxa, measure, source=str(args.truth), cap=cap)
        return write_status_source(
            args.workdir, args.module, _merge_entries(path, [entry], args.workdir, args.module)
        )
    if args.cmd == "pfam-specificity":
        hits = {h["target"] for h in pfam.parse_domtblout(args.domtbl) if h["acc"] == args.family}
        with open(args.members, encoding="utf-8-sig", newline="") as fh:
            members = {
                r["protein_id"]
                for r in csv.DictReader(fh, delimiter="\t")
                if r["pfam_acc"] == args.family
            }
        universe = {p.id for p in read_fasta(args.universe_fasta)}
        outside = members - universe
        if outside:
            raise ValueError(
                f"{len(outside)} member(s) are not in the proteome, for example {sorted(outside)[0]}"
            )
        rep = pfam.specificity_report(hits & universe, members, universe)
        print(f"family\t{args.family}\thits\t{len(hits & universe)}\tmembers\t{len(members)}")
        for key in ("tp", "fp", "fn", "tn", "sensitivity", "specificity"):
            print(f"{key}\t{rep[key]}")
        tm = {}
        if args.tm_table:
            with open(args.tm_table, encoding="utf-8-sig", newline="") as fh:
                tm = {r["protein_id"]: r["pred_hel"] for r in csv.DictReader(fh, delimiter="\t")}
        for pid in rep["nonmember_hits"]:
            print(f"nonmember_hit\t{pid}\tn_tm={tm.get(pid, 'NA')}")
        for pid in rep["missed_members"]:
            print(f"missed_member\t{pid}")
        return None
    if args.cmd == "allergen-lso":
        ids = [p.id for p in read_fasta(args.allergen_fasta)]
        rep = allergen.lso_report(args.blast, ids)
        print(f"sequences\t{rep['n_sequences']}\tspecies\t{rep['n_species']}")
        for r in rep["recall"]:
            print(f"{r['rule']}\t{r['recovered']}/{r['n']}")
        return None
    rows, summary = panel_check(args.calls_long, args.panel)
    for r in rows:
        print(
            "\t".join(
                [r["protein"], r["call"], r["variant"], r["expected"], r["observed"], r["verdict"]]
            )
        )
    print(json.dumps(summary), file=sys.stderr)
    return None


def main(argv=None):
    try:
        run(build_parser().parse_args(argv))
    except (ValueError, KeyError, TaxonError, OSError) as err:
        print(f"cellsurface_sorting_hat_calibrate: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

In `pyproject.toml`, `[project.scripts]`, add:

```toml
cellsurface_sorting_hat_module = "cellsurface_sorting_hat.modules.cli:main"
cellsurface_sorting_hat_calibrate = "cellsurface_sorting_hat.calibration.cli:main"
```

In `tests/surface_glyco/test_package.py`, add the same two entries to the `scripts` dictionary that `test_pyproject_names_match_the_decision` compares (Plan 1 added `cellsurface_sorting_hat` there).

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_module_cli.py tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py -q`
Expected: `28 passed`. Then `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/surface_glyco/test_package.py::test_pyproject_names_match_the_decision -q` Expected: `1 passed`.

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/cli.py src/cellsurface_sorting_hat/calibration/panel.py src/cellsurface_sorting_hat/calibration/cli.py pyproject.toml tests/surface_glyco/test_package.py tests/cellsurface_sorting_hat/modules/test_module_cli.py tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py
git commit -m "feat(sorting-hat): module and calibration commands

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 10: Job scripts

**Files:**
- Create: `scripts/sorting_hat/signalp_gpu.sbatch`, `pfam_hmmsearch.sbatch`, `repeats.sbatch`, `blast_allergen.sbatch`, `tmhmm.sbatch`, `fetch_proteomes.sh`, `submit_modules.sh`
- Test: `tests/cellsurface_sorting_hat/test_data_and_scripts.py`

**Interfaces:**
- Consumes: `load_family_table` (Task 2) and the data files of Tasks 2 and 7.
- Produces: raw tool outputs under `$WORKDIR/raw/{signalp,pfam,repeats,allergen,tmhmm}/`. `pfam_hmmsearch.sbatch` fetches by model name, stops if the model names or accessions in the database differ from the family table, and checks that `PFAM_RELEASE` is a whole directory name (`*Pfam38.2/`). `tmhmm.sbatch` runs in `$TMP` (TMHMM writes `TMHMM_<pid>` into the working directory).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_data_and_scripts.py` with exactly this content:

```python
"""The shipped family table, the Phase C species table and the job scripts."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

from cellsurface_sorting_hat.modules.pfam import load_family_table

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = sorted((ROOT / "scripts" / "sorting_hat").glob("*"))


def test_shipped_family_table_loads_and_starts_with_every_family_inactive():
    families = load_family_table(ROOT / "data" / "sorting_hat" / "family_table.tsv")
    assert len(families) == 15
    assert not any(
        f.active for f in families
    )  # a family is made active only after its specificity test
    assert {f.module for f in families} == {"pfam_adhesion", "pfam_allergen"}


def test_phasec_species_table_has_one_row_per_species():
    with open(ROOT / "data" / "sorting_hat" / "phasec_set_species.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len({r["scientific_name"] for r in rows}) == len(rows) == 6


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_scripts_are_valid_bash_and_follow_the_site_rules(script):
    text = script.read_text()
    assert subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode == 0
    assert "BASH_SOURCE" not in text  # breaks under sbatch
    if script.suffix == ".sbatch":
        assert "SCRATCH:?" in text and "#SBATCH" in text and "set -euo pipefail" in text
        assert "/tmp" not in text


def test_there_is_a_script_for_every_job_that_submit_modules_starts():
    names = {p.stem for p in SCRIPTS if p.suffix == ".sbatch"}
    text = (ROOT / "scripts" / "sorting_hat" / "submit_modules.sh").read_text()
    for job in text.split("for job in ", 1)[1].split(";", 1)[0].split():
        assert job in names


@pytest.mark.skipif(shutil.which("shellcheck") is None, reason="shellcheck is not installed")
@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_shellcheck(script):
    assert (
        subprocess.run(["shellcheck", "-S", "warning", str(script)], capture_output=True).returncode
        == 0
    )


def test_the_pfam_job_checks_model_names_against_the_family_table():
    text = (ROOT / "scripts" / "sorting_hat" / "pfam_hmmsearch.sbatch").read_text()
    assert "differ from the family table" in text and "hmmfetch -f" in text and "--cut_ga" in text


def test_the_pfam_job_fetches_by_model_name_and_the_table_names_are_unique():
    """hmmfetch -f on the pressed Pfam database needs the NAME; an accession without a version is refused."""
    text = (ROOT / "scripts" / "sorting_hat" / "pfam_hmmsearch.sbatch").read_text()
    assert "cut -f2 | sort -u" in text and "cut -f1 | sort -u" not in text
    families = load_family_table(ROOT / "data" / "sorting_hat" / "family_table.tsv")
    names = [f.name for f in families]
    assert len(names) == len(set(names))
    assert "Pfam${PFAM_RELEASE}/" in text  # the release must be a whole directory name


def test_every_sbatch_script_cleans_its_scratch_directory_on_exit():
    for script in SCRIPTS:
        if script.suffix == ".sbatch":
            assert "trap 'rm -rf \"$TMP\"' EXIT" in script.read_text(), script.name
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_data_and_scripts.py -q`
Expected: failures (the script glob finds no script, so the tests that read `submit_modules.sh` or `pfam_hmmsearch.sbatch` cannot open them). The family table and species table tests pass because Tasks 2 and 7 created those files.

- [ ] **Step 3: Write the scripts**

Make the files executable (`chmod +x scripts/sorting_hat/*`). Each job script takes its inputs from exported variables, uses `$SCRATCH`, removes its scratch directory with a `trap`, and copies results back with a temporary name.

Create `scripts/sorting_hat/signalp_gpu.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# SignalP 6 (GPU build, fast mode) on one proteome FASTA. Output: $WORKDIR/raw/signalp/prediction_results.txt
# Submit: sbatch --export=ALL,FASTA=/path/p.faa,WORKDIR=/path/work signalp_gpu.sbatch
# Measured rate on gpu12: about 192 proteins/s (analysis/step1_compare/jobs/j1_features.sh, job 29280458).
# The GPU build refuses to run without a CUDA device. Do not use the CPU build (about 2 s/protein).
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH -J csh_signalp
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_signalp.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$WORKDIR/raw/signalp"
module load signalp/6-gpu
command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH" >&2; exit 1; }
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
signalp6 --fastafile "$TMP/in.fasta" --organism eukarya --output_dir "$TMP/out" \
  --format none --mode fast --write_procs 4 --torch_num_threads 8
cp "$TMP/out/prediction_results.txt" "$WORKDIR/raw/signalp/.tmp.prediction_results.txt"
mv "$WORKDIR/raw/signalp/.tmp.prediction_results.txt" "$WORKDIR/raw/signalp/prediction_results.txt"
signalp6 --version > "$WORKDIR/raw/signalp/version.txt" 2>&1 || true
```

Create `scripts/sorting_hat/pfam_hmmsearch.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# Search the Pfam models of the family table with hmmsearch --cut_ga (as analysis/cys_candidates/01_known_family_hmm.sh).
# Output: $WORKDIR/raw/pfam/domtbl.txt and provenance.json (path, resolved path, sha256 of Pfam-A.hmm).
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=...,FAMILY_TABLE=...,PFAM_RELEASE=38.2 pfam_hmmsearch.sbatch
# PFAM_RELEASE is stated by the person who submits; it is checked against the file name of the resolved path.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J csh_pfam
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
: "${FAMILY_TABLE:?export FAMILY_TABLE (data/sorting_hat/family_table.tsv) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_pfam.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/pfam"
mkdir -p "$TMP" "$OUT"
module load hmmer/3.4
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"Pfam${PFAM_RELEASE}/"*) ;; *) echo "FATAL: $RESOLVED is not in a directory named *Pfam${PFAM_RELEASE}" >&2; exit 1 ;; esac
# hmmfetch -f looks up by the model NAME (column 2) in the pressed database, not by the accession
tail -n +2 "$FAMILY_TABLE" | cut -f2 | sort -u > "$TMP/keys.txt"
hmmfetch -f "$PFAM_HMM" "$TMP/keys.txt" > "$TMP/families.hmm"
N=$(grep -c '^NAME' "$TMP/families.hmm" || true)
[ "$N" -eq "$(wc -l < "$TMP/keys.txt")" ] || { echo "FATAL: hmmfetch found $N of $(wc -l < "$TMP/keys.txt") models" >&2; exit 1; }
awk '/^NAME/{n=$2} /^ACC/{split($2,a,"."); print a[1]"\t"n}' "$TMP/families.hmm" | sort > "$TMP/fetched.txt"
tail -n +2 "$FAMILY_TABLE" | cut -f1,2 | sort > "$TMP/table.txt"
diff -q "$TMP/fetched.txt" "$TMP/table.txt" > /dev/null || { echo "FATAL: model names in the database differ from the family table" >&2; diff "$TMP/fetched.txt" "$TMP/table.txt" >&2; exit 1; }
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
hmmsearch --cut_ga --cpu "${SLURM_CPUS_ON_NODE:-2}" --noali -o /dev/null --domtblout "$TMP/domtbl.txt" "$TMP/families.hmm" "$TMP/in.fasta"
cp "$TMP/domtbl.txt" "$OUT/.tmp.domtbl.txt" && mv "$OUT/.tmp.domtbl.txt" "$OUT/domtbl.txt"
printf '{"pfam_hmm": "%s", "resolved": "%s", "release": "%s", "sha256": "%s", "hmmer": "%s"}\n' \
  "$PFAM_HMM" "$RESOLVED" "$PFAM_RELEASE" "$(sha256sum "$RESOLVED" | cut -d' ' -f1)" "$(hmmsearch -h | sed -n 2p)" > "$OUT/provenance.json"
```

Create `scripts/sorting_hat/repeats.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# Run both repeat detectors on one proteome FASTA. Output: $WORKDIR/raw/repeats/{repeat02,repeat14}.tsv
# Submit: sbatch --export=ALL,PROJ_ROOT=...,FASTA=...,WORKDIR=... repeats.sbatch
# Run time is not measured for a whole proteome; record it from the job (sacct) in the run report.
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J csh_repeats
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_repeats.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/repeats"
mkdir -p "$TMP" "$OUT"
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
/usr/bin/python3.12 "$PROJ_ROOT/analysis/cocci_repeats/02_repeat_profile.py" "$TMP/in.fasta" --out "$TMP/repeat02.tsv"
/usr/bin/python3.12 "$PROJ_ROOT/analysis/cocci_repeats/14_repeat_detect_general.py" "$TMP/in.fasta" --out "$TMP/repeat14.tsv"
for m in repeat02 repeat14; do cp "$TMP/$m.tsv" "$OUT/.tmp.$m.tsv" && mv "$OUT/.tmp.$m.tsv" "$OUT/$m.tsv"; done
```

Create `scripts/sorting_hat/blast_allergen.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# BLASTP of one proteome against the WHO/IUIS fungal allergen sequences.
# Output: $WORKDIR/raw/allergen/blast.tsv (columns: qseqid sseqid pident length qlen slen bitscore evalue)
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=...,ALLERGEN_FASTA=... blast_allergen.sbatch
# -max_target_seqs 200 is larger than the database (116 sequences), so no hit is lost to the
# early-stop behaviour of -max_target_seqs.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH -J csh_blast
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
: "${ALLERGEN_FASTA:?export ALLERGEN_FASTA (made by build_allergen_fasta) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_blast.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/allergen"
mkdir -p "$TMP" "$OUT"
module load ncbi-blast/2.14.0+
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
makeblastdb -in "$ALLERGEN_FASTA" -dbtype prot -out "$TMP/allergens" -logfile /dev/null
blastp -query "$TMP/in.fasta" -db "$TMP/allergens" -evalue 1 -max_target_seqs 200 \
  -outfmt "6 qseqid sseqid pident length qlen slen bitscore evalue" -num_threads "${SLURM_CPUS_ON_NODE:-2}" \
  -out "$TMP/blast.tsv"
cp "$TMP/blast.tsv" "$OUT/.tmp.blast.tsv" && mv "$OUT/.tmp.blast.tsv" "$OUT/blast.tsv"
blastp -version | head -1 > "$OUT/version.txt"
```

Create `scripts/sorting_hat/tmhmm.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# TMHMM 2.0c on one proteome FASTA. Output: $WORKDIR/raw/tmhmm/tmhmm.tsv (columns as analysis/cocci_spherule/01b_tmhmm.sh)
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=... tmhmm.sbatch
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH -J csh_tmhmm
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_tmhmm.$$"
trap 'rm -rf "$TMP"' EXIT
OUT="$WORKDIR/raw/tmhmm"
mkdir -p "$TMP" "$OUT"
module load tmhmm/2.0c
cd "$TMP"  # TMHMM writes a TMHMM_<pid> directory into the current directory
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
tmhmm -short "$TMP/in.fasta" > "$TMP/tmhmm.out"
awk 'BEGIN{OFS="\t"; print "protein_id","len","exp_aa","first60","pred_hel","topology"}
  {for(i=2;i<=NF;i++){sub(/^[^=]*=/,"",$i)} print $1,$2,$3,$4,$5,$6}' "$TMP/tmhmm.out" > "$TMP/tmhmm.tsv"
cp "$TMP/tmhmm.tsv" "$OUT/.tmp.tmhmm.tsv" && mv "$OUT/.tmp.tmhmm.tsv" "$OUT/tmhmm.tsv"
```

Create `scripts/sorting_hat/fetch_proteomes.sh` with exactly this content:

```bash
#!/bin/bash
# Download the two extra A. fumigatus proteomes (decision D13) and write provenance records.
# Run on the login node: WORKDIR=/path/work bash fetch_proteomes.sh
# A1163 (CEA10, FGSC A1163, CBS 144.89): UniProt proteome UP000001699 (9,942 proteins on 2026-10-04).
# W72310: NCBI GCA_040167795.1 (UCR_Afum_W72310_1.0).
set -euo pipefail
: "${WORKDIR:?export WORKDIR}"
DEST="$WORKDIR/proteomes"
mkdir -p "$DEST"
curl -fsSL "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001699&format=fasta" -o "$DEST/.tmp.A1163.faa"
mv "$DEST/.tmp.A1163.faa" "$DEST/Afum_A1163_UP000001699.faa"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/167/795/GCA_040167795.1_UCR_Afum_W72310_1.0"
curl -fsSL "$BASE/GCA_040167795.1_UCR_Afum_W72310_1.0_protein.faa.gz" -o "$DEST/.tmp.W72310.faa.gz"
mv "$DEST/.tmp.W72310.faa.gz" "$DEST/Afum_W72310_GCA_040167795.1.faa.gz"
echo "downloaded; now write provenance (see Task 12 of the plan)"
```

Create `scripts/sorting_hat/submit_modules.sh` with exactly this content:

```bash
#!/bin/bash
# Submit the module jobs for one proteome. Usage:
#   PROJ_ROOT=... FASTA=... WORKDIR=... FAMILY_TABLE=... PFAM_RELEASE=38.2 ALLERGEN_FASTA=... bash submit_modules.sh
# Prints one job ID per line. Wait for the jobs (squeue / sacct), then run the converters
# (module commands of Task 13) and the core command. All paths come from PROJ_ROOT.
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
: "${FASTA:?export FASTA}"
: "${WORKDIR:?export WORKDIR}"
: "${FAMILY_TABLE:?export FAMILY_TABLE}"
: "${PFAM_RELEASE:?export PFAM_RELEASE}"
: "${ALLERGEN_FASTA:?export ALLERGEN_FASTA}"
S="$PROJ_ROOT/scripts/sorting_hat"
mkdir -p "$WORKDIR/logs"
for job in signalp_gpu pfam_hmmsearch repeats blast_allergen tmhmm; do
  sbatch --parsable --export=ALL -o "$WORKDIR/logs/$job.%j.log" -e "$WORKDIR/logs/$job.%j.log" "$S/$job.sbatch"
done
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_data_and_scripts.py -q`
Expected: `13 passed, 7 skipped` (the skips are `shellcheck` tests; run `shellcheck -S warning scripts/sorting_hat/*` on a machine that has it). Then the whole suite: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q`. Expected: `256 passed, 7 skipped`.

Check the key extraction of the Pfam job against the real database (login node; this is the step that failed when keys were accessions):

```bash
module load hmmer/3.4
T="$SCRATCH/pfamcheck"; mkdir -p "$T"
tail -n +2 data/sorting_hat/family_table.tsv | cut -f2 | sort -u > "$T/keys.txt"
hmmfetch -f /bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm "$T/keys.txt" | grep -c '^NAME'   # expected: 15
readlink -f /bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm                                  # expected: a path with Pfam38.2/
```

- [ ] **Step 5: Mutation checks for Tasks 1 to 9 (confirm the tests can fail)**

Make each change in turn, run `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q`, confirm that the named tests fail, then undo the change with `git checkout -- <file>`.

| Change | File | Tests that fail (measured) |
|---|---|---|
| In `status_from_measure`, replace `if spec and enough and narrow:` by `if enough and narrow:` | `calibration/measure.py` | `test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate`, `test_a_set_with_no_specificity_cell_is_at_most_smoke_even_if_labelled_estimate` |
| In `status_from_measure`, delete the loop over `n_clusters_pos`, `n_clusters_neg` | `calibration/measure.py` | `test_a_tight_interval_with_few_recorded_clusters_is_not_an_estimate` |
| In `cluster_bootstrap`, return `min(lo, value)` and `max(hi, value)` without the Wilson bounds | `calibration/intervals.py` | `test_all_positives_called_gives_a_width_not_a_zero_width_interval`, `test_one_cluster_per_class_cannot_give_an_estimate` |
| In the `signalp` command, set `artefact_digest=_file_digest(args.results)` | `modules/cli.py` | `test_the_signalp_module_identity_does_not_depend_on_the_proteome` |
| In `pfam_rows`, delete the two lines of the `no_tm` condition | `modules/pfam.py` | `test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices`, `test_pfam_command_applies_the_no_tm_condition_from_the_tm_table` |
| In `tm_rows`, set `mature = len(starts)` | `modules/lookups.py` | `test_tm_rows_count_helices_after_the_signal_peptide_window`, `test_pfam_command_applies_the_no_tm_condition_from_the_tm_table` |
| In `make_entry`, delete the `if cap == SMOKE and status == ESTIMATED:` block | `calibration/measure.py` | `test_a_leakage_cap_limits_an_estimate_to_smoke`, `test_leakage_other_than_none_caps_an_estimate_at_smoke` |
| In `best_hit_other_species`, delete ` or species[q] == species[s]` | `modules/allergen.py` | `test_other_species_hits_exclude_the_same_species` |
| In `entries_from_phasec`, add `entry["status"] = "estimated"` after the `make_entry` line | `calibration/phasec.py` | `test_the_status_is_the_weaker_of_the_phase_c_label_and_the_rule`, `test_a_set_with_no_specificity_cell_is_at_most_smoke_even_if_labelled_estimate` |
| In `_merge_entries`, replace the condition `any(old_file.get(k) != now.get(k) for k in keys)` by `False` | `calibration/cli.py` | `test_a_changed_module_identity_drops_the_old_status_entries` |
| In `compare_runs`, set `other = pid` | `calibration/stability.py` | `test_only_identical_sequences_are_compared_by_sha256_not_by_id` |
| In `run`, replace `if not active:` by `if False:` | `modules/cli.py` | `test_pfam_command_writes_unavailable_modules_when_no_family_is_active`, `test_pfam_command_applies_the_no_tm_condition_from_the_tm_table` |

- [ ] **Step 6: Commit**

```bash
git branch --show-current   # must print sorting-hat-modules
git add scripts/sorting_hat tests/cellsurface_sorting_hat/test_data_and_scripts.py
git commit -m "feat(sorting-hat): SLURM job scripts for the module tools

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 11: Proteomes and taxonomy (HPCC)

**Files:**
- Create (not committed): `_workdir/sorting_hat/proteomes/*`, `_workdir/sorting_hat/taxdump/*`, `_workdir/sorting_hat/provenance/*.json`

**Interfaces:**
- Consumes: `write_provenance` (Task 8), `fetch_proteomes.sh` (Task 10).
- Produces: six proteome FASTA files with provenance records, and `nodes.dmp`/`names.dmp`.

Use one shell for Tasks 11 to 14. Set the environment first (the package is not installed in the system Python):

```bash
export PROJ_ROOT=$PWD WORKDIR=$PWD/_workdir/sorting_hat PYTHONPATH=$PWD/src
```

- [ ] **Step 1: Download the taxonomy and the two extra proteomes**

```bash
mkdir -p "$WORKDIR/taxdump" && cd "$WORKDIR/taxdump"
curl -fsSL https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -o taxdump.tar.gz && tar xzf taxdump.tar.gz nodes.dmp names.dmp
cd "$PROJ_ROOT" && bash scripts/sorting_hat/fetch_proteomes.sh
ls "$WORKDIR/proteomes"
```

Expected: `Afum_A1163_UP000001699.faa` and `Afum_W72310_GCA_040167795.1.faa.gz`.

- [ ] **Step 2: Check the taxa**

Status applies to a protein only if its taxon is a descendant of a tested taxon (species level). Check the strain taxa you will pass:

```bash
/usr/bin/python3.12 - <<'E'
import os
from cellsurface_sorting_hat.taxonomy import Lineage
from cellsurface_sorting_hat.calibration.phasec import read_names, species_taxid
W = os.environ["WORKDIR"] + "/taxdump"
lin, names = Lineage.from_nodes_dmp(W + "/nodes.dmp"), read_names(W + "/names.dmp")
afum = species_taxid(names, "Aspergillus fumigatus")
for label, taxon in {"Af293": 330879, "A1163": 451804, "W72310": 746128}.items():
    print(label, taxon, lin.is_descendant_or_self(taxon, afum), taxon == afum)
E
```

Expected: `True` in the third column for all three. The fourth column shows whether the taxon is the species ID itself (W72310 has no strain taxon: it is expected to be `True` there). Taxon 451804 for A1163 and 330879 for Af293 were read from the microbiology review of 2026-10-05; confirm them against the UniProt proteome record of UP000001699 (`taxonomy.taxonId`) and of UP000002530. If a taxon is not under the species, stop.

- [ ] **Step 3: Write provenance records and check the counts**

```bash
/usr/bin/python3.12 - <<'E'
import os, datetime
from cellsurface_sorting_hat.proteomes import write_provenance
W = os.environ["WORKDIR"]; P = os.environ["PROJ_ROOT"]
today = datetime.date.today().isoformat()
os.makedirs(W + "/provenance", exist_ok=True)
sets = [
 ("Cimm_RS", "NCBI RefSeq GCF_000149335.2", P + "/_workdir/cocci_spherule/ref/GCF_000149335.2_ASM14933v2_protein.faa.gz", 9910, 9910),
 ("Afum_Af293_UniProt", "UniProt UP000002530 (the Phase C gene models)", P + "/_workdir/step1_compare/downloads/UP000002530.fasta.gz", 9000, 10500),
 ("Afum_Af293_Fungi5k", "Fungi_5k input (a different annotation)", "/bigdata/stajichlab/shared/projects/Fungi_5k/input/Aspergillus_fumigatus_Af293.proteins.fa", 9161, 9161),
 ("Scer_S288C", "SGD orf_trans_all.fasta.gz (includes dubious ORFs)", P + "/_workdir/step1_compare/downloads/orf_trans_all.fasta.gz", 6000, 7500),
 ("Afum_A1163", "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001699", W + "/proteomes/Afum_A1163_UP000001699.faa", 9000, 10500),
 ("Afum_W72310", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/167/795/GCA_040167795.1_UCR_Afum_W72310_1.0/", W + "/proteomes/Afum_W72310_GCA_040167795.1.faa.gz", 9000, 11000),
]
for name, source, fasta, lo, hi in sets:
    r = write_provenance(f"{W}/provenance/{name}.json", name, source, fasta, today, lo, hi)
    print(name, r["n_proteins"], r["n_unique_sequences"], r["n_invalid"])
E
```

Expected: one line per proteome. RS prints 9,910 and the Fungi_5k Af293 file 9,161 (counted 2026-10-05). The ranges for the other proteomes are plausibility ranges chosen from known proteomes (UniProt *A. fumigatus* 9,942 for A1163; S288C 6,722 entries including 683 dubious ORFs and 14 with an internal `*`); they are not measured. Record the counts. A `ProteomeError` means the file is not the expected proteome: stop.

For S288C, make a second FASTA without the dubious ORFs (headers that contain `Dubious ORF`) and run both; report calls with and without them.

- [ ] **Step 4: Commit nothing** (the files are under `_workdir`, which is git-ignored). Copy the provenance JSON files to `docs/reports/data/sorting_hat/provenance/` for Task 15.

### Task 12: First end-to-end run on A. fumigatus Af293 (HPCC)

**Files:**
- Create (not committed): `_workdir/sorting_hat/Afum_Af293_UniProt/{raw,modules,status,out,logs}`, and the same for `Afum_Af293_Fungi5k`

**Interfaces:**
- Consumes: Tasks 1 to 11 and Plan 1 (`cellsurface_sorting_hat`).
- Produces: `out/calls.long.tsv.gz`, `out/report.md`, and the job resource table.

- [ ] **Step 1: Install the package in a virtual environment (Python 3.12)**

```bash
/usr/bin/python3.12 -m venv "$SCRATCH/csh-venv"
"$SCRATCH/csh-venv/bin/pip" install -e . --no-deps && "$SCRATCH/csh-venv/bin/pip" install pyyaml numpy
export PATH="$SCRATCH/csh-venv/bin:$PATH"
cellsurface_sorting_hat --help | head -2 && cellsurface_sorting_hat_module --help | head -2 && cellsurface_sorting_hat_calibrate --help | head -2
```

Expected: three usage lines. (`$SCRATCH` is set in a job shell; on a login node use `/scratch/$USER` if it is empty.)

- [ ] **Step 2: Build the allergen database FASTA and submit the module jobs for the UniProt Af293 proteome**

```bash
export FASTA=$PROJ_ROOT/_workdir/step1_compare/downloads/UP000002530.fasta.gz
export PROT=Afum_Af293_UniProt WORKDIR=$PROJ_ROOT/_workdir/sorting_hat/Afum_Af293_UniProt
export FAMILY_TABLE=$PROJ_ROOT/data/sorting_hat/family_table.tsv PFAM_RELEASE=38.2
export ALLERGEN_FASTA=$PROJ_ROOT/_workdir/sorting_hat/iuis_fungal_allergens.faa
/usr/bin/python3.12 -c "
import os
from cellsurface_sorting_hat.modules.allergen import build_allergen_fasta as b
n, skipped = b('analysis/allergen_scoping/iuis_fungal_isoallergens.tsv', os.environ['ALLERGEN_FASTA'], 'analysis/allergen_scoping/iuis_fungal_allergens.tsv')
print(n, skipped)"
grep -c '>' "$ALLERGEN_FASTA"
bash scripts/sorting_hat/submit_modules.sh | tee "$WORKDIR/jobs.txt"
```

Expected: the builder prints 111 and lists 5 skipped entries (measured 2026-10-05 on the extract of 2026-10-04: Epi p 1.0101 is free text; Asp fl 13.0101, Cur l 1.0101, Tri t 1.0101 and Tri t 4.0101 are peptide fragments of 9 to 29 residues, so two of the four Onygenales entries are not searchable); 97 of the 111 have IUIS evidence text in the meta table; `grep -c` gives the same number as the builder; five job IDs print. Wait until `squeue -u $USER` shows none of them. If a job fails, read its log in `$WORKDIR/logs/` before resubmitting; do not resubmit in a loop.

- [ ] **Step 3: Convert the tool outputs to module tables and run the core command**

```bash
T=$PROJ_ROOT/_workdir/sorting_hat/taxdump
M="cellsurface_sorting_hat_module"
FQ=$PROJ_ROOT/_workdir/sorting_hat/Afum_Af293_UniProt.faa
zcat "$FASTA" > "$FQ"        # the module commands read plain FASTA or .gz; a plain copy keeps IDs identical for the tools
PFAMJ=$WORKDIR/raw/pfam/provenance.json
SHA=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['sha256'])")
HMMER=$(/usr/bin/python3.12 -c "import json;print(json.load(open('$PFAMJ'))['hmmer'].split()[2])")
$M signalp --fasta $FQ --workdir $WORKDIR --results $WORKDIR/raw/signalp/prediction_results.txt --signalp-version "$(head -1 $WORKDIR/raw/signalp/version.txt)"
$M tm --fasta $FQ --workdir $WORKDIR --table $WORKDIR/raw/tmhmm/tmhmm.tsv
$M pfam --fasta $FQ --workdir $WORKDIR --domtbl $WORKDIR/raw/pfam/domtbl.txt --family-table $FAMILY_TABLE --pfam-release 38.2 \
  --pfam-sha256 "$SHA" --hmmer-version "$HMMER" --sp-module step1_rule@R0 --tm-module tm
$M repeat02 --fasta $FQ --workdir $WORKDIR --table $WORKDIR/raw/repeats/repeat02.tsv --script analysis/cocci_repeats/02_repeat_profile.py
$M repeat14 --fasta $FQ --workdir $WORKDIR --table $WORKDIR/raw/repeats/repeat14.tsv --script analysis/cocci_repeats/14_repeat_detect_general.py
$M allergen --fasta $FQ --workdir $WORKDIR --blast $WORKDIR/raw/allergen/blast.tsv --allergen-fasta $ALLERGEN_FASTA --blast-version "$(head -1 $WORKDIR/raw/allergen/version.txt)"
cellsurface_sorting_hat --fasta $FQ --taxon 330879 --taxdump $T/nodes.dmp --workdir $WORKDIR --out $WORKDIR/out
head -60 $WORKDIR/out/report.md
```

Expected: exit codes 0. The report lists `repeat02`, `repeat14`, `step1_rule@R0`, `tm` and `allergen_homology` as `ok`; `pfam_adhesion` and `pfam_allergen` as `unavailable` (no family is active yet) and `wall_family_domain` and `iuis_allergen_homolog` consequently unknown for the domain part; `antigen_lookup` is listed under "Modules the rules need and that were not found" (it is run for RS only). Every status is `unvalidated` until Task 13.

If `pfam` stops with `no '# [ok]' trailer` or `not run with --cut_ga`, the domain table is not from a finished `--cut_ga` run: resubmit the Pfam job.

- [ ] **Step 4: Compare with the numbers already measured, and record the job resources**

```bash
zcat $WORKDIR/out/calls.long.tsv.gz | awk -F'\t' '$2=="iuis_allergen_similarity" && $4=="called"' | wc -l
zcat $WORKDIR/out/calls.long.tsv.gz | awk -F'\t' '$2=="iuis_allergen_homolog" && $4=="called"' | wc -l
sacct -j "$(paste -sd, $WORKDIR/jobs.txt)" --format=JobID,JobName%20,Partition,Elapsed,MaxRSS,AllocCPUS,State
```

The numbers measured on the **Fungi_5k** Af293 file (9,161 proteins) with ncbi-blast 2.14.0+ were 105 proteins at identity 35% and aligned length 80 (the `iuis_allergen_similarity` field) and 41 at 70% identity and 80% coverage; 30 reached 95% identity. For the UniProt proteome (the one used here) no number was measured: record yours and compare it with the Fungi_5k run of Step 5. Of the 30 proteins at 95% or more, most are the *A. fumigatus* allergens that are in the reference set: report the counts split into "hit to an allergen of the same species" and "hit to an allergen of another species" (use `allergen_species`). Write the `sacct` table into the report (Task 15); these are the first measured times for the whole-proteome Pfam, repeat and TMHMM jobs.

- [ ] **Step 5: Run the Fungi_5k Af293 file as a second annotation**

Repeat Steps 2 to 4 with `FASTA=/bigdata/stajichlab/shared/projects/Fungi_5k/input/Aspergillus_fumigatus_Af293.proteins.fa`, `PROT=Afum_Af293_Fungi5k` and `WORKDIR=$PROJ_ROOT/_workdir/sorting_hat/Afum_Af293_Fungi5k` (the allergen FASTA is reused; do not run `phasec` calibration for this run: it is not the annotation that was measured). Compare R0 calls of the two annotations for the genes they share:

```bash
/usr/bin/python3.12 - <<'E'
import os
from cellsurface_sorting_hat.fasta import read_fasta
P = os.environ["PROJ_ROOT"] + "/_workdir"
a = read_fasta(P + "/step1_compare/downloads/UP000002530.fasta.gz")
b = read_fasta("/bigdata/stajichlab/shared/projects/Fungi_5k/input/Aspergillus_fumigatus_Af293.proteins.fa")
sha_a = {p.sha256 for p in a}
print("UniProt", len(a), "Fungi_5k", len(b), "identical sequences", sum(p.sha256 in sha_a for p in b))
E
```

Record: the identical count (4,743 measured 2026-10-05 by the microbiology review), and, using the Task 14 stability command with a best-reciprocal-hit table, the R0 agreement for proteins that share their C-terminus but differ at the start. The bioinformatics review measured a change of the R0 call in 5 of 80 such pairs (6.3%, Wilson 95% interval about 2.7% to 13.8%; small sample).

- [ ] **Step 6: No commit** (outputs are under `_workdir`). Keep `jobs.txt` and the `sacct` output for the report.

### Task 13: Calibrations that can be measured today (HPCC)

**Files:**
- Create (not committed): `$WORKDIR/status/*.json`, `_workdir/sorting_hat/calibration/*`
- Modify: `data/sorting_hat/family_table.tsv` (only for a family that passed review and has a sign-off)

**Interfaces:**
- Consumes: Tasks 6, 7, 9 and the Phase C `metrics.json`.
- Produces: the status source of `step1_rule@R0` (per species), the allergen leave-species-out table, one specificity review table per Pfam family.

- [ ] **Step 1: Status source for rule R0 from the Phase C metrics (UniProt Af293 run)**

```bash
cellsurface_sorting_hat_calibrate phasec --workdir $WORKDIR --metrics _workdir/step1_compare/phasec/metrics.json \
  --set-species data/sorting_hat/phasec_set_species.tsv --names-dmp $T/names.dmp
cellsurface_sorting_hat --fasta $FQ --taxon 330879 --taxdump $T/nodes.dmp --workdir $WORKDIR --out $WORKDIR/out
grep -A12 "Module calibration" $WORKDIR/out/report.md
```

Expected: `step1_rule@R0` for taxon 330879 reports **status `smoke`** (the *A. fumigatus* entry: 19 positives, 45 negatives, Phase C label `smoke test`), calibration set `S3-Eurotiomycetes:Afum_ASPFU`, sensitivity 0.947 [0.833, 1.000] and specificity 1.000 widened with Wilson on 45 negatives (about [0.92, 1.000]); numbers from `metrics.json`, read 2026-10-05. The `signal_peptide_protein[R0]` rows carry status `smoke` with basis `step1_rule@R0:taxon:<A. fumigatus species ID>`. For *A. nidulans* (the one species with Phase C label `estimate`) the status is `estimated`. The pooled Eurotiomycetes value that earlier drafts used (0.727) is mostly *A. nidulans* and is not written. The R0 entries are an UniProt-annotation measurement: the Fungi_5k run (Task 12 Step 5) must not use them.

- [ ] **Step 2: Allergen leave-species-out recall (sensitivity only)**

```bash
module load ncbi-blast/2.14.0+
mkdir -p _workdir/sorting_hat/calibration && cd _workdir/sorting_hat/calibration
makeblastdb -in $ALLERGEN_FASTA -dbtype prot -out alg -logfile /dev/null
blastp -query $ALLERGEN_FASTA -db alg -evalue 1 -max_target_seqs 200 -outfmt "6 qseqid sseqid pident length qlen slen bitscore evalue" -out iuis_vs_iuis.tsv
cd $PROJ_ROOT && cellsurface_sorting_hat_calibrate allergen-lso --blast _workdir/sorting_hat/calibration/iuis_vs_iuis.tsv --allergen-fasta $ALLERGEN_FASTA
```

Expected: the first line `sequences <n> species <k>` with `<n>` equal to the FASTA count of Task 12; then one line per rule with `recovered/<n>`. The command refuses if any sequence has no BLAST line (BLAST masked or dropped it): find the cause, do not remove the check. Each number says how many known allergens a rule finds when the only reference allergens are from other species. It is not a specificity and writes no status entry. The `species` code is the genus and species letters of the allergen name (for example `Asp_f`); isoallergens and variants of one species are therefore never counted as relatives.

- [ ] **Step 3: Pfam family specificity review (one table per family and proteome)**

For each family in `data/sorting_hat/family_table.tsv`, run `pfam_hmmsearch.sbatch` (with the family table) on the proteomes of four clades (*A. fumigatus* Af293 UniProt, *S. cerevisiae* S288C, *C. immitis* RS, *C. neoformans* H99; their FASTA paths are in `analysis/step1_compare/species.tsv` and Task 11) and run TMHMM (`tmhmm.sbatch`) on the same files. Make `members.tsv` (TAB separated, header `pfam_acc`, `protein_id`) from the curated tables in `data/curated/` and the literature: proteins known to belong to the family by function, not because Pfam found a domain in them. Then, for each family and proteome:

```bash
cellsurface_sorting_hat_calibrate pfam-specificity --family PF05730 --domtbl <proteome>/raw/pfam/domtbl.txt \
  --members members.tsv --universe-fasta <proteome FASTA> --tm-table <proteome>/raw/tmhmm/tmhmm.tsv > <proteome>.PF05730.specificity.tsv
```

The numbers `sensitivity` and `specificity` in that output treat every protein that is not a known member as a negative (an unlabelled set). Do not report them as the specificity of the family. Read every `nonmember_hit` line (annotation, length, domain architecture, TMHMM helices) and classify it: false domain hit, uncharacterised true member, or receptor-like (CFEM with helices). Report the three counts per family, and for CFEM and hydrophobins the share of hits that carry the conserved Cys pattern (8 Cys). Write one line per hit class in the `specificity_note` column of the family table.

A family becomes `active` only when the owner signs off (decision C4): set `active = yes`, `active_by`, `active_date` in the family table and commit that change alone, with the review tables as evidence. `second_condition` must be set where the review shows a need (PA14 needs `signal_peptide`).

- [ ] **Step 4: Truth tables for the other modules (blocked on owner decisions C1 to C3)**

When an owner decision provides a truth table `truth.tsv` (`id`, `label` 0 or 1, `cluster`), the calibration command is:

```bash
cellsurface_sorting_hat_calibrate truth --workdir $WORKDIR --module <module> --calls-long <out>/calls.long.tsv.gz \
  --call <call name> --truth truth.tsv --calibration-set <name> --taxa <species-level taxon IDs> --leakage <none|partial|tuned_on_truth|in_reference|unknown> --notes "<what the set is>"
```

The run that produced `calls.long.tsv.gz` must contain the truth proteins. Use `--leakage tuned_on_truth` for `cocci_specificity_rank_top15` against the four anchors, and `in_reference` for an allergen that is in the IUIS set. Do not write an entry for a module whose truth set has no negatives. The command records the sensitivity bound with not-assessable positives counted as missed; report both.

- [ ] **Step 5: Commit the family table change and the review tables only**

Copy the `*.specificity.tsv` files to `docs/reports/data/sorting_hat/specificity/` first.

```bash
git branch --show-current   # must print sorting-hat-modules
git add data/sorting_hat/family_table.tsv docs/reports/data/sorting_hat/specificity
git commit -m "analysis(sorting-hat): Pfam specificity review tables; family table sign-offs

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 14: C. immitis RS and the stability of A. fumigatus strains (HPCC)

**Files:**
- Create (not committed): `_workdir/sorting_hat/{Cimm_RS,Afum_A1163,Afum_W72310,Scer_S288C}/...`

**Interfaces:**
- Consumes: Tasks 12 and 13; `compare_runs` (Task 7); `panel` command (Task 9).
- Produces: runs for the other proteomes, a strain comparison, and the report-only panel check.

- [ ] **Step 1: Run each remaining proteome the way Task 12 does**

Set `PROT`, `WORKDIR`, `FASTA` and the taxon: *C. immitis* RS (RefSeq protein file, taxon 246410), *S. cerevisiae* S288C (`orf_trans_all.fasta.gz`, taxon 559292; also the file without dubious ORFs), A1163 (taxon read in Task 11), W72310 (taxon 746128). Only the Af293 UniProt run got the Phase C status source; the other *A. fumigatus* runs get `calibrate phasec` too, because their taxa lie under the species, but the report must say that their gene models are not the measured annotation. For RS also build the lookup modules. The Cys-rich pipeline of `analysis/cys_candidates/README.md` must first run on the RefSeq FASTA (its IDs must be the `XP_...` IDs; a table made on the FungiDB IDs matches nothing and the command refuses). Use the per-proteome table that the pipeline writes (all tiers), not `candidates.tsv.gz`, which lists only the 304 non-`other` rows:

```bash
$M antigen --fasta $FASTA --workdir $WORKDIR --taxon 246410 --ranking analysis/cocci_antigens/cocci_antigen_ranking.tsv --protein-map _workdir/cocci_spherule/ref/protein_map.tsv
$M expression --fasta $FASTA --workdir $WORKDIR --taxon 246410 --table analysis/cocci_spherule/spherule_surface_table.tsv.gz --protein-map _workdir/cocci_spherule/ref/protein_map.tsv
$M cys --fasta $FASTA --workdir $WORKDIR --taxon 246410 --candidates <the all-tier table of the RefSeq run, for example _workdir/cys_candidates_refseq/<proteome>.tsv.gz>
```

Record, for RS: the number of proteins with `antigen_lookup` state `ok`, `not_in_reference` and `not_applicable` (the match rate of the ID map; 131 genes have more than one protein and 126 RefSeq proteins inherit the values of another transcript); how many of the 14 Tier 1 and 45 Tier 2 candidates have `cocci_specificity_rank_top15` called; how many of the calls have a signal peptide (expected 143 of 1,371, measured by the microbiology review); and how many of the 9,910 RS proteins have a signal peptide under rule R0 (460, counted from `analysis/cocci_repeats/signalp/CimmitisRS_FungiDB/prediction_results.txt` on 2026-10-05; the sequences of the FungiDB and RefSeq files are identical for all 9,910 proteins per the bioinformatics review, so the count should repeat). Print the expected number of false R0 calls per proteome (one minus specificity times the number of non-surface proteins), at stated assumed prevalences (for example 3% and 10%) and mark them assumed.

- [ ] **Step 2: Compare the three *A. fumigatus* runs**

```bash
/usr/bin/python3.12 - <<'E'
import os
from cellsurface_sorting_hat.calibration.stability import compare_runs
W = os.environ["PROJ_ROOT"] + "/_workdir/sorting_hat"
for other in ("Afum_A1163", "Afum_W72310", "Afum_Af293_Fungi5k"):
    r = compare_runs(f"{W}/Afum_Af293_UniProt/out", f"{W}/{other}/out")
    print(other, r["n_a"], r["n_b"], r["n_identical_sequences"], r["n_unique_identical_sequences"])
    for key, c in sorted(r["per_call"].items()):
        print("  ", key, dict(c))
E
```

Expected: the number of identical sequences and, per call, how many agree. This compares only byte-identical sequences, so it tests determinism of the software; it cannot show real differences between strains. Proteins without an identical sequence in the other proteome are not compared: report how many. Also make a best-reciprocal-hit table (`diamond blastp` both ways, or `blastp`, identity and coverage binned) and report, for the non-identical pairs, the agreement of `signal_peptide_protein[R0]` with 2x2 counts and Cohen's kappa. Name the cases CspA (repeat region varies between isolates), RodA, CalA and the Asp f allergens, and classify every difference as sequence difference, N-terminal gene model difference or presence/absence. Strain differences are partly annotation differences: UniProt, RefSeq and GenBank pipelines differ.

- [ ] **Step 3: Report-only panel check, one call per proteome**

Make `docs/reports/data/sorting_hat/panel.tsv` (columns `protein`, `call`, `variant`, `expected`, `source`, `tuning`) from the panel of spec section 4, with the renamed calls: SOWgp (`tandem_repeat_protein`, `cocci_specificity_rank_top15`; tuning `tuned` for the repeat detectors, which were set on it); Als1 and Flo11 (`tandem_repeat_protein`); Ag2/PRA and Rbt5 (`wall_family_domain`); RodA and CalA (`wall_family_domain`); Gel1 (`signal_peptide_protein` only); Asp f 1, Asp f 2 and Asp f 34 (`iuis_allergen_homolog`; tuning `in_reference`: they are in the reference set, so the expected result is guaranteed); at least one cross-species allergen row (a non-*Aspergillus* fungal protein with published IgE reactivity whose species is absent from the reference set); Hsp60 (`known_miss`). Each row needs the protein ID of the proteome it is in and a source (PMID or table row). A protein that is not in the proteome that was run gets `not_in_run`. Run the command once per proteome:

```bash
cellsurface_sorting_hat_calibrate panel --calls-long <out>/calls.long.tsv.gz --panel docs/reports/data/sorting_hat/panel.tsv
```

Expected: one line per panel row and a summary on stderr (`agree`, `disagree`, `not_in_run`, `known_miss`, `excluded_leakage`). Nothing is gated on the result; list every disagreement in the report with the reason.

- [ ] **Step 4: No commit** (outputs under `_workdir`); commit `panel.tsv` with the report.

### Task 15: Calibration and run report

**Files:**
- Create: `docs/reports/2026-10-05-sorting-hat-run-and-calibration.md`, `docs/reports/data/sorting_hat/provenance/*.json`, `docs/reports/data/sorting_hat/panel.tsv`

**Interfaces:**
- Consumes: the numbers from Tasks 11 to 14.
- Produces: the report that Plan 3 (acceptance) and the owner read.

- [ ] **Step 1: Write the report**

Start the report with this header, verbatim: "Research use only. This is not a regulatory allergenicity assessment and not a diagnostic result. No row is supported by an IgE, antibody or T-cell measurement. `signal_peptide_protein` means that SignalP calls a signal peptide and nothing more. Peptide-level calls only; glycan epitopes are not assessed. WHO/IUIS lists no *Coccidioides* allergen; delayed-type hypersensitivity skin-test reactivity is T-cell mediated and is not what the allergen columns address. The antigen ranking covers the *C. immitis* RS reference only (not *C. posadasii*). Non-protein adhesins such as galactosaminogalactan are not detected." Add the WHO/IUIS citation (allergen.org and a recent IUIS publication) and a licence line for every redistributed derived table.

Sections, each filled from a command above (copy the numbers, do not retype them from memory):
1. Proteomes: the provenance records (source, date, sha256, counts, annotation source).
2. Jobs: the `sacct` table (partition, elapsed, MaxRSS) for each module job.
3. Calibration table: one row per module, species and calibration set with positives, negatives, independent clusters where known, Sn and Sp with 95% intervals, per-stratum specificity, status, leakage, and "not measured" where nothing was measured. State the expected number of false calls per proteome at stated assumed prevalences.
4. Allergen: leave-species-out recall per rule (Task 13 step 2), and the Af293 counts split into same-species and other-species hits, with the exposure route and the evidence of the matched allergens. State that Af293 counts are self-recognition and are not module performance.
5. Pfam reviews: the three counts per family, the Cys-pattern share, and the sign-off decision per family.
6. Antigen lookup on RS: match rate, the share of calls with a signal peptide, Tier 1 and Tier 2 coverage; a statement that the cut of 15 was set after the anchors were seen, that Ag2/PRA, PRA2 and PRA3 are one family, and that 15% of the proteome is called. Add the sentence of `data/curated/antigens/README.md` that sequence homology does not transfer epitopes.
7. Strain stability (Task 14 step 2) with the best-reciprocal-hit comparison.
8. Panel check (disagreements and reasons).
9. What is still not measured, by module, and which owner decision (C1 to C4) blocks it.

- [ ] **Step 2: Commit**

```bash
git branch --show-current   # must print sorting-hat-modules
git add docs/reports/2026-10-05-sorting-hat-run-and-calibration.md docs/reports/data/sorting_hat
git commit -m "docs: sorting hat run and calibration report

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

## Self-review against the spec and Plan 1

| Item | Where |
|---|---|
| Module tables and run records per the Plan 1 contract; identity of tool, database, parameters | Tasks 1 to 5, 9 |
| `step1_rule@R0`, `pfam_adhesion`, `pfam_allergen`, `repeat02`, `repeat14`, `allergen_homology`, `antigen_lookup`, `cys_rich`, `expression`, `tm` | Tasks 1 to 5, 9 |
| Names that say what the tools measure; report limits; research-use header | Task 0, Task 15 |
| Family table with `active` sign-off and `signal_peptide`, `no_tm` conditions; unavailable while inactive | Tasks 2, 9, 13 |
| Allergen evidence fields (exposure, evidence text, species); leave-species-out recall; clean FASTA | Task 4, Task 13 |
| Antigen lookup: exact-taxon applicability, protein map, `idmap_method`, refusal on zero matches | Tasks 5, 9, 14 |
| Status rule: negatives, clusters, leakage cap, per-species entries, Phase C label and rule combined, identity check on merge | Tasks 6, 7, 9, 13 |
| `measure` shown in the report ("Module calibration" of Plan 1) | Tasks 6, 7, 13 |
| Strain stability by sequence hash; reciprocal best hits | Task 7 (identical sequences), Task 14 (reciprocal best hits, manual) |
| Proteome downloads with provenance; D13 proteomes; UniProt Af293 as the calibrated annotation | Tasks 8, 10, 11, 12 |
| SLURM job scripts (`$SCRATCH`, trap, atomic copy-back, name-based `hmmfetch`) | Tasks 10, 12 |
| Run-level check on Af293 on HPCC | Task 12 |
| Spec 3.8 GPU out-of-memory retry, job timeout and preemption handling | **Not implemented.** Fixed resources; the person resubmits after reading the log. Add to Plan 3 if jobs fail in practice. |
| Spec 3.2 kind K exact sha256 match and fuzzy match with `idmap_method` identity and coverage | **Not implemented.** The antigen and expression lookups use protein, gene, best transcript. 126 RS isoforms inherit another transcript's values. |
| Spec 3.6 and 4 per-sha256 cache | **Not used by the wrappers.** Adding one protein reruns the whole job. |
| Spec 4.1 rule 2 best reciprocal hits | Manual in Task 14; no code. |
| Spec 5 kappa, expected false calls, PPV | Kappa for the R0 reciprocal-hit comparison by hand (Task 14); expected false calls and assumed-prevalence PPV in the report by hand (Task 14, Task 15); not in code. |
| Sliding 80-residue window of the Codex rule | **Not implemented.** The field is a single local alignment and is named as such. |
| IEDB negative IgE assays; leave-species-out antigen test; SignalP 6 training-set overlap with the Phase C positives; protein-family annotation per allergen call | Owner decisions C2, C3 and a calibration task for later (not code here). |

Placeholder scan: none. Every code step has the complete file. The Pfam review (Task 13 step 3) uses the tested command `pfam-specificity`.

Known gaps, stated so the reviewer can check them:
1. The antigen calibration is circular (cut chosen after the anchors were seen; three anchors are one family). The plan records it and caps the status; it does not remove it.
2. No negatives exist for the allergen module. The leave-species-out run writes no status entry, so `iuis_allergen_homolog` stays `unvalidated`.
3. TMHMM, not Phobius, is the TM tool; neither has a calibration here. The `no_tm` window of 35 residues is a choice, not a measurement.
4. `hmmsearch --cut_ga` uses the Pfam gathering thresholds; they were not tuned for fungal proteins. Two proteins that UniProt names "CFEM domain-containing" were missed at GA in one proteome (microbiology review).
5. Run times of the whole-proteome Pfam, repeat and TMHMM jobs are not measured yet (Task 12 step 4 measures them). Scale to 5,000 proteomes would need batching: five short jobs per proteome conflicts with the 1 to 1.5 h job-size rule.
6. The IUIS set has 111 usable sequences (116 with a sequence; 4 peptide fragments and 1 free-text entry excluded). Onygenales is represented by four *Trichophyton* proteases, two of which are fragments that cannot be searched; none is from *Coccidioides*.
7. R0 status comes from UniProt gene models; a de novo annotation may change the call at the N-terminus (6.3% of sampled pairs).
8. The `tandem_repeat_protein` detectors have no gate against intracellular repeat families (ubiquitin, EF-hand, ankyrin, WD40, zinc finger). The adhesion call is gated by the signal peptide only.
9. The plan does not add Phobius, an ML step 1 variant, R1/R2 freezing, Fungi_5k batch runs, a stored golden `calls` file, or a PredGPI call.
