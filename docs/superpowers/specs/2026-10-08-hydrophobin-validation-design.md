# Hydrophobin call validation and 8-cysteine rescue rule: design

*2026-10-08. Revision 2 (after independent review 1, `2026-10-08-hydrophobin-validation-design-review-1.md`:
4 blockers and 12 major findings, all addressed below). Author: Claude Code (claude-sonnet-5-5). Owner
request of 2026-10-08: make the hydrophobin class validated, using the literature and the Pfam data, and
write the 8-cysteine pattern idea into the paper notes. Status: draft for a second review. No code is
written.*

## 1. Goal and scope

Goal: a hydrophobin call in `cellsurface_sorting_hat` with a measured status, so that the report can say
"this protein is a hydrophobin, with sensitivity X and specificity Y in species Z, with these limits".

**What status can be reached.** Every status file from this plan will be `smoke`, not `estimated`. The reason is
numeric (section 4.5): no species with a full proteome in scope has 20 independent positive clusters. A
`smoke` status adds a measured number and an interval to the report, and a note that the call was tested. It
does not give the reader a claim of reliability across species. The plan states this up front so that the
result is not read as more than it is.

In scope:
1. A hydrophobin call, separate from `wall_family_domain`.
2. A truth set for hydrophobins that does not depend on a Pfam match.
3. A rescue rule based on the 8-cysteine pattern for hydrophobins that no Pfam model finds. Whether the
   rescue is kept in the default call is decided by a rule fixed in advance (section 6.3).
4. Measurement of the Pfam-only call and the Pfam-plus-rescue call.

Out of scope: class I against class II as a measured call; structure prediction; model fine-tuning; changing
the Pfam gathering cutoff.

## 2. Facts that this design depends on

Numbers were measured on 2026-10-08. Items marked "(reviewer)" were re-derived by the independent reviewer.

| # | Fact | Source |
|---|---|---|
| F1 | Pfam 38.2 has seven hydrophobin-class models: PF01185 Hydrophobin, PF06766 Hydrophobin_2, PF22354 Eas, PF28987 DewD, PF29785 Hydrophobin_D, PF29802 Hyd1F, PF29465 Hydrophobin_like. HsbA (PF12296) is a related model. Clan CL0925 holds only PF01185, PF06766 and PF22354. A search of `Pfam-A.hmm.dat` for hydrophobin-related terms found no other model. | `Pfam-A.hmm.dat` (HPCC) |
| F2 | UniProt Swiss-Prot (reviewed), Fungi: 174 entries have "hydrophobin" or "rodlet" in the protein name (189 rows returned, 15 do not match the name). 165 carry one of the seven models **in the UniProt Pfam cross-reference** (PF01185 99, PF06766 47, PF22354 12, PF29802 2, PF28987 2, PF29785 2, PF29465 1; no entry has two). 9 carry none. The cross-reference is not an `hmmsearch` at Pfam 38.2 GA, and the two can differ. UniProt release 2026_03. (reviewer) | `analysis/hydrophobin_truth/` (`README.md` has the query) |
| F3 | The 9 entries without a model all match my screening regex (8 to 11 cysteines; `nopfam9_regex.tsv`). Best full-sequence hit to the 18 searched models, from `hmmsearch -E 1e-3`: HCF2_FULFL Eas 2.1e-8; PSH_FLAVE Hydrophobin 2.4e-8 (domain score 22.1 against a domain GA of 23.0); HYD1_GIBZE Eas 7.6e-7; HFBE_PENEN Hydrophobin 7.4e-6; HYD2_GIBZE Eas 9.9e-4. 4 entries have no hit at E <= 1e-3. | `nopfam9.domtbl`, `nopfam9_regex.tsv` |
| F4 | `categories.yaml` has no hydrophobin call. `wall_family_domain` is `pfam_adhesion.hit`, true for a hit to any active family. The `pfam_adhesion` module output has the columns `hit` and `families` only. | `categories.yaml`, `modules/pfam.py` |
| F5 | PR #75 (open, not merged) activates the hydrophobin models in `pfam_adhesion`. Once merged, a hydrophobin hit sets `wall_family_domain` and, with a signal peptide, `cell_wall_adhesion_candidate`. A hydrophobin is category 2c of the five categories, not a wall adhesin. | derived from F4 |
| F6 | Hits per proteome at `--cut_ga` (discovery search): seven hydrophobin-class models: Af293 7, A1163 6, W72310 8, *C. immitis* RS 2, *B. dermatitidis* ER3 2, S288C 0, *C. albicans* 0, H99 0. With HsbA: Af293 16, A1163 13, W72310 14, RS 2, ER3 4. The eight proteomes are five species (three are *A. fumigatus* strains). In Af293, 9 of the 16 hits are HsbA. (reviewer) | `_workdir/sorting_hat/calibration/hydrophobin_discovery/*.domtbl` |
| F7 | The only Swiss-Prot hydrophobins from the proteome species in scope are RodA to RodG of *A. fumigatus* (7 entries). All 7 carry a hydrophobin-class model in the UniProt cross-reference. None of the 9 entries without a model is from an in-scope species (they are *P. expansum*, *T. asperellum*, *H. virens*, *G. zeae* x2, *F. fulva* x3, *F. velutipes*). In Af293 the screening regex flagged 3 secreted proteins without a model; all three are annotated GPI-anchored proteins. (reviewer) | UniProt query; `analyse.py` output |
| F8 | Prior art for the 8-cysteine pattern, from PubMed abstracts: Yang 2006 used the C-CC-C-C-CC-C pattern with MEME and MAST, kept 9 candidates "after filtering by pattern, domain and length" (PMID 17217508); Kubicek 2008 and Seidl-Seiboth 2011 describe the eight conserved cysteines and Trichoderma hydrophobins that deviate in spacing (PMIDs 18186925, 21424760); Xu 2021 found 40 *P. ostreatus* genes by the pattern (PMID 33636611); a 2026 paper reports a surface-active protein, PAC3, without the motif (PMID 42546942). Yang 2006, Kubicek 2008 and Peñas 1998 have open full text. | PubMed, `docs/paper/05` section 5 |

## 3. Design decisions to make the call measurable

### 3.1 Dependencies

This design uses the per-call status code (PR #76, branch `per-call-status`) and the family-table changes of PR
#75. Neither is merged. H1 starts after both merge, or after the work branch is rebased on `per-call-status`.

### 3.2 A separate module and call

The hydrophobin-class family rows move from `pfam_adhesion` to a new Pfam module name `pfam_hydrophobin`
(`MODULES` in `pfam.py`, the module-set assertion in the tests, and the module list in
`outputs.py` change with it). The call table gets:

```yaml
- name: hydrophobin_domain
  expr: {flag: pfam_hydrophobin.hit}
```

This removes the F5 side effect. `wall_family_domain` keeps CFEM, PA14 and later wall families.

**HsbA and Hydrophobin_like (decision D6, before H1).** In Af293, 9 of 16 hits would be HsbA. The owner counts
both as hydrophobins for cataloging, but there is no hydrophobin truth label for HsbA in this design, so a
combined call would be mostly unlabelled in *A. fumigatus*. Default: put the seven hydrophobin-class models in
`pfam_hydrophobin` and HsbA in its own module `pfam_hsba` with its own call `hsba_domain`. The report and the
category 2c line combine them (`hydrophobin_domain` OR `hsba_domain`) as a display step, not as a measured
call. The owner decides whether HsbA is combined.

### 3.3 Rescue module

New module `cys8_pattern` (tool output, one row per protein). It reads the step 1 R0 result as a condition
table, in the way `pfam_adhesion` reads `second_condition` (`modules/cli.py` `_condition_table`). It does not
need the SignalP cleavage site, so the SignalP module is not changed and the R0 statuses stay valid.

| Field | Meaning |
|---|---|
| `hit` | 1 if the sequence matches the fixed pattern (section 5) **and** the R0 call for the protein is `called`. 0 if it matches and R0 is `not_called` or the pattern does not match. Empty (not assessable) if R0 is missing. |
| `spacing_class` | `I`, `II` or `other`, by which literature spacing matched |
| `n_cys` | cysteine count of the full sequence |
| `length` | full sequence length |

The pattern is run on the full sequence, not on a mature region.
The module `params` hold the pattern string, the spacing source, and any length window, so `params_hash`
changes when the rule changes and a call file for the old rule goes stale.

### 3.4 Call

```yaml
- name: hydrophobin_protein
  expr:
    or: [{flag: pfam_hydrophobin.hit}, {flag: cys8_pattern.hit}]
```

The call is ungated (no step 1 variant), reads two modules, has no `ref`, and is eligible for a call status
file (`call_eligible`). It can also be listed as a mechanism of the `other_*` calls. Its status depends on the
identity of the R0 module through `cys8_pattern`.

The Pfam-only measurement is the **module status of `pfam_hydrophobin`**, and `hydrophobin_domain` is
reported with that status. A one-module call status is not allowed (decision D3).

### 3.5 Effects on existing calls (decision D2)

1. `hydrophobin_protein` is a mechanism of `other_surface_no_mechanism` and `other_not_surface`, so a surface
   protein with a hydrophobin call is not "no mechanism". Both lists change. Default: add it to both.
2. `categories.yaml` changes `config_sha256`. The existing call files for `tandem_repeat_protein`
   (S288C, *C. albicans*) become stale. Task H1b re-runs `calibrate truth --call-status` for both and checks
   the numbers equal the earlier ones (0.435 and 0.462 sensitivity).
3. Both Pfam modules share one artefact digest over the whole family table (`modules/cli.py`, `_family_digest`).
   A change to a hydrophobin row changes the identity of `pfam_adhesion` as well, and the reverse. Task H1c
   gives each module a digest over its own rows only. If H1c is not done, the report says so.

## 4. Truth set

### 4.1 Rule that avoids circularity

UniProt names come partly from Pfam and InterPro rules. A positive therefore needs evidence that does not
depend on a domain match.

| Tier | Evidence | Use |
|---|---|---|
| T1 | A paper reports protein-level evidence: purification or mass spectrometry of the protein, rodlet or film formation, contact-angle or surface-tension change, or a **hydrophobin-specific** deletion phenotype (loss of rodlets, wettable conidia, loss of surface hydrophobicity). A generic growth or development phenotype does not count. The PMID is recorded. | primary positives |
| T2 | Swiss-Prot reviewed entry with a literature-backed function or location (evidence ECO:0000269) | primary positives |
| T3 | Name or family comes only from a rule or sequence similarity (ECO:0000255, ECO:0000250 without a paper, automatic annotation) | secondary set, reported separately, never in the primary result |

**What the labels do not remove.** Hydrophobins are named because they have the 8-cysteine pattern (F8). T1 and
T2 remove Pfam-derived labels, not pattern-derived ones. So the sensitivity of the rescue rule on named
hydrophobins measures how well proteins conform to the chosen spacing. It does not measure discovery of
unknown hydrophobins. The report says this beside the number.

### 4.2 Sources

1. UniProt Swiss-Prot reviewed, Fungi. The query string, release (2026_03) and date are in
   `analysis/hydrophobin_truth/README.md`. Evidence codes per annotation come from the UniProt JSON API.
2. PubMed for papers that purify or knock out a hydrophobin. Every added protein has a PMID.
3. Curated entries in the repository (`data/curated/`).

The counts of T1 and T2 entries and of clusters are not known yet. Task H2 measures them before any
threshold is chosen. A rough proxy from the review: Swiss-Prot entries with a PubMed citation in the function
comment are 14 for *P. ostreatus*, 9 for *B. bassiana* and 7 for *A. fumigatus*. This is not the T1/T2 tier.

### 4.3 Negatives

One negative set serves both variants. No label depends on a prediction.

1. **Per-species population (follows the repeat precedent, `docs/reports/2026-10-08-repeat-call-calibration.md`
   section 2).** Curated proteins plus reviewed, secreted UniProt proteins of the species mapped to the proteome
   by exact sequence, minus any protein with a hydrophobin label of any tier. These negatives are **assumed**
   (absent annotation is not measured absence). `docs/paper/03` section 3 says specificity is never inferred
   from absence in a database. The repeat work set the same exception, and the owner accepted it with the
   label "assumed". The same label is used here. Decision D8 asks the owner to confirm the exception for this
   call.
2. **Cross-species hard negatives** (CFEM, cerato-platanin, killer toxins, small secreted cysteine-rich
   proteins, PIR and Ccw12-like wall proteins, and the ATP8/ATP9 entries that the keyword search returned) are
   used only in a sequence-level report table, because most are not in an in-scope proteome.
3. A protein with a Pfam hit and no label is a negative for both variants. It is listed as a candidate false
   positive and reviewed by hand (section 6.2). It is not removed from the negatives.

### 4.4 Clustering, split and species

1. Cluster all truth proteins with MMseqs2 at 30% identity and 0.5 coverage, as in the repeat work. Record
   species. Bootstrap CIs by cluster.
2. **The development/test split is made in H2, before any measurement, and is unconditional.** Split by
   cluster, 30% development, 70% test, fixed seed, stored in a file. If the test part has fewer than 10
   positive clusters, the report says the split gives no usable test and the result is `smoke` with leakage
   `partial` (section 5.3). The pattern and any threshold are tuned on the development part only.

### 4.5 Species and what status they can reach

| Level | Species | Positives (Swiss-Prot hydrophobins) | Status |
|---|---|---|---|
| Status file, in scope now | *A. fumigatus* (Af293, A1163, W72310; one species) | 7 (RodA to RodG), all with a Pfam model | `smoke` at best |
| Sequence-level table, no status | all Swiss-Prot fungi (174 named entries) | 165 with a model, 9 without | report table |
| Extension E1, not yet in scope | *G. zeae* PH-1, *H. virens* Gv29-8, *P. expansum*, *F. fulva*, *P. ostreatus* PC15, *B. bassiana* ARSEF 2860, *T. asperellum* | the 9 no-model entries come from five of these | `smoke` at best, 0 to 28 named entries each |

In *A. fumigatus* the rescue cannot add a true positive: all 7 labelled entries already have a model (F7). The
rescue gain in a proteome can therefore be measured only in E1 species. Task H5b downloads the E1 proteomes (a
shared-storage and compute cost, listed in decision D7) and runs the modules on them. Without E1, the rescue
gain is reported only at sequence level on the Swiss-Prot set, and the call stays Pfam-only in the default
configuration.

`estimated` needs at least 20 positives, 20 negatives and 20 clusters per species (`docs/paper/03` section 2).
The highest count of named entries in one species is 28 (*P. ostreatus*), and those are from genome
annotation, so many may be T3. `estimated` is not expected anywhere.

## 5. The rescue pattern

### 5.1 Source of the parameters, and the order of tasks

The spacing comes from the literature and is **frozen in a commit before the truth set is built** (H3 before
H2). Papers to read: Wessels 1994, Linder 2005, Sunde 2008 (class I and class II spacing), Kubicek 2008,
Seidl-Seiboth 2011 and Yang 2006 (full text open for the last two and for Kubicek). H3 records for each number
the paper, the page and the sentence.

The two spacings that I know from memory are unverified and are not in any file. I did not run them on any
truth protein.

The regex used in the discovery search (`C.{3,12}CC.{8,45}C.{5,30}C.{3,12}CC.{4,25}C.{0,20}`) is a screening
tool I wrote. It is not the rule.

### 5.2 Other conditions

1. The R0 call must be `called` (section 3.3).
2. No cap on total length (decision D4, default). CFTH1 has three class II domains in one protein (PMID
   10336622). The match is per protein, and `hit` does not count domains. If a count of domains is wanted, it
   is a separate field `n_spans` and does not change `hit`.
3. The pattern itself has eight cysteines. A separate count condition is not added.

### 5.3 Leakage

The honest label is **`partial`**, with a note. Reasons: I read the 9 entries without a model and ran a regex on
them; I ran the regex on the 8 proteomes to be measured and saw the false-positive list; a no-length-cap
default also keeps PSH_FLAVE (515 aa), a protein already seen. Rules that follow:

1. The 9 entries are not used as unit tests. Unit tests use synthetic sequences built from the frozen
   spacing and proteins outside the truth set.
2. The 9 entries and their clusters stay in the measured set. They are marked as "seen" in the table, and
   the result is also reported without them.
3. Any status file records `leakage: partial` with a note naming the 9 entries and the discovery search.
   Status is capped at `smoke`. This costs nothing because every status is `smoke` anyway.

## 6. Measurement and outputs

### 6.1 Metrics

For `pfam_hydrophobin` (module status) and `hydrophobin_protein` (call status), per species group:

- Sensitivity and specificity with cluster-bootstrap 95% CIs, and counts of positives, negatives and clusters.
- Precision of the call. For the rescue: number of **rescue-only calls** (called by `cys8_pattern`, not by
  Pfam), their precision after manual review with a Wilson interval, and how many are labelled true.
- Number of `not_assessable` records (missing R0 or module rows).
- False-positive count per proteome.
- Sensitivity of the frozen spacing pattern on the 165 Pfam-positive Swiss-Prot entries. It shows how often the
  fixed spacing misses known members (see Seidl-Seiboth 2011 on deviating spacing).
- The per-model table: hits and labelled proteins for each Pfam model, including HsbA.
- T1/T2 and T3 results in separate tables.

### 6.2 Manual review

Every labelled protein that is not called (miss), every called protein with no label (candidate false
positive, including Pfam hits without labels and every rescue-only call), and every protein flagged "seen".
Each row records a BLAST top hit, length, pattern spacing and a decision. A label changes only with a recorded
reason and a PMID, as in the repeat work (`owner_decisions.tsv`).

### 6.3 Rule for keeping the rescue (fixed now)

The rescue stays in the default `hydrophobin_protein` call only if, **after manual review**, the rescue-only
calls over all measured proteomes number at least 5 and have a Wilson 95% lower bound for precision of at
least 0.5. If there are fewer than 5 rescue-only calls (including 0), the result is "no evidence that the
rescue helps". Then `cys8_pattern` stays in the tool and is reported, but it is not part of the default call.
This rule uses precision only after review, so it does not count unlabelled true hydrophobins as false.

## 7. Risks

1. **Few clusters.** See section 4.5. Every status is `smoke`.
2. **Assumed negatives hide true hydrophobins.** An unannotated true hydrophobin is counted as a false
   positive. Specificity is a lower bound. Manual review (6.2) bounds it.
3. **Swiss-Prot transfer.** The 28 *P. ostreatus* entries are from genome annotation and may be T3. Tiers
   are assigned per annotation.
4. **HsbA and Hydrophobin_like.** See decision D6.
5. **Non-8-Cys surface-active proteins** (PAC3, PMID 42546942) cannot be found. This is a stated limit.
6. **Basidiomycota and clade bias.** The truth is dominated by *Pleurotus*, *Flammulina*, *Agaricus*,
   *Trichoderma* and *Beauveria*. A leave-one-species-out table at sequence level shows transfer.
7. **Pfam and pattern disagree at the cutoff.** PSH_FLAVE is 0.9 bits below the domain gathering cutoff
   (F3). Pfam-only misses it by a margin that is not meaningful. This is reported, and the cutoff is not
   changed here.

## 8. Novelty and manuscript note

The 8-cysteine pattern is not novel (F8). A method paper cannot claim it. Yang 2006 already combined the
pattern with a domain and length filter, so the comparison with that paper must be made from its full text
before any claim. What could be new, and what this plan can test:
1. A measured gain over Pfam alone. The figure "9 of 174 Swiss-Prot hydrophobin-named entries (5.2%) have no
   Pfam model" is a count on names and UniProt cross-references. It is an **upper bound** on the gain in that
   set, not a measured gain.
2. A measured false-positive cost of the rescue in whole proteomes.
3. A hydrophobin call with a status and a leakage statement inside a tool that reports its own uncertainty.
   This is a property of the tool, not of the rule.

The search so far is PubMed only. Task H3b reads Yang 2006, Kubicek 2008 and Peñas 1998 in full (open access),
and extends the search to bioRxiv, InterPro entry notes and the Pfam release notes. Until then, "abstracts only"
stays in `docs/paper/05`.

## 9. Decisions

| # | Decision | Default |
|---|---|---|
| D1 | May a status entry cover a clade from a pooled measurement? | No. Species entries only. Pooled results are report tables. |
| D2 | Does `hydrophobin_protein` count as a mechanism in `other_surface_no_mechanism` and `other_not_surface`? | Yes, in both lists. |
| D3 | Allow a call status on a one-module call? | No. The module status covers `hydrophobin_domain`. |
| D4 | Cap mature length in the rescue? | No cap. `n_spans` is a separate field and does not change `hit`. |
| D5 | Exclude T3 labels from the primary result? | Yes. |
| D6 | Are HsbA and Hydrophobin_like in `pfam_hydrophobin`? | Seven hydrophobin-class models in `pfam_hydrophobin`. HsbA in `pfam_hsba`. Combine for display only. Owner decides. |
| D7 | Download the E1 proteomes (about seven fungal proteomes) to measure the rescue? | Yes, after the owner agrees. Without E1 the rescue is not measured in whole proteomes. |
| D8 | Confirm the "assumed negative" exception of the repeat work for this call? | Yes, labelled "assumed". |

## 10. Work plan

| Task | Content | Output |
|---|---|---|
| H0 | Correct the PF01185 `specificity_note` (it says 18 hits in 4 proteomes; the discovery tables give 20 hits in 5 proteomes, adding *B. dermatitidis* ER3). Commit `analysis/hydrophobin_truth/`. | data commit |
| H1 | After PR #75 and #76 merge: new module `pfam_hydrophobin` (and `pfam_hsba` by D6), `hydrophobin_domain`, `hydrophobin_protein`, `cys8_pattern` registration; update `MODULES`, `test_data_and_scripts.py` module-set assertion, `outputs.py` report text, golden files and `other_basis` strings in `test_cli.py`, and the D2 mechanism lists | code, tests pass |
| H1b | Re-run `calibrate truth --call-status` for `tandem_repeat_protein` in S288C and *C. albicans*; check the numbers match | status files |
| H1c | Per-module family digest | code, tests |
| H3 | Read the source papers; record the exact cysteine spacing with page and sentence; freeze the pattern in a commit | table in `docs/paper/05`, frozen commit |
| H3b | Read Yang 2006, Kubicek 2008 and Peñas 1998 in full; wider prior-art search | `docs/paper/05` |
| H2 | After H3: build the truth set, fetch evidence codes, assign tiers, add PMID-backed entries, cluster, make the unconditional split, record "seen" proteins | `analysis/hydrophobin_truth/` tables with provenance |
| H4 | Write the `cys8_pattern` wrapper and tests (synthetic sequences and proteins outside the truth set) | module, tests |
| H5 | Re-run `hmmsearch` with the new family table, then the wrappers, on the measured proteomes | module tables |
| H5b | If D7 is yes: download and run the E1 proteomes | module tables |
| H6 | Measure both variants, apply rule 6.3, manual review, write status files | report in `docs/reports/` |
| H7 | Update paper notes (ledger rows, limits, novelty) | docs |

Each of H1 to H4 has its own commit. A plan and an independent review of the plan follow this spec, before
code, as in earlier work.
