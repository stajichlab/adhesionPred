# Hydrophobin call validation and 8-cysteine rescue rule: design

*2026-10-08. Revision 1. Author: Claude Code (claude-sonnet-5-5). Owner request of 2026-10-08: make the
hydrophobin class validated, using the literature and the Pfam data, and write the 8-cysteine
pattern idea into the paper notes. Status of this document: draft for independent review. No code
is written.*

## 1. Goal and scope

Goal: a hydrophobin call in `cellsurface_sorting_hat` with a measured status, so that the report can say
"this protein is a hydrophobin, with sensitivity X and specificity Y in species Z".

In scope:
1. A hydrophobin call, separate from `wall_family_domain`.
2. A truth set for hydrophobins that does not depend on Pfam.
3. A rescue rule based on the 8-cysteine pattern for hydrophobins that no Pfam model finds.
4. Measurement of the Pfam-only call and the Pfam-plus-rescue call, with status files.

Out of scope: fine class assignment (class I against class II) as a measured call (section 4.4 gives it as
an optional later step); structure prediction; other surface-active proteins (HsbA is in the family table as
hydrophobin by owner decision, but is measured only as a part of the Pfam set, see section 7.3).

## 2. Facts that this design depends on

All numbers below were measured on 2026-10-08 unless a source is named.

| # | Fact | Source |
|---|---|---|
| F1 | Pfam 38.2 has seven hydrophobin-class models: PF01185 Hydrophobin, PF06766 Hydrophobin_2, PF22354 Eas, PF28987 DewD, PF29785 Hydrophobin_D, PF29802 Hyd1F, PF29465 Hydrophobin_like. HsbA (PF12296) is a related model. Clan CL0925 holds only PF01185, PF06766 and PF22354. A search of `Pfam-A.hmm.dat` for hydrophobin-related terms found no other model. | `Pfam-A.hmm.dat`, release on HPCC |
| F2 | UniProt Swiss-Prot (reviewed), Fungi: 174 entries have "hydrophobin" or "rodlet" in the protein name. 165 carry one of the seven models (PF01185 99, PF06766 47, PF22354 12, PF29802 2, PF28987 2, PF29785 2, PF29465 1). 9 carry none. | `analysis/hydrophobin_truth/sp_hydrophobin_query.tsv` |
| F3 | The 9 without a model all contain the 8-cysteine pattern (my generous regex). 5 have a sub-cutoff hit (E-value 7e-6 to 1e-3) to Eas or Hydrophobin. 4 have no hit at E <= 1e-3. | `analysis/hydrophobin_truth/nopfam9.domtbl` |
| F4 | `categories.yaml` has no hydrophobin call. `wall_family_domain` is `pfam_adhesion.hit`, which is true for a hit to any active family. The `pfam_adhesion` module output has the columns `hit` and `families` only. | `src/cellsurface_sorting_hat/categories.yaml`, `modules/pfam.py` |
| F5 | Because of F4, activating the hydrophobin models on 2026-10-08 (PR #75) makes a hydrophobin hit set `wall_family_domain`, and so `cell_wall_adhesion_candidate` when the protein has a signal peptide. A hydrophobin is category 2c of the five categories, not a wall adhesin. | derived from F4 |
| F6 | Status entries are per species (one species per entry), and `estimated` needs at least 20 positives, 20 negatives and 20 clusters. A proteome has few hydrophobins (1 to 6 in the eight proteomes searched, more in some Basidiomycota). | `docs/paper/03`; the discovery search of 2026-10-08 |
| F7 | The 8-cysteine pattern is textbook knowledge. PubMed abstracts show: Yang 2006 used the C-CC-C-C-CC-C pattern with MEME and MAST to find 9 new candidates (PMID 17217508); Kubicek 2008 and Seidl-Seiboth 2011 describe the eight conserved cysteines, and Seidl-Seiboth 2011 describes Trichoderma hydrophobins that deviate in cysteine spacing (PMID 18186925, 21424760); Xu 2021 found 40 hydrophobin genes in *P. ostreatus* by this pattern (PMID 33636611). A 2026 paper reports a surface-active fungal protein (PAC3) that has no 8-cysteine motif (PMID 42546942). Abstracts only. | PubMed, 2026-10-08 |

## 3. Design decisions to make the call measurable

### 3.1 A separate module and call

Add a family-table module name `pfam_hydrophobin`. The seven hydrophobin-class rows and HsbA move from
`pfam_adhesion` to `pfam_hydrophobin`. `pfam.py` gets `pfam_hydrophobin` in `MODULES`. The call table gets:

```yaml
- name: hydrophobin_domain
  expr: {flag: pfam_hydrophobin.hit}
```

This removes F5: hydrophobins no longer set `wall_family_domain`. `wall_family_domain` keeps CFEM, PA14 and the
later wall families. `other_surface_no_mechanism` and `other_not_surface` must list the hydrophobin call as a
mechanism or not, and that is decision D2 (section 9).

### 3.2 Rescue module

New module `cys8_pattern` (tool output, one row per protein). Fields:

| Field | Meaning |
|---|---|
| `hit` | 1 if the sequence matches the fixed pattern (section 5) |
| `spacing_class` | `I`, `II` or `other`, by which of the two literature spacing patterns matched |
| `n_cys` | cysteine count of the mature region |
| `mature_length` | length after the signal peptide |

The module uses the SignalP 6 result of step 1 to cut the signal peptide. The module reads that table, so it
needs a dependency declaration like the one `pfam_adhesion` has for `second_condition`. It does not read any
truth file.

### 3.3 Calls

```yaml
- name: hydrophobin_domain        # Pfam only
  expr: {flag: pfam_hydrophobin.hit}
- name: hydrophobin_protein       # Pfam OR rescue
  variants:
    pfam_only: {ref: hydrophobin_domain}
    pfam_or_cys8:
      or: [{ref: hydrophobin_domain},
           {and: [{flag: cys8_pattern.hit}, {ref: signal_peptide_protein}]}]
```

Both variants have at least two modules, as the per-call status rules need (section 6a of
`docs/paper/03`). `pfam_only` reads one module. It is a module status case, not a call status case. The
call status machinery needs an expression of two or more modules, so `pfam_only` is measured as the status
of the module `pfam_hydrophobin`, and only `pfam_or_cys8` gets a call status file. Decision D3 asks
whether to allow a one-module call status. The default is no.

### 3.4 Optional later step: class I against class II

`spacing_class` can drive a split call. It is not measured in this plan. Reason: class II is mostly
Sordariomycetes, and a split truth set would have fewer clusters per class than the total.

## 4. Truth set

### 4.1 Rule that avoids circularity

UniProt protein names come partly from Pfam and InterPro rules. A truth set made from the name
"hydrophobin" can therefore measure Pfam against Pfam. A positive must have evidence that is independent
of a domain match.

Evidence tiers for a positive:

| Tier | Evidence | Use |
|---|---|---|
| T1 | A paper reports protein-level evidence: purification or mass spectrometry of the protein, rodlet or film formation, contact-angle or surface-tension change, or a deletion phenotype (for example loss of rodlets, wettable conidia), with the PMID recorded | primary positives |
| T2 | Swiss-Prot reviewed entry with a literature-backed hydrophobin description (UniProt evidence code ECO:0000269 on function or subcellular location) | primary positives |
| T3 | The name or family comes only from a rule or from sequence similarity (ECO:0000255, ECO:0000250 without a paper, or an automatic annotation) | secondary set, reported separately |

The primary measurement uses T1 and T2. The T3 set shows how much the result depends on this choice. If the
T1 plus T2 set is smaller than 20 clusters, status stays `smoke` and the report says why.

### 4.2 Sources

1. UniProt Swiss-Prot reviewed, Fungi (taxon 4751), with the query in `analysis/hydrophobin_truth/`.
   Evidence codes per annotation are fetched with the UniProt JSON API, not from the TSV.
2. Literature mining by PubMed for papers that purify or knock out a hydrophobin. Each added protein
   carries a PMID, and the abstract only is read unless full text is open access.
3. Curated entries in the repository (`data/curated/`), where a hydrophobin is already listed.

Work not done yet: the number of Swiss-Prot entries in T1/T2 and the number of clusters. The plan (section
10, task H2) measures them before any threshold is chosen.

### 4.3 Negatives

1. Hard negatives, by name and by family, from Swiss-Prot reviewed fungal entries: CFEM proteins,
   cerato-platanin and other small secreted cysteine-rich proteins (killer toxins, small secreted effectors
   with 4 to 12 cysteines), cysteine-rich wall proteins (PIR, Ccw12-like), and the ATP9/ATP8 proteins that
   the keyword search returned.
2. Bulk negatives: all other proteins of the measured proteomes. These are assumed negatives (annotation is
   absent, not a measured absence). The report states this and treats specificity as a bound.
3. Exclusion: a protein with a hydrophobin name, a T3 label, or a hit to any of the eight Pfam models at
   the gathering cutoff is not allowed in the negative set of the Pfam-only variant. For the
   `pfam_or_cys8` variant, a Pfam hit is a prediction, not a label, so such a protein is a negative only if
   it has no hydrophobin label of any tier. This difference is intended and is shown in the report.

### 4.4 Clustering and species

Cluster all truth proteins with MMseqs2 at 30% identity and 0.5 coverage, as for the repeat truth set.
Record the species of every protein. Bootstrap CIs by cluster, as in the repeat work.

Status files are per species. Species with many T1/T2 positives are rare. The plan therefore reports two
things:
1. Per-species measurement for species that have a full proteome and at least 3 T1/T2 positives (candidates:
   *A. fumigatus* Af293, *A. nidulans*, a *Trichoderma* species, *P. ostreatus*, *B. bassiana*).
   These get status files. Most will be `smoke`.
2. A pooled measurement over all species, with leave-one-species-out, to show whether the rule transfers.
   A pooled result is a report table and is not written as a status file, unless decision D1 allows a clade
   status entry.

## 5. The rescue pattern

### 5.1 Source of the parameters

The spacing must come from the literature before any truth protein is looked at. Source papers to read:
Wessels 1994, Linder et al. 2005, Sunde et al. 2008 (class I and class II cysteine spacing), Kubicek 2008 and
Seidl-Seiboth 2011 (deviations). Task H3 records the exact spacing from these papers and the paper,
page and the sentence for each number. Until then the repository holds two numbers from my memory that are
unverified: class I roughly C-X(5-9)-CC-X(11-39)-C-X(8-23)-C-X(5-9)-CC-X(6-18)-C-X(2-13) and class II
roughly C-X(9-10)-CC-X(11)-C-X(16)-C-X(8-9)-CC-X(10)-C-X(4-5). They are not used until checked.

The regex used in the discovery search of 2026-10-08
(`C.{3,12}CC.{8,45}C.{5,30}C.{3,12}CC.{4,25}C.{0,20}`) was written by me as a generous pattern. It is only a
screening tool. It is not the rule.

### 5.2 Other conditions of the rule

1. Signal peptide called by SignalP 6 (the step 1 gate).
2. Mature length window from the literature (hydrophobins are about 100 +/- 25 residues, PMID 9758836 gives
   this range). Multi-domain hydrophobins exist (CFTH1 has three class II domains in one protein, PMID
   10336622), so the window must allow a long protein with repeated domains, or the rule counts one
   hydrophobin domain per protein and does not cap length. Decision D4.
3. At least 8 cysteines in the matched span.

### 5.3 Leakage and tuning

The pattern and the window are fixed from the literature. They are not tuned on the truth. Leakage is
`none` for this rule. The sub-cutoff Pfam hits in F3 do show that a protein can be near the Pfam cutoff. A
second variant that lowers the Pfam E-value cutoff is not part of this design, because the cutoff choice
would need the truth.

If a change to the pattern is made after seeing truth results, then:
1. The truth set is split by cluster into a development part (30%) and a test part (70%) before any
   result is seen. The split is made once and stored with a seed.
2. A change is tested on the development part only.
3. The final measurement uses the test part. The leakage field is `none` with a note that the development
   part was used to tune. The pooled result without the split is also reported, with leakage
   `tuned_on_truth` (status cap `smoke`).

## 6. Measurement and outputs

For each variant (`pfam_only`, `pfam_or_cys8`) and each species group, report:

- Sensitivity and specificity with cluster-bootstrap 95% CIs, counts of positives, negatives and clusters.
- For the rescue rule: how many true positives it adds over Pfam-only, and how many negatives it adds as false
  positives. This is the net effect and decides whether the rescue is kept. A criterion that is set before
  the result: keep the rescue only if the added false positives in the bulk negatives are fewer than the added
  true positives, per species and pooled. The reason: bulk negatives hold many cysteine-rich secreted
  proteins (CFEM-like proteins, GPI-anchored proteins) that the regex found in the discovery search.
- The same table for T1/T2 and T3 positives.
- A confusion list: every protein that is a hydrophobin by evidence and not called (misses), and every called
  protein with no label (candidate false positives). The list is reviewed by hand with the PMID. The
  review result can change a label only with a recorded reason and a PMID, as in the repeat work
  (`owner_decisions.tsv`).

Status files follow the rules in `docs/paper/03` section 6a. Entries record `measure`, `leakage`, `taxa`
and the module identities.

## 7. Risks

1. **Few clusters.** The T1/T2 truth may give fewer than 20 clusters. Then the call stays `smoke` in every
   species. The report says "not measurable to estimated" and gives the count.
2. **Assumed negatives hide true hydrophobins.** Unannotated proteomes can contain hydrophobins that no
   curator has labelled. The rescue rule would count them as false positives. Specificity is then a lower
   bound. A manual review of the rescue-only calls in each measured proteome (listing ID, length, pattern,
   BLAST top hit) bounds this.
3. **Circularity through Swiss-Prot transfer.** Many *P. ostreatus* Swiss-Prot entries (28) are from genome
   annotation. They may be T3. Tier assignment must be per annotation, not per entry type.
4. **Hydrophobin_like and HsbA.** The owner counts them as hydrophobins for cataloging. HsbA has no
   published hydrophobin class assignment in the sources read here. Measurement treats them as part of the
   Pfam set. A per-model table (hits and labelled proteins per model) shows each model's contribution. If a
   model adds only unlabelled hits in a measured proteome, the report flags it for review.
5. **Non-8-Cys surface-active proteins.** PMID 42546942 reports a surface-active protein without the
   motif. The rule cannot find such proteins. This is a stated limit, not an error.
6. **Basidiomycota.** Class I hydrophobins are most numerous there, and the project stopped curation of
   Basidiomycota earlier (docs/HANDOFF-2026-10-03.md). The truth set may be biased to a few species
   (*Pleurotus*, *Flammulina*, *Agaricus*). The leave-one-species-out table shows this.

## 8. Novelty and manuscript note

The 8-cysteine pattern is not novel. It defines the family (F7). A method paper cannot claim the pattern.
What could be new, and what the measurement can test, is:
1. A measured sensitivity gain over Pfam alone (F2 and F3 show 5% of named Swiss-Prot hydrophobins have no
   Pfam model).
2. A measured cost in false positives of the rescue rule across proteomes.
3. A hydrophobin call with a status and a leakage statement inside a tool that reports its own
   uncertainty.

The paper notes record the prior art (PMIDs in F7), the abstract-only reading, and the instruction to read
Yang 2006 in full before any novelty claim. A claim of novelty needs a wider search than the one made on
2026-10-08 (PubMed only; no Google Scholar, bioRxiv, InterPro entry notes or the Pfam release notes were
read). The notes go to `docs/paper/05-literature-verification.md` (section 5) and
`docs/paper/04-limits-and-open-questions.md`.

## 9. Decisions

| # | Decision | Default |
|---|---|---|
| D1 | May a status entry cover a clade (for example Pezizomycotina) from a pooled measurement? | No. Species entries only. Pooled results are report tables. |
| D2 | Does `hydrophobin_domain` count as mechanism evidence in `other_surface_no_mechanism`? | Yes. A surface protein with a hydrophobin hit has a mechanism, so it is not "no mechanism". |
| D3 | Allow a call status for a one-module call? | No. The module status covers it. |
| D4 | Length rule of the rescue: cap the mature length, or count domains? | Count domain-like spans. No cap on total length. |
| D5 | Truth tiers: T3 excluded from the primary result? | Yes. |

## 10. Work plan

| Task | Content | Output |
|---|---|---|
| H1 | Move hydrophobin rows to `pfam_hydrophobin`, add `hydrophobin_domain`, update tests, golden files and the mechanism lists | code on a branch, tests pass |
| H2 | Build the truth set: fetch evidence codes, assign tiers, add PMID-backed literature entries, cluster, record species | `analysis/hydrophobin_truth/` tables with provenance |
| H3 | Read the source papers, record the exact cysteine spacing with page and sentence | table in `docs/paper/05` |
| H4 | Write the `cys8_pattern` module and tests (unit tests with the 9 no-Pfam proteins and class I and II controls) | module, tests |
| H5 | Run the modules on the measured proteomes (HPCC job) | module result tables |
| H6 | Measure both variants, write status files, review misses and rescue-only calls by hand | report in `docs/reports/` |
| H7 | Update paper notes (ledger rows, limits, novelty) | docs |

Each of H1 to H4 has its own commit. H4 comes after the owner accepts the spacing in H3. The spec, plan and
an independent review precede code, as in earlier work.

## 11. Not in this design

- Fine-tuning a model. The rule-based call is measured first.
- Changing the Pfam gathering cutoff.
- Fungal species outside the proteomes already in `_workdir/sorting_hat/` (a later extension for
  *Trichoderma*, *P. ostreatus* and *B. bassiana* needs their proteomes downloaded).
