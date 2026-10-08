# Sorting hat: runs with PA14 active, a check against curated labels, and the overall scale

*2026-10-07. Exploratory. Everything here is descriptive. No call except `signal_peptide_protein[R0]` has a measured accuracy, and that one only where the status column says so. Data: `docs/reports/data/sorting_hat/classifier_check/` and `.../runs/`. Research use only; domain and similarity calls are evidence, not adhesin, allergen or antigen calls.*

## 1. What was run

Nine proteomes were run end to end with the family table of `main` (commit 9a3202c). Only PA14 (PF07691) is active. The other 14 families are inactive, so `wall_family_domain = not_called` means "no hit in an active family". It does not mean "no wall domain".

New in this round: *Blastomyces dermatitidis* ER-3 (Fungi_5k file, taxon 559297, 8,107 proteins; no truth set, no status, so every status is `unvalidated`). Its jobs took 2 to 12 minutes each (Pfam 2.4, repeats 2.2, BLAST 1.9, TMHMM 6.9, SignalP on a GPU 11.9).

## 2. Class breakdown per proteome

Counts of proteins called. `n/a` is `not_assessable` (S288C: 8 invalid sequences; the allergen-homolog, antigen and serodiagnostic calls need a module that is unavailable or limited to *C. immitis*).

| call | Af293 | A1163 | W72310 | Af293 Fungi_5k | S288C | S288C (no dubious) | *C. albicans* | *C. immitis* RS | *B. dermatitidis* |
|---|---|---|---|---|---|---|---|---|---|
| proteins | 9,647 | 9,942 | 10,556 | 9,161 | 6,722 | 6,039 | 6,212 | 9,910 | 8,107 |
| signal peptide (R0) | 753 | 751 | 783 | 694 | 310 | 296 | 367 | 460 | 366 |
| cell wall adhesion candidate | 8 | 6 | 9 | 8 | 13 | 13 | 15 | 7 | 8 |
| surface, no mechanism evidence | 745 | 745 | 774 | 686 | 297 | 283 | 352 | 316 | 358 |
| not surface (no signal peptide) | 8,879 | 9,171 | 9,758 | 8,455 | 6,392 | 5,729 | 5,838 | 8,208 | 7,730 |
| tandem repeat evidence | 23 | 26 | 24 | 20 | 25 | 19 | 22 | 27 | 19 |
| wall family domain (PA14 only) | 0 | 0 | 0 | 0 | 4 | 4 | 0 | 0 | 0 |
| allergen similarity (35%, 80 aa) | 104 | 104 | 115 | 105 | 84 | 84 | 64 | 84 | 69 |
| allergen homolog (70%, 80% coverage) | 41 | 39 | 43 | 41 | 15 | 15 | 13 | 16 | 14 |
| serodiagnostic marker candidate | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 143 | 0 |
| antigen rank top 15% | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1,376 | 0 |

Read these rows as follows.
- Signal peptide, candidate, no-mechanism and not-surface partition the assessable proteins (in Af293: 8 + 745 = 753 with a signal peptide, and 8,879 without).
- Between 4.5% and 7.8% of proteins have a signal peptide (S288C 4.6%, *B. dermatitidis* 4.5%, Af293 7.8%). Between 6 and 15 proteins per proteome are `cell_wall_adhesion_candidate`.
- The repeat call fires for 19 to 27 proteins per proteome. A large share has no signal peptide, so it is not a candidate (Af293: 23 repeat calls, 8 candidates).
- The antigen rows are meaningful for *C. immitis* RS only. The module is limited to that reference, so every other proteome shows 0 or `n/a`.
- Status: `signal_peptide_protein[R0]` is `smoke` for Af293, S288C and *C. albicans*. It is `unvalidated` for the other proteomes and for every other call.

## 3. Check against curated labels (exploratory)

`data/curated/adhesins/adhesins.tsv` classes joined to the calls (`analysis/sorting_hat_run/curated_label_check.py`). Proteins that tuned a module, and their homologs, are flagged by a name pattern of mine (ALS, FLO, SOWgp, PRA, RodA, CalA, CTS1, EPA, HWP1, IFF, HYR); this removes whole families. Some E3 rows are labels inherited from a shared Pfam domain. Two of the S288C "adhesins" (TDA8 and BSC1) say "no adhesion evidence" or "PA14 fragment ORF" in their own evidence field.

Strongest set (not tuned or homolog; adhesin E1 or E2, hard negative N1 or N2):

| proteome | class | n | signal peptide | repeat | PA14 domain | candidate |
|---|---|---|---|---|---|---|
| S288C | adhesin | 3 | 3 | 0 | 0 | 0 |
| S288C | hard negative | 27 | 26 | 5 | 0 | 5 |
| *C. albicans* | adhesin | 5 | 5 | 0 | 0 | 0 |
| Af293 | adhesin | 1 | 1 | 1 | 0 | 1 |
| Af293 | hard negative | 10 | 6 | 0 | 0 | 0 |

All rows, and the set without tuned proteins only, are in `curated_label_check.summary.tsv`. The same numbers for all rows:
- **S288C.** 13 adhesins: 11 with a signal peptide, 5 with the repeat call, 4 with PA14 (FLO1, FLO5, FLO9 and FLO10, all flagged tuned). 28 hard negatives: 27 with a signal peptide, 5 repeat calls, no PA14.
- ***C. albicans*.** 25 adhesins: 24 with a signal peptide (IFF11 has none), 5 with the repeat call. Of 56 indirect regulators, 1 has a signal peptide. 20 of 26 "surface other" proteins have one.
- **Af293.** 6 adhesins, 10 hard negatives (the beta-glucosidases have no signal peptide).

What this shows, and what it does not:
1. **R0 removes the intracellular proteins and keeps nearly every secreted wall protein.** Regulators: 1 of 56 called. Hard negatives: 27 of 28 called in S288C. R0 is a filter for surface candidates. It does not separate adhesins from other wall proteins.
2. **The mechanism calls do not separate them either, on proteins that did not tune a module.** In S288C the repeat call fires for 5 of 27 hard negatives (PIR1, HSP150, TIR1, TIR4 and FIT1, which the curated table marks as hard negatives) and for 0 of 3 strong-evidence adhesins. PA14 fires for none of the strong-evidence adhesins. The proteins the repeat call and PA14 do catch are the ones the detectors were tuned on (FLO1, FLO5, FLO9, FLO10).
3. **The samples are small.** 3 and 5 proteins in the strong adhesin sets, and 27 hard negatives in one species. These numbers are not an estimate of sensitivity or specificity. They are a reason not to claim one.
4. **The check is partly circular.** The detectors and domains were built from known adhesins, and the curated table was built from the same literature.

## 4. A genome with no truth set: *Blastomyces dermatitidis* ER-3

The distribution in section 2 is typical: 366 signal peptides (4.5%), 8 candidates, 19 repeat calls, no PA14, 69 and 14 allergen similarity and homolog calls. The BAD1 ortholog (`F00FD2C2_006066-T1`; 88.6% identity over 1,196 aa to UniProt A4D962, `data/.../bad1_vs_bder.tsv`) is called by both repeat detectors (period 24, about 27 copies, coverage 0.55), has a signal peptide, and is a `cell_wall_adhesion_candidate`. BAD1 was also a control in earlier detector work (`analysis/cocci_repeats/38_unit_bad1_control.py`), so this is a sanity check and not an independent test.

## 5. Overall scale

**The framing.** The tool was meant to take a proteome and return a breakdown of the classes it can identify.

**What it returns today.** The configuration defines ten calls (section 2). They are:
- a partition of the proteome by signal peptide: not surface, surface without mechanism evidence, and cell wall adhesion candidate (repeat or domain evidence plus a signal peptide);
- evidence calls: tandem repeat, wall family domain, allergen similarity, allergen homolog;
- *Coccidioides* only: antigen rank top 15% and serodiagnostic marker candidate.

Each call carries a status. Only `signal_peptide_protein[R0]` has a measurement. This is a labelling tool with evidence columns. It is not a classifier that assigns each protein to one mechanism class.

**Where it stands against that framing.**

| need | state |
|---|---|
| Run on a proteome, get per-class counts | Works. Nine proteomes in minutes per job. |
| Remove non-surface proteins | Works, with a measured status (`smoke` in three species, `estimated` in *A. nidulans*, recall 0.48 to 1.00 depending on species; the high values come from 7 to 19 positives). |
| Identify adhesin mechanism classes (repeat, domain families, hydrophobin, CFEM, Als-like) | Partial. The repeat call exists and is unvalidated. One domain family is active (PA14). 14 are waiting for review or sign-off. |
| Separate adhesins from other wall proteins | Not shown. Section 3. |
| Per-class accuracy | Not available for any mechanism call. Needs per-call status entries (spec ready, no code), truth tables from curation tasks 03, 04 and 07, and more independent positives. |
| Antigen and allergen classes | Allergen similarity works as a similarity flag. Antigen ranking works for one genome only. |

**What stands between today and a tool that reports class counts with stated reliability.**
1. Per-call status entries (plan, tests, code).
2. Truth tables: repeat clusters (19 against a floor of 20), reviewed Pfam member and non-member lists (issue #69), hard negatives (issue #14).
3. A decision on the remaining 14 Pfam families (CFEM and Hydrophobin are reviewed on 8 rows; the Candida families need more).
4. A rescue design for the proteins R0 does not call (recall 0.48 in *C. albicans*).
5. A re-run of all proteomes on the final table and a report.
The model-based step 1 (issue #15) is not on this list because decision 1 of 2026-10-07 defers it.

## 6. Limits of this document

- The curated-label check uses a name pattern for tuned proteins. A different pattern changes the counts.
- The `wall_family_domain` row depends on the number of active families (one).
- The Pfam, repeat and allergen modules use the data in the stored runs. The Af293 Fungi_5k and A1163 and W72310 runs have no status source, so their statuses are `unvalidated`.
- *B. dermatitidis* has no truth set, so section 4 is a distribution only.
