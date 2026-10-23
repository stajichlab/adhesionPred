# PR 6 hyd1_like build: held-out recall of a custom hyd1_like model

*2026-10-23. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_hyd1_like_2026-10-23/` (README there lists the jobs and files). Freeze: `analysis/wall_fusion_truth/pr6_hyd1_freeze.json` (signed 2026-10-23, #132). Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md`. Earlier ecm33 record: `2026-10-23-pr6-ecm33-build.md`.*

Status: every cell is `descriptive`. A custom model is report-only. No call, tier, panel, status or `narrow_ok` changes with this record.

## 1. Pre-check (job 30321490, 3 min 24 s, code `5fb6b48`)

Query A0ACE4KFX7 (KS0017) against the shared hard negatives, the four kept proteomes and the curated truth set. All numbers below were recomputed from the committed tables in `negcheck/`.

- 0 family members among the 161 shared negatives (phmmer E <= 1e-5, the freeze's labelled rule). `negatives_ignored` stays empty.
- One negative, B6QUR1 (MP1, Talaromyces marneffei, `laccase_etc`), hits the HsbA panel model at GA (score 168.7). It is not called family by phmmer, so it is a panel hit and not a family member.
- 0 proteome proteins are called family by the Kiss-protein phmmer (`proteome_labelled_phmmer.tsv` has only a header).
- The panel (eight hydrophobin-class models) hits 19 rows in the proteomes (`proteome_panel_ga.tsv`): 18 distinct proteins, because RODB_ASPFU hits both Hydrophobin and Eas (16 distinct Aspergillus fumigatus proteins and 2 Coccidioides immitis proteins). With B6QUR1 that is 19 `panel_ignored` ids, which are the 19 ids in the freeze.
- 7 curated truth proteins are related to the Kiss protein at E <= 1e-3 (`related.tsv`, scores 13.2 to 34.2). Removing the training member A0A6V8R0V1 leaves 6, the held-out external set. No weak hit lies between E 1e-5 and 1e-3 among the negatives.

## 2. Extract (job 30325384, 3 min 37 s, code `9f4843b`)

30 members: the query plus 29 UniRef50 representatives (drift-supported hits, query cover >= 0.5, reading `total`), in 4 named clades: Sordariomycetes 19, Dothideomycetes 7, Leotiomycetes 3, Eurotiomycetes 1 (the query, the Kiss protein). 3 hits without a named clade are listed in `excluded.tsv` and not counted. Units: Sordariomycetes `unit`, Dothideomycetes `unit`, Leotiomycetes `too_few_members`, Eurotiomycetes `kiss_single`. Two cluster units (30 percent MMseqs2 clusters, 10 held-out members each: C7ZQZ4 and A0A2K0UGM8) and the `external` line (6 curated truth proteins) were added by the freeze.

## 3. Build (array 30335767, 10 tasks, all COMPLETED exit 0, 11 to 37 s each)

`sbatch --array=0-9%4`, partition `batch`, 8 cores, 16 GB, 4 h, code `9f4843b`, freeze sha256 `65ed9bb6...552ff2`. 9 unit/mode rows plus the `all` model. Independence: `cluster30` is the headline (no training member in the same 30 percent / 50 percent coverage cluster as a held-out member), `none` removes only the held-out members.

**Sordariomycetes under `cluster30` is not buildable.** After the independence rule 4 training sequences remain (Dothideomycetes 3, Leotiomycetes 1) in 2 clades, below the freeze rule (5 sequences, 2 clades); 7 neighbours were excluded. Its `none` column was built.

Settings: held-out searches start only after `cutoff.json` exists; the cutoff takes no positive score; extra-call limit 16.2455 over 32,491 proteins; hard-negative limit 0.05.

## 4. Result per unit and mode

Columns: held out n; training size; excluded 30 percent neighbours; hits of the custom model (score at or above the cutoff), Pfam GA, relaxed hydrophobin panel (`data/sorting_hat/hydrophobin_relaxed.hmm`), seed-only phmmer, jackhmmer ceiling; u_indep (held-out members missed by every independent comparator); recovered of u; cutoff; extra calls over the four proteomes at that cutoff (panel_ignored ids left out; in brackets the calls at the same cutoff with them counted); the cutoff that would apply if the panel_ignored ids were counted.

| Unit | Mode | Held out | Train | Excl. | Custom | Pfam GA | Panel | phmmer seed | jackhmmer | u_indep | Rec. of u | Cutoff | Extra calls (counted) | Cutoff if counted |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Dothideomycetes | cluster30 | 7 | 15 | 8 | 7 | 0 | 5 | 0 | 0 | 2 | 2 | 10.2 | 6 (11) | 14.4 |
| Dothideomycetes | none | 7 | 23 | 0 | 7 | 0 | 5 | 0 | 0 | 2 | 2 | 12.3 | 3 (4) | 6.6 |
| Eurotiomycetes (Kiss) | cluster30 | 1 | 20 | 9 | 1 | 0 | 1 | seed in unit | 0 | 0 | 0 | 0.1 | 0 (2) | 7.5 |
| Eurotiomycetes (Kiss) | none | 1 | 29 | 0 | 1 | 0 | 1 | seed in unit | 0 | 0 | 0 | 0.1 | 0 (6) | 7.0 |
| Sordariomycetes | cluster30 | 19 | 4 | 7 | not buildable | | | | | | | | | |
| Sordariomycetes | none | 19 | 11 | 0 | 19 | 0 | 19 | 3 | 13 | 0 | 0 | 13.2 | 1 (1) | 6.4 |
| cluster C7ZQZ4 | cluster30 | 10 | 20 | 0 | 10 | 0 | 10 | 2 | 10 | 0 | 0 | 0.2 | 1 (3) | 0.2 |
| cluster A0A2K0UGM8 | cluster30 | 10 | 20 | 0 | 9 | 0 | 8 | 0 | 0 | 2 | 1 | 15.0 | 1 (5) | 7.9 |
| external (6 curated) | cluster30 | 6 | 8 | 22 | 5 | 3 | 5 | 0 | 0 | 1 | 1 | 17.0 | 1 (1) | 17.0 |
| all (30 members) | none | 0 | 30 | 0 | | | | | | | | 0.1 | 0 (6) | 8.0 |

Exact two-sided 95 percent Clopper-Pearson intervals (full set in `report/recall_by_unit.tsv`): Dothideomycetes cluster30 custom 7/7 0.590-1.000, Pfam GA 0/7 0.000-0.410, panel 5/7 0.290-0.963, recovered of u 2/2 0.158-1.000; Sordariomycetes none custom 19/19 0.824-1.000, phmmer seed 3/19 0.034-0.396, jackhmmer 13/19 0.434-0.874; cluster A0A2K0UGM8 custom 9/10 0.555-0.997; external custom 5/6 0.359-0.996, Pfam GA 3/6 0.118-0.882; Kiss unit 1/1 0.025-1.000. With these n a wide interval is the honest statement of how little each row shows.

## 5. Findings

(a) **No unit is testable.** u_indep is at most 2 in every row (Dothideomycetes u = 2, recovered 2; cluster A0A2K0UGM8 u = 2, recovered 1; external u = 1, recovered 1; all others u = 0), below `min_u` 3. The rule result is `not_testable` in all 9 rows.

(b) **The existing relaxed hydrophobin comparator already recovers most held-out members**: 19 of 19 Sordariomycetes (`none`), 10 of 10 in cluster C7ZQZ4, 8 of 10 in cluster A0A2K0UGM8, 5 of 7 Dothideomycetes, 5 of 6 external. The custom model adds little beyond it on these units. Seed-only phmmer is weak (0 of 7 Dothideomycetes, 3 of 19 Sordariomycetes). Pfam GA hits 0 of 7 Dothideomycetes and 3 of 6 external.

(c) **The Kiss unit (n = 1)**, the Kiss protein held out and the rest used for training, is recovered (score 57.7 under cluster30, 73.7 under `none`). n = 1 is a single protein.

(d) **Several cutoffs are at or near the floor** (0.1 bits for the Kiss unit and the `all` model, 0.2 for cluster C7ZQZ4), because the false-hit side allows almost any cutoff. This does not inflate recall here. Recomputed from `hits_custom.tsv` and `heldout_ids.txt` at the higher cutoff that applies when the panel_ignored ids are counted (`cutoff_if_panel_counted`), recall is unchanged in every buildable row: 7/7, 7/7, 1/1, 1/1, 19/19, 10/10, 9/10, 5/6 (same as at the frozen cutoff). Held-out members score well above both cutoffs: the Kiss protein 57.7 against 7.5; the lowest Dothideomycetes member 18.5 against 14.4 (cluster30); the lowest Sordariomycetes member 20.0 against 6.4. The two misses are one cluster-A0A2K0UGM8 member with no hit (score 0) and one external protein scoring 6.8 below its cutoff 17.0.

(e) **Every cell is descriptive.** Members are iterated-search hits, not curated hyd1_like proteins. The held-out external line is 6 curated hydrophobin truth proteins related to the Kiss protein by E <= 1e-3: tiny n, wide intervals. The custom model is report-only.

## 6. Open owner questions (not decisions)

1. Is a custom hyd1_like model worth shipping when the existing relaxed hydrophobin model recovers most of these members?
2. Should the inside-clade drift check be run for the Eurotiomycetes-clade hits and the 3 no-clade hits? The clade units are tiny and Eurotiomycetes has only the Kiss protein.
3. Should the model be integrated report-only, as Ecm33 was?

## 7. Not in this record

`train.aln.fa`, `model.hmm` and the member FASTA files stay on HPCC. No SLURM job was submitted for this record; the report ran locally.
