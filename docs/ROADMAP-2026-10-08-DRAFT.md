# Roadmap draft: cell wall protein identification from a proteome (2026-10-08)

*Draft. Written 2026-10-08 by a read-only survey (Claude Code, claude-opus-5-5). No file was edited and
nothing was committed. The owner has not reviewed it. Each statement cites a file, a commit or a
command. "Unverified" means I did not check it myself.*

Survey basis: branch `hydrophobin-validation` at `1201bc4` (it moved from `75a6c47` during the survey:
a second session committed task L6a at 22:08 and has untracked L7a job files in
`analysis/hydrophobin_truth/run/`). `origin/main` is `484247e` (PR #74).

## 0. Headline

1. The tool `cellsurface_sorting_hat` runs end to end. It ran on 9 proteomes, plus 5 proteomes for the
   hydrophobin work (`_workdir/sorting_hat/`, git-ignored; `docs/reports/2026-10-07-classifier-runs-and-scale.md`).
2. Only one call has a measured status on `main`: `signal_peptide_protein[R0]`. It is `estimated` in
   *A. nidulans* only and `smoke` in five other species (`docs/reports/2026-10-07-sorting-hat-run-and-calibration.md` §3).
3. Two more calls have `smoke` status, but only on unmerged branches: `tandem_repeat_protein`
   (S288C and *C. albicans*, PR #76) and `hydrophobin_domain` (8 proteomes, local branch only).
4. No call has any status in *Coccidioides*. R0 is `unvalidated` in *C. immitis* RS (no Onygenales truth).
5. SOWgp itself is flagged in RS today (`tandem_repeat_protein`, `cell_wall_adhesion_candidate[R0]`,
   `cocci_specificity_rank_top15`; `_workdir/sorting_hat/Cimm_RS/out/calls.long.tsv.gz`, protein
   `XP_001245172.2`). SOWgp tuned these modules, so this is not a test (panel row, `data/sorting_hat/panel.tsv`).
   The tool has no call that says "this looks like SOWgp" in another species.
6. About 57 commits of work (per-call status, hydrophobin, CFEM activation) exist only on a local,
   unpushed branch (`git rev-list --count origin/main..hydrophobin-validation` = 57 at `1201bc4`).

## 1. Git state (facts)

Commands: `git branch -a -vv`, `git rev-list --left-right --count origin/main...<b>`, `gh pr list`, `gh issue list`, `git worktree list`.

| branch | behind / ahead of origin/main | PR | state |
|---|---|---|---|
| `hydrophobin-validation` (checked out) | 0 / 57 | none | local only, not pushed. Contains `per-call-status` (merge `f005e20`) and 3 of the 4 commits of `signoff-cfem-hydrophobin` |
| `per-call-status` | 2 / 20 | #76 OPEN | CI green (`gh pr checks 76`). Contains `repeat-truth-curation` (merge `ddc9cbc`) |
| `signoff-cfem-hydrophobin` | 0 / 4 | #75 OPEN | CI green. Commit `1466f90` (PF01185 note fix) is not in `hydrophobin-validation` and conflicts with it in `data/sorting_hat/family_table.tsv` (`git merge-tree --write-tree hydrophobin-validation signoff-cfem-hydrophobin`) |
| `repeat-truth-curation` | 2 / 2 | none | in worktree `../adhesionPred-curation`; already merged into `per-call-status` |
| `main` (local) | 23 behind | | stale; fast-forward only |
| `pfam-proposals-and-curation-leads`, `sorting-hat-classifier-runs`, `sorting-hat-run-report`, `verify-literature-claims` | 0 ahead | #71-#74 merged | can be deleted |
| `origin/fig-repeat-viz-styling`, `origin/handoff-2026-10-03`, `origin/paper-notes`, `origin/task-08-repeat-other-species` | 0 ahead | merged | can be deleted |

Open issues: #9, #10, #12, #13, #14, #15, #16, #19, #26, #69 (`gh issue list --state all`).

## A. The pipeline as built: input FASTA to outputs

The five categories of version 1 are: surface glycoprotein, cell wall and adhesion candidate, antigen,
allergen, other (`docs/superpowers/specs/2026-10-04-orchestrator-design.md` §1, decision D2).
The calls are in `src/cellsurface_sorting_hat/categories.yaml`. Command:
`cellsurface_sorting_hat --fasta P.faa --taxon <id> --taxdump <dir> --workdir W --out O`. Module tables come
from `cellsurface_sorting_hat_module <name>` and the SLURM scripts in `scripts/sorting_hat/`
(`submit_modules.sh` submits SignalP, TMHMM, Pfam, repeats, allergen BLAST).

Outputs: `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`, `report.md`,
`run.json` (spec §3.6; seen in `_workdir/sorting_hat/Cimm_RS/out/`).

Key to columns: Impl = code exists; Tests = unit tests exist; Runs = proteomes run; Status = measured status;
Where = merged to `main` or branch only. "9 proteomes" = Af293 UniProt, Af293 Fungi_5k, A1163, W72310,
*C. albicans* SC5314, *C. immitis* RS, *B. dermatitidis* ER-3, S288C, S288C without dubious ORFs.

| # | Step | Tool, module | Call | Category | Impl | Tests | Runs | Status (numbers) | Where |
|---|---|---|---|---|---|---|---|---|---|
| 0 | Read FASTA, check IDs and residues, taxon | core `cli.py` | (`na_invalid`) | all | yes | yes | 9 | n/a (software) | main |
| 1 | Signal peptide | SignalP 6 GPU, `step1_rule@R0` (`modules/signalp.py`) | `signal_peptide_protein[R0]` | surface glycoprotein | yes | yes | 9 | `estimated` *A. nidulans* sens 0.688 [0.592, 0.779], spec 0.988; `smoke` S288C 0.848, *C. albicans* 0.477, *A. fumigatus* 0.947 (19 pos), *C. neoformans* 0.857 (7), *U. maydis* 1.000 (9). Leakage: SignalP 6 training overlap not measured. `unvalidated` for RS, A1163, W72310, Fungi_5k, *B. dermatitidis* | main |
| 1b | R1, R2, ML step 1 | `analysis/step1_compare/` | `[R1]`, `[R2]`, `[ML]` variants | surface glycoprotein | analysis only | yes (637 collected in `tests/step1_compare`) | Phase C | `unavailable` in the tool. Owner decision 2026-10-07: R0 only, ML deferred | main (analysis) |
| 2 | TM helices | TMHMM 2.0c, `tm` (`modules/lookups.py`) | none (evidence; `no_tm` condition for CFEM) | evidence | yes | yes | 9 | `unvalidated` | main |
| 3 | Wall family domains | `hmmsearch --cut_ga`, Pfam 38.2, `pfam_adhesion` + `family_table.tsv` | `wall_family_domain` | cell wall and adhesion | yes | yes | 9 | `unvalidated`. Active on main: PA14 only. Active on branch: PA14, CFEM (`no_tm`) | main (CFEM: PR #75 and branch) |
| 4 | Hydrophobin domain | same search, `pfam_hydrophobin` (7 models) | `hydrophobin_domain` | cell wall (not part of `cell_wall_adhesion_candidate`; only in the `other_*` mechanism list) | yes | yes | 9 + 5 | `smoke` in 8 proteomes: sens 1.00 (Af293 7 pos, A1163 6, W72310 5, *B. bassiana* 7, *P. ostreatus* 3), 0.80 *P. expansum*, 0.60 *F. graminearum*, 0.25 *F. fulva*; spec 0.9994-1.0000 on assumed negatives; leakage `partial` | branch only |
| 5 | HsbA | `pfam_hsba` | `hsba_domain` | cell wall (other_* list) | yes | yes | 9 + 5 | `unvalidated`, no truth set | branch only |
| 6 | Hydrophobin relaxed level | `hmmsearch --nobias`, 2.6 bits, `hydrophobin_relaxed.hmm`, >= 8 Cys, R0 (`modules/hydrophobin_relaxed.py`, commit `1201bc4`) | none yet (`hydrophobin_extended` waits for the ship rule, plan L6c) | cell wall | yes | yes | L7a in progress (untracked sbatch) | L5: recovers 6 of 6 Pfam-missed T2 clusters, Wilson [0.61, 1.00]; 130 of 131 T2; 0 hard-negative calls; `partial` leakage | branch only |
| 7 | 8-Cys pattern | `cys8_pattern` (`modules/cys8.py`, `cys8_spacing.yaml`) | none | cell wall | yes | yes | 12 | rule 6.3: no evidence the rescue helps (1 labelled hydrophobin gained in 8 proteomes); spacing matches 78 of 131 T2 | branch only |
| 8 | Tandem repeats | `02_repeat_profile.py` and `14_repeat_detect_general.py` via `repeat02`, `repeat14` | `tandem_repeat_protein` | cell wall and adhesion | yes | yes | 9 | `smoke` (owner definition, 2026-10-08): S288C sens 0.435 [0.067, 0.667], spec 0.987 (23 pos / 16 clusters); *C. albicans* 0.462 [0.000, 0.788], 0.976 (13 / 6). Leakage `tuned_on_truth`. Negatives assumed. No measurement in filamentous fungi | PR #76 and branch |
| 9 | Composite | engine | `cell_wall_adhesion_candidate[R0]` = (repeat OR wall domain) AND R0 | cell wall and adhesion | yes | yes | 9 | derived (weakest deciding status), not a measurement (`docs/paper/04` §2) | main |
| 10 | Allergen similarity | BLASTP vs 111 WHO/IUIS sequences, `allergen_homology` | `iuis_allergen_similarity` (35% / 80 aa), `iuis_allergen_homolog` (70% / 80% cov or allergen Pfam) | allergen | yes | yes | 9 | `unvalidated`. Leave-species-out recall 65/111 and 41/111; no negatives. `pfam_allergen` families inactive | main |
| 11 | Antigen ranking lookup | precomputed `analysis/cocci_antigens/cocci_antigen_ranking.tsv`, `antigen_lookup` | `cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate[R0]` | antigen | yes | yes | RS only | `unvalidated`; ranking prints NOT CALIBRATED (3 of 4 anchors). Calls 1,376 RS proteins (15%); 143 with R0 | main |
| 12 | Cys-rich tier, spherule expression | lookups `cys_rich`, `expression` | none (evidence columns) | evidence | yes | yes | RS only | no accuracy claim | main |
| 13 | Other | engine | `other_not_surface[R0]`, `other_surface_no_mechanism[R0]` | other | yes | yes | 9 | derived | main |

Not in the tool (checked in `src/cellsurface_sorting_hat/modules/`): GPI anchor calls (PredGPI was run only
in `analysis/step1_compare/jobs/j1_features.sh`), Phobius, any ESM model, a Bys1 call (PF04681 inactive), ALS,
Flo11, GLEYA, Hyr/Iff, PIR and flocculin families (inactive rows), any SOWgp-type call, enzyme classes.

Tests (collect-only, `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest --collect-only -q`, on `1201bc4`):
`tests/cellsurface_sorting_hat` 818, `tests/step1_compare` 637, `tests/hydrophobin_truth` 83,
`tests/cys_candidates` 67, `tests/calibration_truth` 18, `tests/surface_glyco` 25 (+3 collection errors,
missing `esm` and others on this node). CI passes on PRs #75 and #76 (`gh pr checks`). The 2026-10-06
handoff said the package CI job had never run; that is no longer true.

Per-class counts on real proteomes (main config, PA14 only active; `docs/reports/2026-10-07-classifier-runs-and-scale.md` §2):
R0 4.5% to 7.8% of proteins; `cell_wall_adhesion_candidate` 6 to 15 per proteome; repeat 19 to 27.
With the branch config (CFEM active), RS has 14 candidates instead of 7, and Ag2/PRA (`XP_001240075.1`) becomes
`wall_family_domain` and `cell_wall_adhesion_candidate` (`_workdir/sorting_hat/Cimm_RS/out_hyd/calls.long.tsv.gz`).

## B. SOWgp and look-alikes

### B1. Original question and intended method

- Goal: "find antigenic candidates in Coccidioides building on SOWgp and PRA1/PRA3", and test presence and
  sequence variation across the pangenome (`docs/reports/2026-09-27-coccidioides-antigen-findings.md` §2).
- Methods used: orthology and prevalence in 488 *Coccidioides* proteomes; absence in four confounder fungi and
  human; spherule against mycelium RNA-seq (Carlin et al. 2021, PMID 34067070); a period-based repeat detector;
  structure surveys for small Cys-rich antigens (same report §3; `docs/TOOL-ARCHITECTURE.md` §3).
- Design insight that still holds: SOWgp and BAD1 are tandem-repeat proteins like FLO11, but Pro/Cys-rich, not
  Ser/Thr-rich. A repeat detector can find them; a composition-trained classifier cannot
  (`docs/TOOL-ARCHITECTURE.md` §2, table "A finding that refines class 2a").

### B2. Inventory

| Item | Location | Key numbers |
|---|---|---|
| Antigen ranking (genus specificity) | `analysis/cocci_antigens/`, `cocci_antigen_ranking.tsv` (9,139 rows) | Tier 1: 14, Tier 2: 45. Acceptance 3 of 4 anchors (PRA2 fails) -> NOT CALIBRATED. SOWgp prevalence 0.92, 0% confounder identity; excluded from Tiers 1 and 2 (no Fungi_5k SignalP match) |
| SOWgp repeat structure | `analysis/cocci_repeats/REPORT_2026-09-29_sowgp_repeat_structure.md` | 47 aa unit from anchor PTDCYGDC; length = base + 47 x units; immitis 2-5 units, posadasii 2-6 |
| Anchored family search | `REPORT_2026-09-30_anchored_family_search.md` | present in 92.3% of 493 proteomes (union); 81% of member hits are short alleles with whole-unit deletions |
| Presence/absence by read depth (PR #43) | `REPORT_2026-10-01_sowgp_depth_vs_repeats.md`, `CIMG_04613_coverage_summary.md` | 39 of 488 strains with no gene model have normal depth (0.92-1.43): gene-model gap. 7 strains with depth ratio < 0.5 at >= 10x: not tested at read level |
| Depth audit (PR #60) | `REPORT_2026-10-04_sowgp_depth_audit.md` | reproduces tables; low array-depth slope is not a pipeline artefact (0.110 immitis, 0.016 posadasii). Cause of the low slope: not investigated (`docs/HANDOFF-2026-10-07.md` §5) |
| Long-read locus check | `REPORT_2026-10-01_sowgp_longread_annotation.md` | 3 of 6 loci split into two gene models; miniprot joins them |
| SOWgp unit HMM, relatives | `sowgp_unit.hmm`, `34_sowgp_unit_hmm.py`, `REPORT_2026-09-30_sowgp_architecture_relatives.md` | 6,313 proteomes: all 1,293 hits with score >= 26 are *Coccidioides*; best other hit 19.7. 74 Onygenales Pro+Cys-rich look-alikes (none with a 47 aa period). No ortholog outside *Coccidioides* at unit level |
| Repeat surface proteins (class 2a curation) | `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md`, `class2a_candidates.tsv` (41), `class2a_candidates_general.tsv` (58) | 20 Pro/Cys-rich, 10 Ser/Thr-rich, 11 other. Only CIMG_04613 of the Pro/Cys-rich RS loci is spherule-induced. Unidentified Ser/Thr period-17 family; a Pro-rich family with no RS gene model |
| BAD1 control | `38_unit_bad1_control.py`; BAD1 ortholog in *B. dermatitidis* called by both detectors (classifier-runs report §4) | Whether script 38 was run: unverified |
| Spherule table | `analysis/cocci_spherule/spherule_surface_table.tsv.gz` (9,757 genes); report `2026-10-03-cocci-spherule-surface-table.md` | 564 up at 48 h; SOWgp +9.09 log2FC (top); only 23 of 564 with SignalP; 41 Cys-rich; no CFEM among the 564 |
| Spherule follow-up | `2026-10-03-cocci-spherule-followup.md`, `followup_genes.tsv` (45) | 19 genes with no hit outside *Coccidioides*; first leads CIMG_12808, CIMG_03569, CIMG_13082; gene models weak |
| Exon support (STAR) | `20_star_align.sh`, `21_exon_support.py` | STAR job 29386162 FAILED (BAM sort memory; log `_workdir/cocci_spherule/logs/star.29386162.log`). `21_exon_support.py` never run |
| Cys-rich (PRA3-like) finder | `analysis/cys_candidates/` (PR #36) | RS: 21 `cys_rich_sp_unassigned`. SOWgp fails its length and window rules by design |
| PRA3 structure, PF28404 | `docs/reports/2026-10-02-pra3-fulllength-structure.md`, `2026-10-02-pf28404-family.md` | PRA3 core 67 aa, no fold match; PF28404 in 283 of 831 proteomes, not genus-specific |
| CFEM (Ag2/PRA, PRA2) | `docs/reports/2026-09-29-class2b-structure.md` | CFEM hemophore fold (TM 0.807 to 4Y7S). 7 CFEM in RS |
| Truth/controls | `data/sorting_hat/panel.tsv` (SOWgp rows, `tuned`); `data/curated/antigens/antigens.tsv` (86 rows); `docs/agent-tasks/06-serodiagnostic-antigen-controls.md` | No independent truth for the antigen calls. Task 06 not started |

### B3. Where it stands

- Last decisions (owner, 2026-10-07; `docs/HANDOFF-2026-10-07.md` §2): antigen layer stays *Coccidioides*-only
  (decision 10); the STAR job and the SOWgp low-slope cause are side projects (decision 12); curation order puts
  task 06 and task 01 (Onygenales truth) last (decision 9).
- Unfinished: STAR rerun and exon support; rename or define the "specific" column (open since
  `docs/HANDOFF-2026-10-03.md` §4 item 1; I did not see it fixed); *C. posadasii* expression (no data in the repo);
  the 7 low-depth strains; the Onygenales truth spec (task 01); task 06 counts.

### B4. Is there a path from user proteins to "this looks like SOWgp"?

No, not as a named call. What exists:

1. In *C. immitis* RS, SOWgp is called by `tandem_repeat_protein`, `cell_wall_adhesion_candidate[R0]`,
   `cocci_specificity_rank_top15` and `serodiagnostic_marker_candidate[R0]` (RS run). All are `unvalidated`, and
   SOWgp tuned the repeat and antigen modules.
2. In other species, a SOWgp-like protein can only appear as a generic `tandem_repeat_protein` plus signal
   peptide. Nothing separates a Pro/Cys-rich array (SOWgp, BAD1 type) from a Ser/Thr-rich array (FLO11 type).
3. The antigen calls are `not_assessable` outside RS (lookup by ID).

Missing pieces:

| Gap | What exists to build from | Size |
|---|---|---|
| A definition: SOWgp ortholog, SOWgp-architecture (secreted Pro/Cys-rich tandem array), or spherule-surface protein | `docs/TOOL-ARCHITECTURE.md` §2-3 property table | owner decision |
| Ortholog call: `sowgp_unit.hmm` as a module (score >= 26 separates *Coccidioides* in 6,313 proteomes) | `analysis/cocci_repeats/34_sowgp_unit_hmm.py` | small |
| Architecture call: repeat call AND composition (Pro+Cys, Cys fraction) AND R0 | repeat module output already has the region; the 74 look-alikes list | medium (spec, freeze, review) |
| Truth for the architecture call | SOWgp alleles (tuned), BAD1 (tuned), 74 Onygenales look-alikes (unreviewed) | medium-large; best status `smoke` |
| R0 status in Onygenales | task 01 | large (curation) |

## C. Other cell wall classes: scope and state

Scope source: orchestrator spec §1, §3.7, D2, D11; `docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md`.

| Class | Scope decision | Family / tool | State |
|---|---|---|---|
| Secreted and wall proteins in general | in (category 1) | R0 | measured (section A row 1) |
| GPI-anchored wall proteins as a class | no `cell_wall_protein` call in v1 (D11) | PredGPI only in step 1 analysis | no module; `curated_gpi.tsv` has no literature rows |
| CFEM (Ag2/PRA, PRA2, Rbt5, Csa1/2, CfmA-C) | in (`wall_family_domain`) | PF05730, `no_tm` | active on branch and PR #75 (owner 2026-10-08); status `unvalidated` |
| PA14 (Flo1/5/9/10, Epa) | in | PF07691, `signal_peptide` | active on main (PR #73); 25 of 29 hits are false domain hits removed by the condition |
| Hydrophobins | in (category 2c) | 7 Pfam models in `pfam_hydrophobin`; relaxed level | `smoke`; module split only on branch; relaxed level waits for L7a/L8/ship rule |
| HsbA | in as separate evidence (`hsba_domain`) | PF12296 | active; no truth (owner decision 4 of the hydrophobin handoff open) |
| Cerato-platanin | out | PF07249 | dropped from table (PR #75 body) |
| Bys1 / CalA | in, inactive | PF04681 | 17 hits, 16 `needs_expert`; calB/calC control not identified (`docs/reports/2026-10-07-pfam-family-proposals.md`) |
| ALS, Flo11, GLEYA, Hyr/Iff (PF11765), PIR | in, inactive | PF11766, PF05792, PF10182, PF10528, PF11765, PF00399 | proposals written; PIR non-members and HWP1/EAP1 miss need an expert |
| Flocculin, Flocculin_t3, Hyr1, PIR1-like_C, ALS_M | candidates, inactive | PF00624, PF13928, PF15789, PF22799, PF30910 | added on branch, no review |
| Tandem-repeat adhesins (FLO11, ALS, AGA1, HWP1, Iff/Hyr) | in (`tandem_repeat_protein`) | repeat02/14 | `smoke`, sens about 0.44; 20 misses are arrays the detectors do not find |
| Small Cys-knot (PRA3) | evidence only (`cys_rich` tier, RS only) | Cys-rich finder | not a classifiable class (`docs/TOOL-ARCHITECTURE.md`) |
| Moonlighting (Hsp60, gp43) | out | none | `known_miss` in panel |
| Glucanases, chitinases, GEL/GAS, yapsins, synthases, GPI biosynthesis | later (spec §3.7) | GH18, GH72, GH16, Asp, chitin synthase families | not started |
| Signaling, septation, non-protein components | out | none | out of scope |
| Allergens | in (category 4) | BLAST + PF16541, PF25312 (inactive) | `unvalidated`; issue #19 |
| Biofilm | out (no phenotype data) | none | not built |

## D. Open work and decisions (de-duplicated)

Size: small = under one session; medium = one to three sessions; large = more. Sizes are guesses (section F).

### D1. Branches and PRs

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 1 | Merge PR #76 (per-call status, repeat-call `smoke`) | owner go-ahead | none; CI green | small |
| 2 | PR #75: as pushed it puts hydrophobin models in `pfam_adhesion`, so hydrophobins set `wall_family_domain` (`docs/paper/04` §4 item 7). The branch supersedes 3 of its 4 commits | owner: close #75 and carry `1466f90` into the hydrophobin branch, or merge #75 and resolve the `family_table.tsv` conflict later | #76 | small |
| 3 | Push `hydrophobin-validation` (57 commits local only) | owner go-ahead | none | small |
| 4 | PR for `hydrophobin-validation` | owner; after hydrophobin owner stop 3 or as is | 1, 2 | small |
| 5 | Delete merged branches (local and remote, list in section 1); remove worktree `../adhesionPred-curation`; fast-forward local `main` | owner | 1 | small |
| 6 | Issue #26: untracked files at repo root | owner | none | small |

### D2. Hydrophobin extended level

Plan `docs/superpowers/plans/2026-10-08-hydrophobin-extended-level.md` rev 2 (not re-reviewed); spec
`docs/superpowers/specs/2026-10-08-hydrophobin-custom-hmm-design.md` rev 4; freeze `f1986e2`; L5 done
(`docs/reports/2026-10-08-hydrophobin-extended-level-L5.md`); L6a done (`1201bc4`).

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 7 | Owner stop 2: HMM decision. L5 says "not testable" (u = 0); HMM stays out unless the owner decides | owner | L5 | small |
| 8 | L7a proteome runs and cost per 10,000 proteins | none (in progress, untracked sbatch) | 7 | small |
| 9 | L8 evidence sheets and extra-call clusters | none | 8 | medium |
| 10 | Owner stop 3: decisions on the sheets (`owner_decisions.tsv`) | owner | 9 | medium (owner time) |
| 11 | L8b ship rule, L6c call `hydrophobin_extended` (only if ship), L6b re-measure repeat call files, L7b status, L9 docs | none | 10 | medium |
| 12 | Earlier hydrophobin decisions: rescue left out (default) or data-derived spacing; review 22 unlabelled Pfam calls; HsbA truth set or not | owner | none | small each |
| 13 | Not done: cross-species hard-negative table, Wessels 1994 / Linder 2005 / Sunde 2008 | none | none | small |
| 14 | Should `hydrophobin_domain` feed `cell_wall_adhesion_candidate`? Today it does not (`categories.yaml`) | owner | none | small (decision), medium (code and golden tests) |

### D3. Repeat call

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 15 | ALS7 override E1 -> E3 | confirmed by owner (`docs/HANDOFF-2026-10-08.md` §2a) | done | none |
| 16 | Derived composite status marker in `status_basis` for every run | owner (decision 4 of 2026-10-08) | #76 | small |
| 17 | Detector work for short-unit arrays (EAP1, PGA18, AGA1) and S/T arrays (FLO11) | owner decision 5; held-out truth set first | 18 | large |
| 18 | Held-out repeat truth; curation tasks 03 and 08; reach 20 clusters (about 60 for half-width 0.10) | none | task 07 | large |
| 19 | Repeat call in filamentous fungi (no array positives in *A. fumigatus*) | truth source | 18 | medium |

### D4. Pfam families (issue #69, task 04)

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 20 | Bys1: identify calB/calC; expert pass on 16 `needs_expert` | owner/expert | none | medium |
| 21 | ALS, Flo11, GLEYA, PF11765, PIR, new candidate rows: expert pass, sign-off | owner | task 07 | medium |
| 22 | AltA1: an *Alternaria* proteome or leave out | owner | none | small |
| 23 | Reviewed member lists to replace draft lists | expert | none | medium |
| 24 | Truth tables so CFEM and PA14 can get a status (today `unvalidated` by decision 4) | curation | 23 | medium |

### D5. Step 1 gate (R0 vs ML)

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 25 | R0 only for v1 (decided 2026-10-07). Inspect the 92 Phase C positives R0 misses before any rescue | none | none | medium |
| 26 | Rescue module per clade (e.g. GPI, Phobius-only SP) | spec | 25 | large |
| 27 | Onygenales R0 truth (task 01) | owner go (spec kept, decision 12) | none | large |
| 28 | More positives in Phase C species (task 02) | none | none | medium |
| 29 | SignalP 6 training overlap (leakage) | none | none | small-medium |
| 30 | ML model work: #9, #10, #15, #16 | deferred | 25 | large |
| 31 | Basidiomycota truth: #50 closed "not planned"; `curated_basidiomycota.tsv` (58 rows) read by no step | owner | none | stale |

### D6. Allergen and antigen

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 32 | Task 05: count IEDB IgE-negative records; decide if a specificity is possible | none | none | small first step |
| 33 | Panel rows missing: Asp f 34, Als1, Rbt5, a cross-species allergen, *Histoplasma* Hsp60 | none | none | small |
| 34 | Task 06: counts of *Coccidioides* proteins with negative assays; leave-species-out feasibility | none | none | small first step |
| 35 | Exact sequence sha256 ID mapping for lookups (126 RS isoforms inherit values) | spec | none | small |
| 36 | Antigen outside RS (*C. posadasii*, other genera) | owner (decision 10 says no for now) | confounder genomes | large |

### D7. SOWgp / Coccidioides side thread

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 37 | Define "SOWgp-like" for the tool (section B4) | owner | none | small |
| 38 | Spec for a SOWgp ortholog call and a Pro/Cys repeat architecture call | spec, review | 37 | medium |
| 39 | Rerun STAR with more memory; run `21_exon_support.py`; check by hand | none | none | small-medium |
| 40 | "specific" column wording in the spherule table and first report | owner | none | small |
| 41 | 7 low-depth strains; low-slope cause | none (side project) | none | medium |
| 42 | *C. posadasii* spherule expression | data | none | unknown |

### D8. Hard negatives and docs

| # | Item | Needs | Depends on | Size |
|---|---|---|---|---|
| 43 | Task 07 / issue #14 hard negatives (first in the curation order) | none | none | medium |
| 44 | Issue #12, #13 (validation gaps, positive set) | | 43 | large |
| 45 | Refresh stale docs (section E) | none | 1-4 | small |
| 46 | Re-check "copied" numbers in `docs/paper/02`; read the three cell wall reviews in full | none | none | medium |
| 47 | Enzyme classes (spec §3.7) | owner, later | family tables | large |

## E. Risks, gaps and contradictions

| # | Statement in a document | What the code and data show | Which is right |
|---|---|---|---|
| 1 | Orchestrator spec header: "No code, data or job exists for this spec" | Package exists on main since PRs #63, #64 | code. Header is stale |
| 2 | `docs/model-review/STATUS.md` (2026-10-02): orchestrator "proposal only, no spec"; families "all inactive" | spec, code and runs exist; PA14 active on main; CFEM and hydrophobins active on branch | code. STATUS.md is the most stale overview, and it is the file the handoffs tell readers to read first |
| 3 | `docs/TOOL-ARCHITECTURE.md` §5: "2c hydrophobin solved by existing HMMs; no work needed" | Pfam misses 9 of 131 T2 hydrophobins and RodD; relaxed level was needed | the hydrophobin reports |
| 4 | `TOOL-ARCHITECTURE.md`: "2a repeat/avidity adhesin works, PR-AUC 0.94-0.98" | those numbers belong to the old ESM classifier (orchestrator spec §2 row "Repeat detectors"); the repeat call measures sens 0.44 / 0.46 | the 2026-10-08 report |
| 5 | `docs/paper/03` §7: "No other module has a status above unvalidated" | repeat call and `pfam_hydrophobin` are `smoke` on branches | branch data (not yet on main) |
| 6 | Orchestrator spec §3.4 table: hydrophobins inside `wall_family_domain` | branch: separate `hydrophobin_domain`; main: hydrophobin rows inactive; PR #75: inside `wall_family_domain` | three states exist at once. The branch is the newest owner intent |
| 7 | `docs/HANDOFF-2026-10-03.md`: STAR job running | job failed (log) | log; corrected in `HANDOFF-2026-10-07.md` |
| 8 | `REPORT_2026-09-30_pfam_confirmation.md`: FLO11 tandem array, coverage near 1.0 | detector coverage 0.045 for FLO11 | 2026-10-07 run report §8 |
| 9 | Hydrophobin handoff lists the plan as next steps | L6a was committed after it (`1201bc4`); a second session is active | git log |
| 10 | Memory note `project_sorting_hat_plan2_implemented.md` | out of date (MEMORY.md says so) | repo |

What a user would be misled by today:

1. **`cell_wall_adhesion_candidate` is not an adhesin call.** In S288C the repeat call fires for 5 of 27 hard
   negatives and 0 of 3 strong-evidence adhesins that did not tune a module
   (`2026-10-07-classifier-runs-and-scale.md` §3). With CFEM active, Ag2/PRA becomes a "cell wall adhesion
   candidate", though its fold is a hemophore and adhesion is not shown (`TOOL-ARCHITECTURE.md` §2).
2. **"Surface glycoprotein" is a signal peptide call.** Glycosylation and surface exposure are not assessed
   (spec known limit 2).
3. **`not_called` for `wall_family_domain` means "no hit in an active family"**, not "no wall domain".
4. **Specificities near 1.0 rest on assumed negatives** (repeat and hydrophobin reports). They are lower bounds of unknown size.
5. **Every status in *Coccidioides* is `unvalidated`**, including R0. The tool's main target species has no measured call.
6. **The antigen call selects 15% of the RS proteome** (1,376 proteins); the cut was set after the anchors were seen.
7. **Results depend on which branch is installed.** Main, PR #75 and the hydrophobin branch give different
   `wall_family_domain` and `other_*` counts.

Gaps in measurement for categories the tool reports: allergen (no negatives), antigen (no independent truth), CFEM and
PA14 (no truth table), HsbA (no truth), composite calls (derived only), repeat call outside yeasts (none).

Process risk: the thread is spread over 8 handoffs, 3 overview files (STATUS.md, TOOL-ARCHITECTURE.md,
PLAN-2026-09-30) and the paper notes. None of the three overview files reflects the state after 2026-10-06.

## F. Updated roadmap

### F1. Target

"A user gives a proteome and gets a cell-wall-protein breakdown with measured sensitivity and specificity per
class, including SOWgp-like proteins." Under the project's own rules (`docs/paper/03` §2), most classes can reach
only `smoke` with existing truth. `estimated` needs at least 20 independent positive clusters and a half-width of
0.10 or less. So the realistic target for v1 is: every reported class has a status entry (`smoke` or better) in
at least two species, with stated leakage, and the report says where it has none.

### F2. Milestones in dependency order

Size in working sessions (one session is about one assistant day with owner check-ins). These are guesses. The
basis is the pace in the handoffs: per-call status (spec to tested code) took about one day
(`HANDOFF-2026-10-07` to `-08`); hydrophobin L0 to L5 took about one day (`git log` timestamps on 2026-10-08).
Curation took longer than planned in the Basidiomycota pilots (`HANDOFF-2026-10-03` §2a). Owner review time is not counted.

| M | Deliverable | Acceptance check | Owner decisions | Size |
|---|---|---|---|---|
| M0 | Clean base: #76 merged; #75 closed or merged; hydrophobin branch pushed and in a PR; merged branches deleted; STATUS.md and spec header refreshed | `git rev-list --count origin/main..per-call-status` = 0; `gh pr list --state open` shows only the hydrophobin PR; `pytest tests/cellsurface_sorting_hat -q` passes on main | merge order; fate of #75 | 1 |
| M1 | Hydrophobin class finished: L7a to L9; `hydrophobin_extended` shipped or rejected by the ship rule | `analysis/hydrophobin_truth/ship_decision.json` exists; L6b repeat numbers reproduce (0.435 / 0.462) | stop 2 (HMM), stop 3 (sheets), HsbA truth, whether hydrophobins count in the cell wall composite | 2-3 |
| M2 | Class list fixed for v1: which calls form the "cell wall breakdown" (repeat, CFEM, PA14, hydrophobin, HsbA, Bys1, ALS, SOWgp-like) and how they combine | a spec revision of `categories.yaml` design with independent review | the class list; whether `cell_wall_adhesion_candidate` is renamed (e.g. to a neutral "wall mechanism evidence") | 1 |
| M3 | Shared hard negatives (task 07, #14) | `counts.md` in the task PR; owner accepted | accept the set | 2 |
| M4 | Pfam families: Bys1, ALS, Flo11, GLEYA, Hyr/Iff, PIR reviewed and signed off or left out; reviewed member lists; truth tables for CFEM and PA14 | `family_table.tsv` rows with `active_by`; `calibrate truth --module pfam_adhesion` status files in two or more species | each family sign-off | 3-5 |
| M5 | SOWgp-like calls: (a) ortholog call from `sowgp_unit.hmm`; (b) Pro/Cys repeat architecture call; truth from alleles, BAD1 and reviewed look-alikes | a status file for each call (expected `smoke`, leakage `tuned_on_truth`); RS run calls SOWgp; *B. dermatitidis* run calls BAD1 under (b) only | definition (D7 #37); whether BAD1 counts as a positive for (b) | 3-4 |
| M6 | Repeat call improved on held-out truth (tasks 03, 08) | status file from a truth set no detector saw; sens interval reported per species | whether detector work is worth it | 4-6 |
| M7 | Onygenales R0 truth (task 01), so *Coccidioides* gets a status | `status/step1_rule@R0.json` with an Onygenales taxon | go/no-go on curation | 4-8 |
| M8 | Allergen and antigen: task 05 and 06 counts; status where possible, else a stated "cannot be calibrated" | `counts.md` per task; report text | whether to continue these categories in v1 | 2 |
| M9 | Final runs on all proteomes with the final config; one report with a per-class table (value counts, status, sensitivity, specificity, leakage, species) | `cellsurface_sorting_hat` on the 14 proteomes; a script that checks every report number against a stored table (as L9 plans) | release name and scope | 2 |
| Later | R0 rescue (92 misses), GPI module, ML step 1 (#9, #10, #15, #16), enzyme classes, antigen outside RS, biofilm | | | large |

Where the hydrophobin work sits: it is M1. It is the most complete class (truth set of 174 Swiss-Prot entries,
frozen procedure, L5 measured). It should finish before new classes start, because L6b re-measures the repeat
call files after the last `categories.yaml` edit.

Critical path to the target: M0 -> M2 -> (M3 -> M4) and M5 in parallel -> M9. M6 and M7 raise status quality; they
are not needed to report a `smoke` status per class, but without M7 nothing in *Coccidioides* has a status.

### F3. Merge and branch plan

1. Merge PR #76 (`per-call-status`) into main. The PR #76 body says it merges cleanly with main (checked with `git merge-tree`); CI is green.
2. Do not merge PR #75 as is. Its hydrophobin rows sit in `pfam_adhesion`, which the branch replaced. Recommended:
   cherry-pick the note fix `1466f90` onto `hydrophobin-validation` (resolve `family_table.tsv` once), then close #75
   as superseded. Alternative (the hydrophobin handoff's order #76, #75, branch): merge #75, then merge main into the
   branch and resolve the same conflict. Both give the same end table; the first avoids a period where main counts
   hydrophobins in `wall_family_domain`. This is the owner's choice.
3. Push `hydrophobin-validation` now as a backup (it is the only copy of 57 commits). Open the PR after owner stop 3,
   or now as a draft.
4. After step 1: delete `repeat-truth-curation` and its worktree, `per-call-status`, the four merged local branches
   and the four merged remote branches; fast-forward local `main`.
5. New work (M2, M5) starts from main after the hydrophobin PR merges, one branch per milestone.

### F4. The five most important next actions

1. **Owner: decide the merge order** (merge #76; close or merge #75) and allow a push of `hydrophobin-validation`.
   This removes the single-copy risk and the three-way difference in `wall_family_domain`.
2. **Finish hydrophobin M1** through owner stop 3 (L7a is already running in another session).
3. **Owner: define "SOWgp-like"** (ortholog, Pro/Cys repeat architecture, or both), then write a spec for M5. Today
   no call answers the SOWgp question outside RS.
4. **Fix the class list (M2)** and refresh `docs/model-review/STATUS.md` and the spec header, so that one current
   file states what the tool reports and what each status means.
5. **Start task 07 hard negatives** (M3). It is first in the owner's curation order and unblocks the Pfam sign-offs
   and any specificity that is not built on assumed negatives.

### F5. Decisions the owner must make first

1. Merge order and fate of PR #75.
2. Whether hydrophobins (and HsbA) are part of the cell wall composite call.
3. What "SOWgp-like" means for the tool.
4. Whether to fund Onygenales R0 truth curation (task 01). Without it, no call has a status in *Coccidioides*.

## G. Unverified in this survey

- Whether `38_unit_bad1_control.py` and `39_unit_distribution.py` were run.
- Whether the "specific" column was renamed after 2026-10-03 (no commit found by `git log` on the report since `dac3ce9`; not checked further).
- `cfem_onygenales.tsv` (415 CFEM proteins in 71 Onygenales genomes) is untracked; not checked that it exists.
- The full text of issues #10, #13, #15, #16, #19.
- Test counts are collect-only; tests were not run in this survey.
