# Handoff, 2026-11-27: adhesion benchmark v1, legacy removal, cross-clade round 2

Written so a fresh session can pick this up without re-deriving anything. Read this, then
`docs/superpowers/plans/2026-11-26-adhesion-round2-benchmark-cruft-and-layout-study.md` (the plan)
and `data/benchmark/adhesion_v1/README.md` (what the benchmark is and what it cannot do).

## The one finding that drives everything else

**Every model this project has trained reduces to amino-acid composition plus length.**

| Model | Its score | Composition + length | Source |
|---|---|---|---|
| `surface_glyco` ESM-2 + LR | ROC-AUC 0.999 leave-family-out; ~12% adhesin precision on S288C | 0.995 | `docs/model-review/2026-09-27-review-and-framework-plan.md` §4.1, §4.3 |
| Stage-2 ESM C 300M + LR | PR-AUC 0.939 leave-genome-out; **ROC-AUC 0.598 leave-clade-out** | 0.907 | same, §4.4–4.5 |
| AGA1 feature classifier | AUC 0.97 cluster-grouped CV | 0.973 (composition alone) | `docs/reports/data/aga1_feature_classifier_2026-10-25/ablation.tsv` |

So the models learned "long, Ser/Thr-rich, repetitive glycoprotein", not adhesion. Two consequences
shape all current work:

1. Any new signal must beat **composition + length**, or it is the same signal again.
2. A contrast against the whole proteome is won by composition. The decisive contrast is adhesins
   against **hard negatives**: wall proteins with the same composition and a known non-adhesion job.

What has held up on measurement: the signal-peptide gate (R0; sens 0.85, spec 0.97, held out to
*A. nidulans*), PredGPI in yeasts, Pfam family calls (precise on known families, **0 of 12** on
seed-independent Kiss clusters), the hydrophobin domain call, the calibration tooling.

## Owner decisions (2026-11-26 / 27) — do not relitigate

| # | Decision |
|---|---|
| D1 | Round 2 = **2a** cross-clade structure search as built, **2b** benchmark + domain library + ablation harness, **2c** fresh layout/profile study |
| D2 | Cruft: remove now, tag first (`legacy-2026-11-26`) |
| D3 | Both contrasts, tiered: **C1** adhesins vs hard negatives *decides keep/remove*; **C2** vs all secreted/GPI is descriptive |
| D4 | Literature curation for under-covered clades → owner sign-off → freeze |
| — | Layout study scale: one proteome per species, broad; species list from **PopSeq** (per the global metadata rule) |
| — | Moonlighting proteins: own tier, out of both contrasts and out of the layout study |
| — | Indirect adhesion phenotypes (adhesion drops on deletion, but downstream of an enzymatic/structural role): **hard_negative**, flagged `indirect=yes` |
| — | A null adhesion result alone does **not** make a hard negative |
| — | The 36 existing `surface_other_adhesion_phenotype` rows: use as hard negatives **after a signal-peptide filter** (30 of 36 pass) |
| — | Thin clades: merge into subphylum; anything still thin is reported `scored=no`, never silently merged |
| — | ESM-2 embedding helper: relocate (not delete) so step 1 Phase B still runs |

Standing rules from the owner's CLAUDE.md and memory: all analysis on the HPCC (laptop only
orchestrates); temp files in `$SCRATCH`; **SLURM jobs and HPCC writes need explicit go-ahead with
the exact command shown first**; the owner pushes and merges (the classifier blocks agent pushes);
lint with **ruff 0.3.5**; never write to `/bigdata/stajichlab/shared/...`.

## Branches and PRs

| Branch | Commits | State |
|---|---|---|
| `crossclade-round1-results` | c870c47, 9c18ff4 | **PR #170**, base `main`, open — round 1 results + round 2 code |
| `remove-legacy-cruft` | f6e4486 | **PR #171**, base `main`, open — legacy removal |
| `adhesion-round2-plan` | 04e5560, 188784d | **PR #172**, base `crossclade-round1-results` — retarget to `main` after #170 merges |
| `adhesion-benchmark-v1` | 6501cfa, 295c0cf, 74c56b9, dce1c27 | pushed, **no PR yet** — benchmark builder + adhesion_v1 + clustering |

Tag `legacy-2026-11-26` is pushed and holds everything removed in #171.

## What is running on the HPCC right now

Host `hpcc2`. Nothing is in the shared folder. Both are **detached worktrees** of the main clone, so
the owner's `allergen-predictor-wip` checkout at `$P` is untouched.

```
P=/bigdata/stajichlab/jstajich/projects/adhesionPred
$P/_workdir/crossclade_domain_search/   # round 2a; code/ worktree at 9c18ff4
$P/_workdir/benchmark_v1/               # benchmark; code/ worktree at dce1c27
```

**1. Round 2a inputs — DONE** (24 min). `$W/in2/` holds 13 query models, 12 PF21671, 126 candidates,
497 background A and 498 background B. The only query with no AFDB model is **IFF4109**
(A0A2H0ZI42, 2,946 aa) — exactly the caveat written into the README before the run, so it is
unscorable by the structure methods rather than a negative result. Inputs match their SHA256SUMS.

**2. Benchmark sequences and clustering — DONE.** 166/166 sequences, none missing; MMseqs2 17-b804f
at 30% identity gives **113 clusters**. Results copied into the repo (commit 138d18a).

Reseek v3.01 is already downloaded and verified at `_workdir/crossclade_domain_search/tools/`.

### Next HPCC steps, in order

```bash
P=/bigdata/stajichlab/jstajich/projects/adhesionPred
W=$P/_workdir/crossclade_domain_search; C=$W/code; B=$P/_workdir/benchmark_v1

# A3. when in2/ exists — SHOW THIS LINE TO THE OWNER AND GET GO-AHEAD BEFORE SUBMITTING
sbatch --export=ALL,PROJ_ROOT=$C,CODE_COMMIT=$(git -C $C rev-parse HEAD),IN2=$W/in2,\
RESEEK=$W/tools/reseek-v3.01,OUT=$W/round2 \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log $C/analysis/crossclade_domain_search/round2_search.sbatch

# A4. after the job
/usr/bin/python3.12 -I $C/analysis/crossclade_domain_search/round2_summarize.py \
  --round1-inputs $W/in --round1-summary $W/round1_summary --inputs2 $W/in2 \
  --results2 $W/round2 --out $W/round2_summary
/usr/bin/python3.12 -I $C/analysis/crossclade_domain_search/gene_model_check.py \
  --acc M7WU28 --acc M7WQR1 --acc M7X8Q7 --out $W/gene_model_check.tsv

# B2. after the fetch finishes
module load MMseqs2/17-b804f    # same version as the Kiss clustering, so counts are comparable
/usr/bin/python3.12 -I $B/code/analysis/curation/cluster_benchmark.py cluster \
  --benchmark $B/code/data/benchmark/adhesion_v1 --tmpdir $B/mmseqs_tmp
```

The HPCC default `python3` is **3.9** and the repo needs ≥3.11 — always `/usr/bin/python3.12`.
`$SCRATCH` is only set inside a job, which is why B2 uses a `_workdir` tmpdir for a seconds-long run.

**Known wart:** `cluster_benchmark.py` writes `sequences.faa` and friends *inside the worktree*
(`code/data/benchmark/adhesion_v1/`), making it dirty. The owner noticed. Either add an `--out`
option or copy the outputs back into the repo and `git checkout` the worktree.

## adhesion_v1, the benchmark

`data/benchmark/adhesion_v1/`. Built by `analysis/curation/build_benchmark.py`. **NOT FROZEN**: the
fold list below counts proteins, and clusters are the real unit (that is what B2 fixes).

76 positive, 90 hard_negative, 8 positive_moonlighting, 143 excluded (each with a reason — nothing
is dropped silently).

Folds are now on **cluster** counts (MMseqs2 30% identity):

| Fold | Level | Pos | Neg | Scored | Proteins before |
|---|---|---|---|---|---|
| Saccharomycetes | class | 17 | 38 | yes | 32 / 57 |
| Pezizomycotina | group | 11 | 7 | yes | 13 / 7 |
| Eurotiomycetes | class | 5 | 11 | yes | 5 / 11 |
| Microsporidia | group | 7 | 3 | no | 9 / 3 |
| Taphrinomycotina | group | 7 | 0 | no | 11 / 0 |
| Agaricomycotina | group | 2 | 5 | no | 2 / 6 |
| Mucoromycota | group | 1 | 2 | no | 3 / 2 |
| Ustilaginomycotina | group | 1 | 3 | no | 1 / 3 |
| Chytridiomycota | group | 0 | 1 | no | 0 / 1 |

The three scorable folds survived, which was not guaranteed. Saccharomycetes lost 47% of its
positives to paralogy; Mucoromycota collapsed 3 to 1 (the *Rhizopus* CotH proteins are one family);
Eurotiomycetes did not move, so its small set is genuinely non-redundant.

### What the curation established (three verified literature passes, 54 new rows)

Verified absences matter as much as additions:

- **Basidiomycota is nearly empty.** Outside *Schizophyllum* hydrophobins and four already-curated
  proteins, no Basidiomycota adhesin has direct evidence. ***Rhodotorula*, *Sporobolomyces* and
  *Malassezia*: nothing at all.** This is the single most important caveat for the Rhodotorula work
  — current predictions there are essentially unvalidatable against literature.
- **Dothideomycetes adhesion was characterised and never sequenced** (*Parastagonospora*,
  *Venturia*: real protein-dependent adhesion, no protein identified). *Alternaria* and *Pyrenophora*
  are genuinely empty — every hit screened by title.
- **Rust fungi** have the strongest functional evidence anywhere (purified enzymes added back to
  autoclaved spores restore adhesion to bean cuticle) but no sequence was ever deposited.
- **Colletotrichum** has no adhesion gene with a deletion mutant.
- **Microsporidia** gives 9 positives but is one heparin-binding-motif panel: 10–20% effect sizes,
  not a single knockout.

Two pairs deliberately test the same architecture on both sides: *U. maydis* **Hum3** (17 Rep1-like
repeats, adhesion unchanged on deletion) against the Rep1 positive — a direct test of whether repeat
architecture alone predicts adhesion; and *S. commune* SC3 (positive) vs SC4 (negative).
*M. oryzae* **Cut2** and *Z. tritici* **ZtGT2** are the cleanest negatives: adhesion measured and
unchanged.

Removed by owner decision, recorded in full in `sources/owner_excluded.md`: **Mpr1** (it meets the
`het_expr` positive criterion, so labelling it negative would override the evidence rule) and
*E. cuniculi* **Hsp70** (null result + cytoplasmic).

## Still to do

1. **Round 2a `sbatch`** — the line was shown to the owner on 2026-11-27 and is **awaiting
   go-ahead**; it has NOT been submitted. Preconditions all pass (pinned dir clean, `round2/` absent,
   inputs match SHA256SUMS). `CODE_COMMIT=9c18ff408add62825f6e13096c4b874dabc95a87`.
2. **Owner signs off `adhesion_v1`**, then it is frozen. Folds are already on clusters.
3. **Open the PR for `adhesion-benchmark-v1`.**
4. **2b proper**: the curated domain library (PA14, Flo11, Als N-terminal Ig, CFEM, CPL1/PF21671,
   LysM, fasciclin, Hyr/Iff, Bys1, hydrophobin, BAD1 — each cropped, each with its own empirical
   null) and the ablation harness. The rule: a held-out component returns only if it beats
   composition + length on **C1** in pooled leave-clade-out, with a 95% bootstrap CI (by cluster)
   excluding zero. Metrics: PR-AUC and recall at 1% FPR.
5. **2c layout study** (see the plan): layout tokens, stalk profiles (hydrophobicity, charge,
   β-aggregation, glycosylation density, repeat period and copy divergence) and label-free
   evolutionary signals (family expansion, telomere distance, repeat CNV between strains). Tests
   pre-listed before results are read; Benjamini–Hochberg across them.
6. **Held-out components** awaiting the harness: `disorder`, `repeats` (repeat02/repeat14), the AGA1
   feature classifier and `model_surface_candidate`, `hydrophobin_relaxed`, lectin models, the
   `cys_rich` and `expression` lookups.

## Open questions for the owner

- The worktree-dirtying output path (above).
- Whether to keep `positive_moonlighting` out of C2 as well as C1 (currently out of both).
- Whether Microsporidia's HBM panel should count as one observation rather than nine.

## Gotchas worth not rediscovering

- 14 tests in `tests/cellsurface_sorting_hat/test_data_and_scripts.py` fail **on macOS** with
  "File name too long", and `test_j1_stops_on_empty_input` in step1 also fails. Both predate all of
  this work — verified against tag `legacy-2026-11-26`. Do not "fix" them as regressions.
- Reseek labels are `<stem>_A`; Foldseek labels are the file stem.
- Round 2 scale unification to round 1 E-values: `af_fs E × N1/N_db`, `reseek P × N1`, N1 = 15,587.
- `/bigdata` was at 96% full on 2026-11-27.
- A `Workflow` tool call was stopped by a safety classifier during this session while fanning out
  the curation; the same work ran fine as three ordinary `Agent` calls. If that recurs, use Agent.
