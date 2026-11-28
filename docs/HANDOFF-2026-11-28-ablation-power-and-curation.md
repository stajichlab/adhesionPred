# Handoff 2026-11-28 — ablation power calibration, and what it says about curation

Continues `docs/HANDOFF-2026-11-27-adhesion-benchmark-and-round2.md`. Read that one first for the
benchmark's construction and the round 2 structural search; this one covers only what came after.

---

## State in one paragraph

`adhesion_v1` is frozen. The cross-clade structural search is finished and returned a negative. The
2b ablation harness was **not** built, because a power calibration run first (PR 0, merged as #179)
showed the benchmark cannot adjudicate the components the harness was going to test. The
calibration results and the curation analysis that follows from them are committed on
`ablation-power-calibration` and **not yet pushed**. There is nothing running on the HPCC.

## Where the code is

| | |
|---|---|
| `origin/main` | `9796b83` (PRs #175–#179 all merged) |
| current branch | `ablation-power-calibration`, **2 commits ahead of main, unpushed** |
| | `287e172` calibration report |
| | `25997d0` PF21671 correction + curation requirement |
| open PRs | none |

**First action on resume: push `ablation-power-calibration` and open the PR.** The owner pushes;
the agent does not.

Untracked and unrelated to this work, left alone: `m`, `to_import/kiss2019_jgi/`, `uv.lock`.

---

## The result

**`adhesion_v1` cannot adjudicate the components the 2b harness was designed to test.**

Job 32519688, `exfab`/gpu12, 16 cores, 1:15:15, commit `5bbd85a`, 1000 replicates and 1000
permutations, both declared training pools. Report:
`docs/reports/data/ablation_calibration_2026-11-27/`.

A domain call needs **66% sensitivity** on positives (at 5% FPR on hard negatives) before the
benchmark detects it at 80% power — and that is the *optimistic* training pool. No candidate
component has a measured sensitivity, so that is the standard each must clear sight-unseen.

Effective sample size is **16.4 positives**, not 50: 33 positive clusters, 27 of them singletons,
with ALS alone contributing 9 proteins.

### Three findings that were not in the plan

**1. The training pool was an undeclared assumption and it dominates.** It decides which bar even
wins — `composition_length` under `all`, bare `length` under `scored` — and moves the MDE by nearly
a factor of two. Under `scored` the largest fold trains on 32 proteins and a single noise column
costs 0.145 concordance. That figure has the same signature as the C-pinning bug fixed during
development, so it was checked: the per-fold C is recorded in `baseline.json` for both pools and
matches. It is genuine small-sample behaviour. It only moves the answer one way — `all` is the
optimistic pool.

**2. The decision rule is conservative, not anti-conservative.** The review feared a nominal 5%
might really be 9.5%. Measured: 0.006 at one noise dimension, 0.003 at five, 0.000 at twenty.
Adding features costs accuracy, so a noise component starts below zero while the permutation null
is symmetric about zero. Validity is safe; **power is being given away**. Centring the null on the
measured penalty would recover it — but that changes the decision rule and must be declared before
any component is tested.

**3. A wall that curation cannot move.** See `curation_requirement.md`. Fitting Δ(s) = 0.1064·s −
0.0217 (which reproduces the measured MDE: 0.655 predicted vs 0.66 observed), the intercept is the
measured dimensionality penalty, and Δ(s) = 0 at **s ≈ 0.20**. Below ~20% sensitivity a component
is not hard to detect — it is *not an improvement at all*. No quantity of curation changes that.

### What curation would actually buy

Extrapolated as 1/√N in independent positive families, **validated only at the current size**:

| to detect sensitivity | families needed | vs today's 33 |
|---|---|---|
| 0.50 | 77 | 2.3× |
| 0.40 | 175 | 5.3× |
| 0.30 | 728 | 22× |
| 0.25 | 3,172 | 96× |

Doubling the curated truth moves the detectable sensitivity from 0.66 only to 0.52.

---

## Correction carried forward

An earlier claim in this session — repeated in commit `287e172` — said PF21671 "recovers 12 of 50
positives = 24% sensitivity". **That is wrong.** PF21671 is a Pucciniomycotina (*Rhodotorula*)
family; the benchmark's three scorable folds are entirely Ascomycota and contain no *Rhodotorula*
row at all. Its sensitivity on these 50 positives is 0. Corrected in `25997d0`; the calibration
numbers themselves were never affected.

The corrected reading is sharper and should drive planning: **the project can discover in clades it
cannot evaluate in.** That taxonomic gap is independent of statistical power and is not fixed by it.

---

## Recommended direction (not yet an owner decision)

1. **Curate for coverage, not for power.** Eurotiomycetes has 5 positive clusters and cannot support
   a real held-out fold. Pucciniomycotina has no truth at all, which is why the structural result
   cannot be assessed. Both are fixable at realistic scale. Making a 25%-sensitivity component
   testable is not.
2. **Reframe the component question.** "Does PF21671 improve a global classifier?" needs ~1,000
   families. "Of the 12 proteins PF21671 flags, how many are genuine adhesins, and does it miss a
   known one in its own clade?" needs about a dozen curated entries and is answerable now.
   Component-ablation is an instrument for general-purpose features; most of what this project has
   found is clade-restricted.
3. **Settle the null-centring question before any component is tested**, so the choice cannot be
   made after seeing a result.

## Open measurements, in priority order

| | what | cost |
|---|---|---|
| 1 | **Learning curve.** Replace the 1/√N assumption: re-run the calibration on nested subsamples of the positive clusters (12/18/24/33, several draws each), fit the real exponent. This is the single most load-bearing assumption in the curation table. | one SLURM job, ~1 h |
| 2 | **`--attenuation 0` rerun.** `plant_heterogeneous` keeps the signal's sign outside the strong fold and the metric is rank-based, so that row is a second optimistic bound, not the realistic case. The genuinely non-transferring case is `--attenuation 0`. Flag already exists; no code change. | one pool, ~47 min |
| 3 | Deferred review items, before any PR 3: length-stratified reporting, label-associated missingness, reuse of `analysis/step1_compare/phasec/bootstrap.py`, pre-registration artefact, `indirect=yes` stratified report. | — |

---

## Running it again

```bash
# worktree is already detached at 5bbd85a
W=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/ablation_power
sbatch --export=ALL,PROJ_ROOT=$W/code,CODE_COMMIT=<sha>,BENCH=$W/code/data/benchmark/adhesion_v1,OUT=$W/<NEW> \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/ablation_power/power.sbatch
```

`$OUT` must not exist; the script refuses to overwrite. It also refuses a dirty or unpinned commit
and a benchmark that fails its `SHA256SUMS`.

## Gotchas added this session

- **Do not run the power code with `python3 -I`.** Isolated mode implies `-s`, which drops user
  site-packages, and numpy/scipy/sklearn live in `~/.local/lib/python3.12` on this cluster. The
  earlier failure was masked because the ssh wrapper swallowed the exit code. A capability check
  must use *the flag the run will use*; mine did not, which is how this got through.
- **`ssh hpcc2 'bash -s'` is not a login shell** — `sbatch`, `squeue` and `module` are missing from
  `PATH`. Use `bash -ls`.
- **`pkill -f <pattern>` matches its own ssh command line** and will kill the connection. Bracket a
  character (`run[_]power`) or use the job id.
- `scp` of a remote path containing `*` is glob-expanded by zsh before it leaves the laptop. Quote
  it, or use `/.`.
- The bar and the component-augmented model **must be fitted at the same per-fold C**. Fitting the
  bar with a tuned C while components ran at C=1.0 made pure noise appear to cost 0.15 concordance.
  `choose_bar` returns its per-fold C for this reason and two tests pin it.
- The power run writes nothing until a pool finishes, so a wall-clock timeout loses everything
  computed. Worth a per-grid-point checkpoint if the grid grows.

## Leftovers

- `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/ablation_power/` — 614 MB: `code/`
  worktree at `5bbd85a`, `logs/`, and `calibration/`. **Purgeable**: `calibration/` is now copied
  into the repo in full.
- Older loose ends still open: `input_tests` collects zero tests; `aga1_class` needs pandas; the
  analysis CI jobs are still `continue-on-error`, so a green tick there does not prove much.
- Contamination is systemic and undeclared per component:
  `data/controls/repeat-mechanism/curation_table.tsv` overlaps 105 of the 125 scored proteins and
  underpins the `repeats` component. AGA1 is worse — 66 of 166 (40%) of the benchmark was in its
  training set, and the leakage-free remainder (27 pos / 22 neg, Eurotiomycetes at 4/1) is too thin
  to test. Any PR 3 needs per-component provenance declarations before it means anything.
