# 2b, part 1: the ablation harness (plan, 2026-11-27)

> **SUPERSEDED IN PART by an independent review, 2026-11-27.** The review implemented this design
> and measured it rather than reasoning about it. Ten findings are blocking. Four were re-verified
> independently before being accepted:
>
> | Finding | Status |
> |---|---|
> | **log length ALONE scores AP 0.690; the specified 21-feature baseline scores 0.551** | verified — the bar is beaten by one of its own inputs |
> | **Zero mixed-label clusters inside the scorable folds** | verified — P2's stated justification below is **wrong** (D-B is still right, for the paralogue reason) |
> | **The D-A table only reproduces if the training pool is all 166 proteins** | verified — Pezizomycotina drops 5 vs 2; the plan never says which pool it means |
> | **Length flips sign between folds** (Sacch length-only AUC 0.912, Eurotiomycetes 0.327) | verified — pooling scores from models that learned opposite rules is incoherent, not merely uncalibrated |
>
> **The headline is that this plan's own PR-2 gate is close to tripped before PR 1 is written.**
> Measured minimum detectable effect at 80% power: ΔAP ≈ 0.085 for a clade-homogeneous planted
> signal, ≈ 0.13 for one that transfers only partly across clades — and that is from an
> *uncorrected* test measured to be anti-conservative (true type-I 9.5% at nominal 5%). Correct the
> type-I rate and apply multiplicity across six components and the per-component MDE exceeds the
> plan's own 0.15 abort threshold. In usable terms: **a domain model would have to hit ~45% of
> adhesins at 5% FPR against hard negatives before this benchmark could see it.**
>
> See `docs/reports/` (pending) for the review's measurements. The design sections below are kept
> for the record; do not implement them as written.

Builds the apparatus that decides whether any component of this project earns its place. The domain
library (2b part 2) comes after, because the harness is what makes a claim about the library
checkable.

Parent plan: `docs/superpowers/plans/2026-11-26-adhesion-round2-benchmark-cruft-and-layout-study.md`
(owner decisions D3, D4 of 2026-11-26). Benchmark: `data/benchmark/adhesion_v1/`, frozen 2026-11-27.

## Why

Every model this project has trained reduces to amino-acid composition plus length: ESM-2 0.999 vs
0.995 leave-family-out, ESM C 0.598 leave-clade-out, the AGA1 feature classifier 0.97 vs 0.973
composition-only. Round 2a then found nothing beyond what Pfam already knew. Components are
currently held out on suspicion rather than on measurement, and there is no instrument that could
let one back in.

This plan builds that instrument, and — before testing anything with it — establishes what it is
capable of detecting.

## What already exists

`data/benchmark/adhesion_v1/`, frozen: 76 positives, 90 hard negatives, 8 moonlighting, 143 excluded
with reasons; 166 sequences in 113 MMseqs2 clusters at 30% identity; `folds.tsv` naming three
scorable folds. Nothing of the harness exists — no baseline implementation, no CV loop, no metrics.

## Three problems found in the pre-declared design

These were found by checking the frozen benchmark against the plan's rules, before writing code.
All three need a decision; none requires changing the benchmark.

### P1. Five clusters span folds, so leave-clade-out leaks as written

| Cluster | Folds it touches |
|---|---|
| P16933 | Agaricomycotina, Pezizomycotina, Ustilaginomycotina |
| Q59NP5 | Pezizomycotina, Saccharomycetes |
| P87020 | Eurotiomycetes, Saccharomycetes |
| P0C7S9 | Eurotiomycetes, Saccharomycetes |
| P38993 | Eurotiomycetes, Saccharomycetes |

Holding out Eurotiomycetes leaves three of its clusters represented in training by their
Saccharomycetes members, at ≥30% identity. That is precisely the leakage leave-clade-out exists to
prevent, and it inflates exactly the transfer claim this project is trying to test honestly.

**Fix, and it is cheap.** When fold F is held out, drop from *training* every protein whose cluster
contains any member of F. Measured cost:

| Held-out fold | Test pos / neg | Train pos / neg | Proteins dropped for leakage |
|---|---|---|---|
| Saccharomycetes | 32 / 57 | 44 / 29 | 4 |
| Pezizomycotina | 13 / 7 | 62 / 79 | 5 |
| Eurotiomycetes | 5 / 11 | 71 / 75 | 4 |

Four or five proteins per fold. **Decision D-A: adopt the strict rule.** The cost is negligible and
the alternative is a number we could not defend.

### P2. One cluster contains both a positive and a hard negative, so clusters cannot be the scoring unit

> **CORRECTED 2026-11-27.** The justification below is **wrong**, verified independently: there are
> **zero** mixed-label clusters inside the three scorable folds. P16933's members sit in
> Agaricomycotina and Ustilaginomycotina (both `scored=no`) plus one Pezizomycotina positive, so the
> SC3/SC4 and Rep1/Hum3 showcase pairs are **entirely outside C1 and decide nothing**. D-B is still
> the right call, but the real reason is the FLO/ALS/EPA paralogue families, which would otherwise
> count as independent evidence in the interval. The "pairs that make the benchmark hard" framing
> must not be carried into any harness report.

P16933 holds both tiers. This is deliberate and documented: the benchmark keeps *S. commune* SC3
(positive) against SC4 (negative), and *U. maydis* Rep1 (positive) against Hum3 (negative), because
they test the same architecture on both sides. They are hydrophobins and cluster together.

So a cluster has no single label, and collapsing to one row per cluster is impossible without
destroying the pairs that make the benchmark hard.

**Decision D-B: score at the protein level, resample by cluster.** Metrics are computed over
proteins; the bootstrap resamples *clusters* and carries all their members. This keeps the
SC3/SC4-style pairs intact while preventing a paralogue family from counting as independent
evidence in the confidence interval.

### P3. "Recall at 1% FPR" is below the resolution of the data

Pooled over the three scorable folds the benchmark has **50 positives and 75 hard negatives**. The
smallest non-zero false-positive rate expressible with 75 negatives is 1/75 = **1.3%**. A 1% FPR
target cannot be estimated; any implementation would silently round to 0 false positives, which is
a different and much noisier quantity.

**Decision D-C needed from the owner.** Options, in my order of preference:

1. **Average precision (pooled out-of-fold) as the primary metric, plus recall at 5% FPR as
   secondary.** 5% FPR is ~4 negatives of 75 — coarse, but estimable, and its CI will be honest
   about that. *Recommended.*
2. Keep a fixed false-positive *count* instead of a rate: recall at 1 FP and at 3 FP. Interpretable
   on a small set, but not comparable to any published figure.
3. Partial AUC over FPR ∈ [0, 0.1], normalised. Statistically cleaner, harder to explain.

I recommend (1) and will write the plan against it unless told otherwise.

### P4. The AGA1 feature classifier was fitted on 40% of the benchmark, and cannot be scored on it

Checked 2026-11-27 against `docs/reports/data/aga1_feature_classifier_2026-10-25/features.tsv`. The
classifier was fitted on 596 proteins (role `class_positive` tier A, plus `hard_negative` and
`background_assumed_negative`). **66 of the benchmark's 166 scored proteins are among them**, with
the labels agreeing — benchmark `positive` ↔ AGA1 `class_positive/A`, benchmark `hard_negative` ↔
AGA1 `hard_negative`. That is direct label leakage, not incidental overlap.

Inside the three scorable folds, removing both the fitted proteins and any protein whose *cluster*
was fitted on (homology leakage at 30% identity, the same standard the splits use):

| Fold | pos / neg | pos fitted on | neg fitted on | leakage-free pos / neg |
|---|---|---|---|---|
| Saccharomycetes | 32 / 57 | 17 | 28 | 13 / 15 |
| Pezizomycotina | 13 / 7 | 2 | 0 | 10 / 6 |
| Eurotiomycetes | 5 / 11 | 1 | 10 | 4 / 1 |
| **pooled** | **50 / 75** | **20** | **38** | **27 / 22** |

Eurotiomycetes falls to 4 positives and 1 negative — not a fold. So there is no leakage-free way to
score this component on `adhesion_v1` that preserves the design.

**Decision D-D: the AGA1 feature classifier and `model_surface_candidate` (which consumes it) are
not evaluable on `adhesion_v1`.** They should be reported as *untestable on this benchmark* rather
than given a number. A number computed on 40% of its own training set would be the exact mistake
this whole exercise exists to stop making. Evaluating them needs a fresh set with no overlap, which
is a v2 question, not a harness question.

**The other held-out components are not affected in the same way**, checked the same day:

- `disorder` and `repeats` compute deterministic properties from sequence. Nothing is fitted, so
  there is no training set to leak. Clean.
- `hydrophobin_relaxed` uses the seven Pfam hydrophobin HMMs *unchanged*; only its bit-score cutoff
  was tuned, and that tuning was frozen with its own held-out folds and a written record of what
  the author had seen (`analysis/hydrophobin_truth/freeze.json`). One tuned scalar is a far smaller
  exposure than a fitted classifier, but its provenance file says "partial leakage" in as many
  words, so the overlap between its tuning set and the benchmark's hydrophobins must be checked
  before it is scored. Tracked as a task in PR 3, not a blocker for PR 1.
- The lectin models were not located in this pass and still need the same check.

## Design

### Data contract

The harness reads only the frozen benchmark and never writes to it:

- `proteins.tsv` — `accession`, `tier`, `clade_class`, `length`
- `clusters.tsv` — `accession` → `cluster_representative`
- `folds.tsv` — fold, members, `scored`
- `sequences.faa` — one sequence per scored protein

It refuses to run if `SHA256SUMS` does not verify, and records the benchmark `VERSION` and digest in
every output. A result that cannot name the benchmark it came from is not a result.

### Contrasts

- **C1, positives vs hard negatives** — the only contrast that decides keep or remove.
- **C2, positives vs other secreted/GPI proteins** — descriptive, assembled at scoring time, not
  part of any keep/remove decision. Out of scope for this plan; the harness leaves a hook for it.

### Splits

Leave-one-clade-out over the three scorable folds, with the strict cluster rule from D-A. Folds
marked `scored=no` are reported with their counts and never silently merged.

Note what the table above shows: holding out Saccharomycetes leaves 44 positives and **29** hard
negatives to train on. That fold is a genuinely different experiment from the other two, and its
result should be read separately as well as pooled.

### Baseline

Logistic regression on 20 amino-acid frequencies + log length, refit inside every split, with the
scaler fit on training data only. This is the bar. It is not a component and is never ablated.

### Components

Each component is an adapter exposing `features(accessions) -> array` plus a name and a version, so
the harness does not care how a component is computed. First adapters, from the held-out list:
`disorder`, `repeats` (repeat02/repeat14), `cys_rich`, the AGA1 feature classifier,
`hydrophobin_relaxed`, the lectin models.

A component that needs an external binary or a precomputed table declares it and is **skipped with a
reason** rather than silently scoring zero.

### Metrics and statistics

Out-of-fold predictions from all three folds are **pooled**, then scored once. Per-fold metrics are
also reported, but the pooled figure is primary: average precision on a fold with 5 positives is
noise, not a measurement.

*(Caveat to state in the output: pooling concatenates scores from three separately fitted models,
which are not calibrated to each other. Per-fold figures are reported alongside so this is visible.)*

The test is a **paired bootstrap by cluster**: resample clusters with replacement, recompute both
models' metrics on the same resample, take the difference. 10,000 resamples, percentile interval.
Paired, because the question is "does adding this help", not "are these two numbers different".

**The rule, unchanged from the parent plan:** a held-out component returns only if, on C1, it beats
the baseline (or baseline + already-kept components) in pooled leave-clade-out with a 95% CI
excluding zero.

### Step 0, before any component is tested: calibrate what the instrument can detect

With 50 positives and 75 negatives, the harness will often return "cannot tell". That outcome is
only interpretable if we know in advance how large a difference it could have detected.

So the first deliverable after the harness runs is a **null and power calibration**:

- a **permuted-label component** (shuffled within fold) — the CI must contain zero;
- a **random-noise component** of the same dimensionality as a real one — the CI must contain zero;
- a **planted-signal sweep**: inject a synthetic feature correlated with the label at varying
  strength, and report the smallest ΔAP whose CI excludes zero.

That last number is the headline of the first report: **the minimum detectable improvement.** Every
later "no improvement" is read against it. Without it, a null result is unreadable, and there is a
standing temptation to treat "no evidence against" as "evidence for keeping".

## Deliverables

One PR each, in order.

**PR 1 — harness core.** `analysis/ablation_harness/`: `data.py` (loads and verifies the frozen
benchmark, builds the strict splits), `features.py` (baseline features, component adapter
protocol), `metrics.py` (AP, recall at a given FPR, paired cluster bootstrap), `harness.py` (the CV
loop), `run_ablation.py` (CLI). Tests in `tests/ablation_harness/`, added to the no-torch CI job.

**PR 2 — calibration report.** Step 0 above, run on the HPCC, written to
`docs/reports/data/ablation_calibration_<date>/` with the minimum detectable improvement stated in
the README.

**PR 3 — first real components.** The cheap ones with no external dependency (`disorder`, `repeats`,
`cys_rich`), each reported against the bar, with the keep/remove/cannot-tell verdict.

## Tests

The statistics are the part that can be wrong while looking right, so they get tests with known
answers:

- splits: no test protein's cluster appears in training, for every fold (the D-A rule, asserted);
- splits: the dropped-for-leakage counts match the table above exactly;
- metrics: AP against a hand-computed small example; recall at FPR against a hand-computed example;
- bootstrap: a component identical to the baseline gives a CI containing zero; a perfectly
  predictive component gives a CI excluding zero;
- bootstrap: resampling is by cluster, verified by checking that members of one cluster always
  appear together;
- determinism: a fixed seed reproduces the interval exactly;
- the harness refuses to run on a benchmark whose `SHA256SUMS` does not verify.

## Order

1. Owner decides **D-C** (the metric). D-A and D-B I consider settled by the measurements above
   unless objected to.
2. PR 1, harness core, with tests.
3. PR 2, calibration — **this is the gate**. If the minimum detectable improvement turns out to be
   large (say ΔAP > 0.15), that is itself the finding, and it changes what is worth doing next:
   growing the benchmark would then outrank testing components.
4. PR 3, the cheap components.
5. Then 2b part 2, the domain library, each domain with its own empirical null.

## Open questions

- **C2 is specified but has no data yet.** It needs proteomes and a secreted/GPI call. Out of scope
  here; worth deciding whether it is wanted at all, since it decides nothing.
- **~~Does the AGA1 classifier's training data overlap the benchmark?~~ Checked 2026-11-27: yes,
  disqualifyingly.** See P4 below.
- **Moonlighting proteins** stay out of both contrasts, per the 2026-11-27 decision. Unchanged.
