# Ablation power calibration — can `adhesion_v1` adjudicate anything?

**Date:** 2026-11-27 · **Job:** 32519688 (`exfab`, gpu12, 16 cores, 1:15:15) · **Code:** `5bbd85a`
**Benchmark:** `adhesion_v1`, SHA256SUMS digest `ae9a6995…3893` · **Settings:** 1000 replicates, 1000
permutations, seed 0, α = 0.05, both declared training pools.

This ran **before** the 2b ablation harness was built, to answer one question first:

> How large an improvement would a component have to make before this benchmark could tell?

**No real component was tested.** Every signal here is planted.

---

## The answer

**`adhesion_v1` cannot adjudicate the components the 2b harness was designed to test.**

A domain call must reach **66% sensitivity** on positives, at 5% false-positive rate on hard
negatives, before this benchmark detects it at 80% power — and that is the *optimistic* training
pool. No candidate component has a measured sensitivity on this benchmark, so that threshold is
the standard every one of them has to clear sight-unseen.

PF21671, the only domain family the cross-clade structural search confirmed, **cannot be evaluated
here at all**, for a reason more basic than power: it is a Pucciniomycotina (*Rhodotorula*) family,
the three scorable folds are entirely Ascomycota, and there is no *Rhodotorula* row anywhere in the
benchmark. Its sensitivity on these 50 positives is 0, not low.

Had PR 1 been built first, it would have produced six null results, not one of them informative.

### Minimum detectable effect at 80% power

| signal | pool `all` | pool `scored` |
|---|---|---|
| homogeneous — same effect in every clade (optimistic bound) | d = 1.48 (Δ +0.041) | d = 2.47 (Δ +0.081) |
| heterogeneous — see caveat 1 below | d = 1.62 (Δ +0.029) | not reached on the grid |
| **domain call at 5% FPR (sensitivity needed)** | **0.66** (Δ +0.048) | **not reached (> 0.80)** |

### Power against the actionable shape (pool `all`)

| sensitivity | 0.05 | 0.15 | 0.25 | 0.35 | 0.45 | 0.60 | 0.80 |
|---|---|---|---|---|---|---|---|
| power | 0.013 | 0.073 | **0.159** | 0.297 | 0.434 | 0.730 | 0.956 |

The bolded column is the level at which a component contributes almost nothing detectable: a call
that finds a quarter of the positives is seen one time in six.

---

## Why the benchmark is smaller than it looks

- 50 positives and 75 hard negatives, in 84 clusters across the scorable folds.
- **Effective sample size: 16.4 positives** (Kish ESS over 33 positive clusters). ALS alone is 9 of
  those 33 clusters. A paralogue family is one observation, not nine.

Quote the effective number, not 50.

---

## The training pool was an undeclared assumption, and it dominates

The original 2b plan never stated which proteins train the bar. It turns out to decide which bar
even wins, and to move the MDE by a factor approaching two:

| | pool `all` | pool `scored` |
|---|---|---|
| bar chosen by inner cluster-grouped CV | `composition_length` (21 features) | `length` (1 feature) |
| largest fold's training set | 44 pos / 29 neg = **73** | 18 pos / 14 neg = **32** |
| mean Δ under a 1-dim noise component | −0.007 | **−0.145** |
| type-I rate, 1 / 5 / 20 dims | 0.006 / 0.003 / 0.000 | 0.000 / 0.000 / 0.000 |

Under `scored`, the largest fold trains on 32 proteins, and the winning bar is a single feature.
Adding one noise column to a one-feature model at that sample size costs 0.145 concordance — twenty
times the penalty under `all`.

That figure has the same arithmetic signature as the C-pinning bug fixed during development (a noise
component appearing to cost 0.15), so it was checked rather than assumed: the bar and the augmented
model are fitted at the same per-fold C, which `baseline.json` records under `bar_C_per_fold` for
both pools. The effect is genuine small-sample behaviour.

**This only moves the answer one way.** `all` is the optimistic pool and even it demands 66%
sensitivity. The conclusion is robust to the assumption that was never declared.

---

## The decision rule is conservative, not anti-conservative

The review's concern was that a nominal 5% might really be 9.5%. The measurement says the opposite:

| component dims | rejection rate, pool `all` | mean Δ under the null |
|---|---|---|
| 1 | 0.006 | −0.0067 |
| 5 | 0.003 | −0.0206 |
| 20 | 0.000 | −0.0582 |

This is explainable. Adding features *costs* accuracy, so a pure-noise component starts below zero
and must climb back before it can even reach zero. The permutation null is symmetric about zero by
construction, while the true null sits below it — so the rule under-rejects, and under-rejects
harder as the component gets wider.

**Validity is safe; sensitivity is being given away.** The MDE figures above already carry this
penalty, so they are honest about what it takes to clear the rule *as written*. Centring the null on
the measured dimensionality penalty would recover real power at no cost to validity — but that is a
change to the decision rule, and it must be declared before any component is tested.

---

## Caveats

**1. The heterogeneous row does not measure what it claims.** `plant_heterogeneous` scales the
signal by `attenuation` (default 1/3) outside the strong fold but **keeps the same sign everywhere**.
The primary metric is rank-based, so a model trained on the attenuated folds still learns the right
direction and ranks the held-out clade nearly as well as under a uniform signal. That is why the
heterogeneous and homogeneous curves nearly coincide (d = 1.62 vs 1.48) when the docstring predicts
a large gap.

A genuinely non-transferring signal is `attenuation = 0`: present in one clade, absent elsewhere.
The flag already exists (`--attenuation`), so this needs a rerun, not a code change. It matters
because clade-restricted families are the normal case in this project — PF21671 is one, and the
expected result, power collapsing to the type-I rate at every effect size, would mean a
narrow-clade family is undetectable by leave-one-clade-out at **any** effect size. That is a
sharper and more general constraint than "needs d = 1.62".

**Until that rerun, treat the heterogeneous row as a second optimistic bound, not as the realistic
case.** The binary row is unaffected and remains the actionable one.

**2. Power estimates near 0.8 carry a standard error of about 1.3%** at 1000 replicates. The MDE is
read off a 9-point grid by interpolation, so quote it to two significant figures at most.

---

## How to read a null result from this benchmark

A component that fails to clear this bar has **not** been shown to be useless. It has been shown to
be **undetectable at this benchmark's size**. Report the interval's upper bound — "rules out an
improvement above X" — never the bare phrase "no improvement".

Treating "no evidence against" as "evidence for keeping" is how the cruft removed in PR #171
accumulated in the first place. The same error in the other direction would discard working
components on the strength of six uninformative nulls.

---

## What follows

1. **Curating more truth outranks building the ablation harness** — but see
   `docs/reports/data/ablation_calibration_2026-11-27/curation_requirement.md` for how much growth
   the ablation design would actually need, because the answer is not "somewhat more".
2. **Rerun the heterogeneous shape at `--attenuation 0`** before the heterogeneous row is cited.
3. **Decide the null-centring question before any component is tested**, so the choice cannot be
   made after seeing a result.

## Files

| file | what it holds |
|---|---|
| `provenance.txt` | commit, node, interpreter, package versions, benchmark digest, settings |
| `{all,scored}/baseline.json` | bar, per-fold C, splits and leakage drops, per-fold report, effective sample sizes, every setting |
| `{all,scored}/type_i.tsv` | rejection rate and mean null Δ at 1, 5, 20 noise dimensions |
| `{all,scored}/power.tsv` | power and mean Δ for each shape × effect |
| `{all,scored}/summary.md` | the run's own generated summary |

Reproduce with `analysis/ablation_power/power.sbatch` at commit `5bbd85a`. Note that the run must
**not** use `python3 -I`: isolated mode drops user site-packages, where numpy, scipy and
scikit-learn live on this cluster.
