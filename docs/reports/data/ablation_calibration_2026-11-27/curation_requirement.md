# How much curation growth would the ablation design need?

Companion to `README.md`. The calibration says a component needs **66% sensitivity** before
`adhesion_v1` can detect it. The obvious follow-up is "so how much bigger must the benchmark be?"
This note answers that, and the answer is not "somewhat bigger".

**Status: extrapolation, not measurement.** Everything below assumes the detectable effect shrinks
as 1/√N in the number of independent positive families. That is the standard asymptotic rate and it
reproduces the measured MDE at current size (0.655 predicted vs 0.66 measured), but it has not been
verified at other sizes. The measurement that would settle it is a learning curve — subsample the
benchmark to several sizes and re-measure the MDE at each. That job has not been run.

## What the benchmark is made of

| | value |
|---|---|
| positives (scorable folds) | 50 proteins in **33 clusters** |
| Kish effective sample size | **16.4** |
| singleton clusters | 27 of 33 |
| largest cluster | **ALS, 9 proteins** (18% of all positives) |
| positive clusters per fold | Saccharomycetes 17 · Pezizomycotina 11 · Eurotiomycetes 5 |

The clustering is already doing its job — 27 of 33 families are singletons. The problem is not
redundancy inflation. It is **absolute scarcity of independent families**.

## The scaling

Fitting the measured relationship between a domain call's sensitivity and the concordance gain it
produces (over s ≥ 0.25, pool `all`):

> Δ(s) = 0.1064·s − 0.0217

| independent families | growth | effective N | detectable Δ | **detectable sensitivity** |
|---|---|---|---|---|
| 33 | 1× | 16 | 0.0480 | **0.66** |
| 50 | 1.5× | 25 | 0.0392 | 0.57 |
| 66 | 2× | 33 | 0.0339 | 0.52 |
| 99 | 3× | 49 | 0.0277 | 0.46 |
| 132 | 4× | 66 | 0.0240 | 0.43 |
| 330 | 10× | 164 | 0.0152 | 0.35 |
| 1386 | 42× | 689 | 0.0074 | 0.27 |

Read the other way:

| to detect a component at sensitivity | needs growth | independent positive families |
|---|---|---|
| 0.50 | 2.3× | **77** |
| 0.40 | 5.3× | **175** |
| 0.30 | 22× | **728** |
| 0.25 | 96× | **3,172** |

Doubling the curated truth — 33 families to 66, a very large undertaking — moves the detectable
sensitivity only from 0.66 to 0.52.

## There is a wall that curation cannot move

The fitted intercept is negative, and that is not noise: it is the **dimensionality penalty**
measured directly in `type_i.tsv` (a pure-noise component costs −0.007 at one dimension, −0.058 at
twenty). A component must first pay for the features it adds before any of its signal counts.

Setting Δ(s) = 0 gives **s ≈ 0.20**.

Below roughly 20% sensitivity, a domain call is not "hard to detect" — it is **not an improvement
at all** under this bar and this rule. It makes the model worse. No quantity of curation changes
that, because the limit is not sampling error. Only reducing the penalty can: fewer features per
component, stronger regularisation, or centring the null on the measured penalty instead of on zero.

This is why the 0.25 row above reads 3,172 families. The design is asking for an impossible number
because it is approaching an asymptote, not because the truth set is slightly too small.

## What this means for the project

**1. Curation is still necessary, but not sufficient, and not primarily for the ablation harness.**
At 5 positive clusters, Eurotiomycetes cannot support a meaningful held-out fold at all; that is a
coverage failure curable by curation at a realistic scale. Making a 24%-sensitivity component
testable is not.

**2. There is a taxonomic gap between where discovery happens and where evaluation is possible.**
The cross-clade structural search confirmed PF21671 in *Rhodotorula* (Pucciniomycotina). The
benchmark's three scorable folds are entirely Ascomycota and contain no *Rhodotorula* at all. The
project can currently *discover* in clades it cannot *evaluate* in. Curation aimed at the ablation
harness and curation aimed at closing that gap are different target lists, and the second is the
one that makes the structural work assessable.

**3. The ablation question is the wrong question for a narrow component.** "Does adding PF21671
improve a global adhesin classifier?" needs on the order of a thousand independent families to
answer. "Of the proteins PF21671 flags, how many are genuine adhesins, and does it miss any known
one in its own clade?" needs only the twelve proteins it flags, plus curated truth in
Pucciniomycotina. The second question is answerable now, at a curation cost of perhaps a dozen
entries, and it is the question a biologist would actually ask of a domain family.

A component-ablation harness is the right instrument for asking whether a *general-purpose* feature
earns its place in a *general-purpose* model. It is the wrong instrument for a clade-restricted
domain family, and most of what this project has found is clade-restricted.

## Recommended next measurement

Replace the 1/√N assumption with a learning curve before any curation target is committed to:
re-run the calibration on nested subsamples of the positive clusters (say 12, 18, 24, 33 families,
several draws each) and fit the observed MDE against family count. If the real exponent is shallower
than 1/√N the table above is pessimistic; if steeper, the wall is nearer than it looks. Either way
the number in a curation plan should be measured, not assumed.
