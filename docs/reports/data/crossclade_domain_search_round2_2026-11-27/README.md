# Cross-clade domain search, round 2 results (2026-11-27)

Design, scales, cut-off and the 2-of-3 rule were declared before the run in the docstring of
`analysis/crossclade_domain_search/round2_summarize.py` and in the round 2 plan
(`docs/superpowers/plans/2026-11-26-adhesion-round2-benchmark-cruft-and-layout-study.md`). This note reports against them.
Nothing below changes a threshold after seeing the data.

## Run

| Step | Where | Record |
|---|---|---|
| inputs | HPCC login node, `round2_inputs.py` at 9c18ff4; AlphaFold DB models trimmed to pLDDT ≥ 50 | `search_provenance.txt` (input `SHA256SUMS` digest) |
| search | SLURM job 32510688, `exfab` (gpu12), **31 s**, 16 cores; Reseek v3.01, Foldseek 10-941cd33 (CPU) | `search_provenance.txt` |
| summary | `round2_summarize.py` at 3f4c5a4 with `/usr/bin/python3.12` | `summary.tsv`, `targets.tsv`, `pairs.tsv`, `nulls.tsv`, `background_hits.tsv`, `scales.txt` |
| gene models | `gene_model_check.py` at 3f4c5a4 (UniProt REST) | `gene_model_check.tsv` |

13 query domains (round 1's 14 minus *C. auris* IFF4109, 2,946 aa, which has no AlphaFold DB model), 1,138 target models,
756 scored query–target pairs. Raw search tables stay on the HPCC in `_workdir/crossclade_domain_search/round2/`.

**Scales** (`scales.txt`): all three methods are put on the round 1 E-value, the expected chance hits in a search of the
N1 = 15,587 round 1 target proteins. `af_fs` E × 15,587/1,138 = E × 13.7; `reseek` P × 15,587. Cut-off 1e-3, and a method
only supports a pair if it also beats that query's own null.

## Result

**12 targets confirmed, and they are exactly the 12 known PF21671 proteins of *R. toruloides* NP11.** No new protein was
confirmed in either panel.

| Panel | Method | Targets supported | NP11 PF21671 | Other NP11 | Other H99 | Background half B |
|---|---|---|---|---|---|---|
| basidio | prostt5 | 13 | 12/12 | 1 | 0 | 0/500 |
| basidio | af_fs | 12 | 12/12 | 0 | 0 | 0/500 |
| basidio | reseek | 1 | 1/12 | 0 | 0 | 0/500 |
| basidio | **2-of-3** | **12** | **12/12** | **0** | **0** | **0/500** |
| asco | prostt5 | 13 | 0/12 | 5 | 8 | 0/500 |
| asco | af_fs | 0 | 0/12 | 0 | 0 | 0/500 |
| asco | reseek | 0 | 0/12 | 0 | 0 | 0/500 |
| asco | **2-of-3** | **0** | **0/12** | **0** | **0** | **0/500** |

`background_hits.tsv` is empty: not one of the 500 half-B proteins passed any method, so the measured false-pass rate is
0/500 per method and 0/500 for the rule. That is an upper bound of about 0.6% at 95% confidence, not a demonstration of zero.

## What this overturns

**Round 1's new candidates do not survive real structures.** Round 1 proposed eight NP11 proteins beyond the Pfam reference
(M7WU28, M7WYW5, M7XG46, M7WND4, M7WT56, M7WHW1, M7XWK5, M7WS82). With AlphaFold DB models in place of ProstT5-predicted 3Di,
**none of them reaches 2-of-3 support on any of its 6–8 query pairs.** They fail by wide margins, not marginally: scaled
`af_fs` 0.067 to 13.3 and `reseek` 10.7 to 560, against a cut-off of 1e-3.

**The whole `asco` discovery panel collapses.** ProstT5 supported 13 targets there; real structures support none. Those 13
were the Hyr/Iff, Flo11 and PA14 queries hitting NP11 and H99 proteins at E 1e-4 to 1e-3 — exactly the regime round 1 flagged
as having an unreliable null. `nulls.tsv` shows why: for the `asco` queries the background best `af_fs` score is 0.13 to 3.94,
so random proteins score about as well as the candidates, whereas for the PF21671 queries the background best is 23 to 98.

The honest reading is that **ProstT5-predicted 3Di generated hits that real structures do not support.** Round 1's clean
decoy null did not catch this, because composition-matched sequence decoys test a different failure than real structures do.
That is the main methodological result of round 2.

## The 2-of-3 rule is carried by two non-independent methods

**All 12 confirmations rest on `prostt5` + `af_fs`, which both score with Foldseek.** Reseek — the one genuinely independent
method — supports a single pair in the whole run (M7WVC0 × J9VVN7, scaled P 2.2e-4), so 11 of the 12 are confirmed by two
methods that share a scoring engine. The rule's stated caveat about non-independence is therefore not hypothetical; it is
doing the work in this run.

Reseek nonetheless **separates the two classes cleanly**, it just does not cross the pre-declared line: scaled P for the 12
PF21671 targets is 2.2e-4 to 0.34 (10 of 12 below 0.1), versus 10.7 to 1,540 for everything else — three to four orders of
magnitude. Lowering the cut-off to 0.1 would confirm 10 of 12 by Reseek alone. **This threshold is not being changed**; it was
declared before the run and the number is reported as it fell. It is recorded here as an observation for the next design.

## Gene-model check (`gene_model_check.tsv`)

Three accessions were checked for a truncated N terminus. **None of the three is a gene-model artefact.**

| Accession | NP11 | UniRef50 relatives | Verdict |
|---|---|---|---|
| M7WU28 | 209 aa, no SP, RHTO_01851 | **none** — a singleton cluster (`UniRef50_M7WU28`) | no relatives to compare; cannot be rescued or refuted |
| M7WQR1 | 985 aa, no SP | 11 *Rhodotorula*, 951–1,081 aa, **all with no signal peptide** | the missing SP is conserved across 11 species: real, not an annotation error |
| M7X8Q7 | 574 aa, **SP 1–20** | 6 *Rhodotorula*, 569–580 aa, several with an SP | already secreted; model is fine |

**M7WU28 is the one genuinely unresolved protein in the run.** It is the single basidio target outside PF21671 that ProstT5
supported, and it did so from **6 independent queries** at E 1.9e-5 to 8.4e-4. But its AlphaFold DB model is `too_short`
(mean pLDDT 36.6), so neither structure method could score it, and its UniRef50 cluster has no other member. It is neither
confirmed nor refuted, and no further sequence or structure evidence is currently available for it.

## Reseek v3.01 input bug (affects anyone re-running this)

Reseek's documentation lists a directory as a valid `STRUCTS` argument for both `-search` and `-db`. **For `-search`, neither
the directory nor the documented `.files` list works in v3.01.**

- directory → `ReadStdioFile failed, attempted 4 bytes, read 0 bytes, errno=21` (EISDIR). This killed job 32502912.
- `.files` list → `Bad magic 6769622f, invalid .bcx file` — it parses the list as a binary structure file (`6769622f` is
  ASCII `gib/`, the first bytes of a `/bigdata/...` path).
- `.bcb` (converted with `-convert`) → **works.** This is the form the run used, recorded in `search_provenance.txt`.

`-convert` accepts a directory happily, so the inconsistency is within one binary. The job script now tries all three forms,
requires a non-empty result (Reseek exits 0 on some failures), and records which form won. The binary needs CPU features the
login node lacks, so no query form can be probed outside a SLURM job.

## Files

| File | Content |
|---|---|
| `summary.tsv` | per panel and method: targets supported, PF21671 recovery, background half-B pass rate |
| `targets.tsv` | one row per target: panel, model status, pLDDT, SP/GPI flags, best score per method, confirmed yes/no |
| `pairs.tsv` | all 756 query–target pairs with each method's scaled score and support call |
| `nulls.tsv` | each query's null per method (the bar a hit must beat) |
| `background_hits.tsv` | half-B proteins that passed — **empty** |
| `scales.txt` | N1, N_db and the scale factor per method |
| `gene_model_check.tsv` | UniRef50 relatives, lengths and signal peptides for the three checked accessions |
| `search_provenance.txt` | job provenance: commit, node, tool versions, Reseek query form, row counts |

## What this means for the project

The *Rhodotorula* CPL1-like family is solid: 12 proteins, recovered by every round and by two methods with a zero background
pass rate. **But round 2 found no new adhesion-relevant protein, and it removed eight that round 1 had proposed.** The
cross-clade structure search, run honestly against real structures and a real null, returns the family that was already known
from Pfam and nothing else.

That is consistent with the project's standing finding — that the signal available to these methods is weaker than it looks
once the null is matched properly — and it is the reason the 2b benchmark work (composition-matched hard negatives, clades
held out) is the right next step rather than more search rounds.
