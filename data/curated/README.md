# Curated label tables (DRAFT — all need expert review)

Five label categories, built by scripts in `analysis/curation/`. They share one population
(fungal cell-surface and secreted proteins) and are deliberately kept separate, because a
protein can be several things at once: an adhesin that is also a biofilm factor and an antigen.

| folder | table | what it labels | rows (2026-09-27) |
|---|---|---|---|
| `surface/` | `surface.tsv` | **stage 1**: surface glycoproteins, each marked adhesin / non-adhesin / unknown | 3,540 |
| `adhesins/` | `adhesins.tsv` | **stage 2**: adhesins and hard negatives, with evidence levels | 250 |
| `biofilm/` | `biofilm.tsv` | biofilm involvement, split into surface factors vs regulators | 207 |
| `antigens/` | `antigens.tsv`, `coccidioides_candidates.tsv` | immune-recognized proteins (IEDB) and ranked *Coccidioides* candidates | 86 + 1,069 |
| *(planned)* | allergens | fungal allergens (issue #19) | — |

## Why the split

The shipped classifier detects the **whole surface class**, not adhesins specifically
(~12% adhesin precision on S288C; see `docs/model-review/`). Stage 1 is therefore a
legitimate predictor in its own right and is already close to solved. The hard problem is
stage 2: telling adhesins apart from the other surface glycoproteins. That needs the
non-adhesive members as explicit negatives, which is what `surface.tsv` supplies.

## Current label counts for stage 2 (from `surface.tsv`)

| status | n | meaning |
|---|---|---|
| `adhesin` | 75 | curated positives |
| `non_adhesin` | 63 | curated negatives |
| `non_adhesin_putative` | 1,257 | proposed: annotated catalytic activity, no adhesion evidence |
| `unknown` | 2,145 | unlabeled — **not** negatives (use positive-unlabeled learning) |

## Rebuilding

```bash
python analysis/curation/build_adhesins.py           # QuickGO + literature + Pfam families
python analysis/curation/build_biofilm.py            # QuickGO biofilm terms
python analysis/curation/build_antigens.py           # IEDB immune-assay evidence
python analysis/curation/build_surface.py            # UniProt surface keywords + cross-links
python analysis/curation/build_antigen_candidates.py # ranks Coccidioides candidates (needs MMseqs2)
```

Scripts query QuickGO, UniProt and IEDB live, so the row counts above drift as those
databases are updated. Review decisions go in each folder's `manual_overrides.tsv`, never
by hand-editing a generated `*.tsv`, so the tables stay reproducible.

## Shared conventions
- `needs_review=yes` on any row a human has not confirmed (nearly all of them).
- Evidence levels: adhesins `E1`>`E2`>`E3`; hard negatives `N1`>`N2`>`N3`; antigens `A1`>`A2`>`A3`.
- `pmids` and `evidence_codes` carry provenance for every automatic call.
- GO `involved_in` does not encode direction: anti-adhesive and anti-biofilm proteins are
  annotated to the same terms as promoting ones. This is the single most common source of
  wrong labels here, and is why every GO-derived row is flagged for review.
