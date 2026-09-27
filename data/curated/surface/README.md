# Surface glycoprotein table (stage 1) — DRAFT

`surface.tsv`, built by `analysis/curation/build_surface.py`. Every UniProt entry of a
reference genome with a surface/secretion keyword (Signal, GPI-anchor, Cell wall, Secreted),
cross-linked to the adhesin, biofilm and antigen tables.

**This is the table the whole review points at.** The shipped model already finds this
population; the open problem is splitting it.

| `adhesion_status` | n | use in training |
|---|---|---|
| `adhesin` | 75 | positives |
| `non_adhesin` | 63 | curated negatives |
| `non_adhesin_putative` | 1,257 | proposed negatives: an annotated catalytic activity (EC number or enzyme keyword) accounts for the protein and there is no adhesion evidence |
| `unknown` | 2,145 | unlabeled — treat as positive-unlabeled, never as negatives |

Genomes: the six tier-T5 reference proteomes plus *C. immitis* RS and *C. posadasii* C735
(for the antigen work).

## Caveats
- **Keyword-based localization.** UniProt keywords are largely automatic annotation for these
  genomes. Running SignalP 6.0 and NetGPI over the actual proteomes would give a better
  stage-1 population, and is the planned replacement (issue #15).
- **`non_adhesin_putative` is proposed, not curated.** A few adhesins are enzymatically
  active (*Paracoccidioides* gp43 is a glucanase), so this rule has real exceptions. Review
  before using these as hard negatives in a published model.
- **No adhesin labels at all for *Coccidioides*** (both species show 0 adhesins/0 non-adhesins).
  Stage 2 cannot be evaluated there; that is a genuine coverage gap, not a result.
- The `labels` column is multi-label on purpose, so a protein can be adhesin + biofilm + antigen.
