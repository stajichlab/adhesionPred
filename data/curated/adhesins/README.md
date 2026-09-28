# Curated adhesin and hard-negative candidates (DRAFT: needs expert review)

`adhesins.tsv` is a **draft** label set for training and evaluating the adhesion classifier
(issues #13, #14). Build it with `python analysis/curation/build_adhesin_candidates.py`, which queries QuickGO
and UniProt live, so re-running can change rows as those databases update.

## Inputs
| file | what |
|---|---|
| QuickGO (live) | experimental (ECO:0000269 descendants: IDA, IMP, IGI, EXP, IEP) annotations in Fungi to *cell adhesion*, *adhesion of symbiont to host*, *flocculation*, *cell adhesion involved in single-species biofilm formation*, *agglutination involved in conjugation*, *cell-cell adhesion* (with is_a/part_of descendants) |
| `literature_seeds.tsv` | adhesins with literature evidence but no GO experimental term (BAD1, CalA, CspA, Cfl1, gp43, Msg, Mad1/2, CpAls7, Als4112) |
| `hard_negative_seeds.tsv` | non-adhesive surface/secreted proteins that share adhesin architecture |
| UniProt Pfam (live) | members of the six reference proteomes carrying an adhesin-defining Pfam domain found in an E1 adhesin (Flocculin, PA14, Flo11, Candida_ALS/_N, Hyr1) |
| `manual_overrides.tsv` | documented judgment calls, applied last (e.g. YWP1 is anti-adhesive) |

## Columns (main ones)
- `cls`: putative class
  - `adhesin`: train as positive once confirmed
  - `hard_negative`: train as negative
  - `surface_other_adhesion_phenotype`: surface enzyme or structural protein whose mutant affects adhesion. Likely indirect; hold out or use as hard negatives after review.
  - `indirect_regulator`: non-surface protein; exclude from training
- `evidence_level`:
  - adhesins: **E1** direct experimental adhesion or binding evidence; **E2** family/domain membership or weaker evidence (biofilm-only, antibody-blocking, no deletion); **E3** domain hit contradicted by function or pseudogene-like (check)
  - hard negatives: **N1** characterized non-adhesive function; **N2** adhesin-like architecture and no adhesion evidence
- `pmids`, `evidence_codes`, `go_terms`, `source`: provenance
- `needs_review`: `no` only for reviewed UniProt entries with E1 adhesin evidence
- `in_reference_genome`: protein belongs to one of the tier-T5 reference proteomes (S288C, *C. albicans* SC5314, *N. glabratus* CBS138, *C. auris* B8441, *S. pombe*, *A. fumigatus* Af293)

## Known limitations of this draft
- GO annotations do not encode direction of effect (`involved_in` covers anti-adhesive proteins such as YWP1). Every GO-derived row needs a human check.
- Taxonomic coverage is heavily skewed to Saccharomycotina. There are **no curated positives for Chytridiomycota, Mucoromycota beyond CotH, or most Pezizomycotina and Basidiomycota**. The model cannot be evaluated in those lineages yet.
- The CalA and CspA accessions were assigned from locus tags recalled from memory (flagged in `evidence_summary`). Cfl1 has no resolved UniProt accession.
- PA14 is shared with non-adhesive glycosidases. PA14-only family hits are E3 or hard negatives.
- Hard negatives are almost all *S. cerevisiae*. Orthologs in the other reference genomes still need adding.

## Review guidance
Reviewers should confirm or change `cls` and `evidence_level` for rows with `needs_review=yes`,
prioritizing (1) all `adhesin` E1/E2 rows, (2) `surface_other_adhesion_phenotype`, (3) `hard_negative` N2.
Record changes in `manual_overrides.tsv` (not by editing `adhesins.tsv`), so the table stays reproducible.
