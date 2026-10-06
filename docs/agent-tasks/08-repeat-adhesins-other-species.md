# Task 08. Repeat-mediated adhesins and look-alikes in species outside Saccharomycotina

*Read `COMMON-RULES.md` and `03-repeat-mechanism-controls.md` first. Written 2026-10-06.*

## Goal

Find published, protein-level evidence of adhesins (and of look-alike proteins that are **not**
adhesins) in the fungal groups below. Add them as new rows in the format of
`data/controls/repeat-mechanism/curation_table.tsv`, with a mechanism class label. The aim is to
widen the positive and negative controls of the repeat call beyond *Candida* and *S. cerevisiae*.

## Why this matters

- The current table is skewed. The 76 rows that agent-A labels `2a` are 23 *C. albicans*, 13
  *C. glabrata*, 12 *C. auris*, 11 *S. pombe*, 8 *S. cerevisiae* and 2 *C. immitis*.
- Only 6 clusters of `2a` rest on a paper that states repeats for that specific protein. 68 rows rest
  on a family inference (agent-A table, 2026-10-06).
- The paper must say in which clades the repeat call was tested. Today the answer is
  Saccharomycotina only.
- A clade-specific model collapses outside its clade (ROC-AUC 0.60 across clades against 0.94 within).
  The same may hold for the detector. Without controls from other clades we cannot tell.

## Groups to examine (named by the owner)

*Pneumocystis* (limited number), *Cryptococcus*, *Rhodotorula*, *Penicillium* (also *Talaromyces*),
*Histoplasma*, *Exophiala*, *Wallemia*, *Hortaea*, *Knufia*, *Botrytis*, and plant-pathogenic fungi
(for example *Magnaporthe*/*Pyricularia*, *Colletotrichum*, *Fusarium*, *Zymoseptoria*, *Ustilago*/
*Mycosarcoma*, *Blumeria*).

## Size of the literature (PubMed, measured 2026-10-06)

Query form: `<genus>[Title/Abstract] AND (adhesin OR adhesion OR adherence)`. These counts are the
number of abstracts that contain the words. **They are not counts of adhesins.** Many hits are about
adhesion to host cells or to surfaces without a protein, and many adhesin papers do not use these
words in the abstract. The count only shows where to expect little.

| Group | PubMed hits | What it tells us |
|---|---|---|
| *Pneumocystis* (also "major surface glycoprotein") | 253 | Large literature; the Msg family is already one row in the table (`adhesin`, `other`) |
| *Cryptococcus* | 193 | Large; includes host-cell adherence. Adding "repeat" or "tandem" to the query gave **0** hits, so repeat evidence is not in abstracts |
| *Penicillium* | 65 | Moderate |
| *Botrytis* | 62 | Moderate |
| *Histoplasma* | 41 | Moderate; 2 *Ajellomyces* adhesin rows and 1 hard negative are in the table |
| *Rhodotorula* | 29 | Small |
| *Exophiala* | 3 | Very small |
| *Hortaea* | 1 | One paper (1996 record). Probably no protein-level adhesin work |
| *Knufia* | 0 | Nothing found by this query |
| *Wallemia* | 0 | Nothing found by this query |
| Plant pathogens (the seven genera above, with "fungus" or "fungal") | 148 | Moderate; includes appressorium and spore adhesion to leaves, much of it not a protein adhesin |

Already in the table for these groups (`data/controls/repeat-mechanism/curation_table.agent-A.tsv`):
*Pneumocystis* 1 adhesin; *Ajellomyces* (Histoplasma/Blastomyces) 2 adhesins and 1 hard negative;
*Pyricularia* 2 adhesins (for example the hydrophobin MPG1); *Mycosarcoma* 1; *Metarhizium* 2.
None for *Cryptococcus* with an organism (the *CFL1* row has no accession), *Rhodotorula*,
*Penicillium*, *Exophiala*, *Wallemia*, *Hortaea*, *Knufia* or *Botrytis*.

**Expect gaps.** For *Knufia*, *Wallemia*, *Hortaea* and probably *Exophiala* the literature may hold
no protein-level adhesin. Report "none found" with the queries that you ran. Do not fill the gap
from another species.

## What to do

1. **Search.** For each group, search PubMed (and PMC full text where open) for protein-level
   evidence. Vary the words: adhesin, adhesion, agglutinin, flocculation, hydrophobin, cell wall
   protein, surface protein, mucilage protein, attachment, appressorium, biofilm. Keep a log of every
   query in `search_log.tsv` (query, database, date, hits, hits read).
2. **Qualify each protein.** A protein qualifies as an adhesin row only if the paper shows adhesion or
   binding for **that protein**: deletion, heterologous expression, antibody blocking, purified
   protein binding. Mark `evidence_level`: `E1` direct experimental, `E2` weaker or family-level.
   Same definitions as `data/curated/adhesins/README.md`.
3. **Label the mechanism** with the class codes of `docs/TOOL-ARCHITECTURE.md` section 2 (`2a`,
   `2b-i`, `2b-ii`, `2b-iii`, `2c`, `2d`, `other`, `unknown`), `label_confidence`, and
   `evidence_basis` (`stated_in_paper`, `family_inference`, `background_knowledge`). The quote must
   support the label. For `2a`, the paper must state tandem repeats, repeat units or a repeat-rich
   region **for that protein** to be `stated_in_paper`. Do not label with a repeat detector.
4. **Get the sequence.** UniProt accession if one exists. If not, record the locus tag, the genome
   and annotation source, and take the sequence from that source. A protein without a retrievable
   sequence is listed in `rejected.tsv`.
5. **Find look-alike negatives in the same groups** (task 07 strata: `gpi_wall_enzyme`,
   `wall_hydrolase`, `mucin_sensor`, `st_linker_enzyme`, `repeat_non_adhesin`). A negative needs a
   paper that shows a non-adhesive function. Write them to `hard_negatives_other_species.tsv` with the
   task 07 columns.
6. **Cluster** new rows together with the existing 186 rows (MMseqs2, 30% identity, 50% coverage).
   Report, for each group, how many new clusters of each class you added, and how many new clusters of
   `2a` rest on `stated_in_paper`.
7. **Candidates are not controls.** If you find proteins in these genomes that look like repeat
   adhesins (secreted, repeat-rich) but have **no** adhesion paper, list them in `candidates.tsv`.
   They never enter `controls.tsv`. They are leads for the owner.

## Where to start in the repository

- `docs/model-review/2026-09-27-review-and-framework-plan.md` section 7.1 (positive families and
  anchor references) and 7.3 (Onygenales and Eurotiales).
- `analysis/curation/BASIDIOMYCETE_NOTES.md` (Basidiomycete adhesion notes; Cfl1 has adhesion evidence
  but the paralog Cpl1 is a secreted effector, so family membership does not imply adhesion).
- `analysis/step1_compare/curated_basidiomycota.tsv` (58 rows; GO-derived location labels for
  *Cryptococcus* and *Ustilago*, not adhesion labels).
- `data/curated/adhesins/literature_seeds.tsv`, `eurotiomycetes_seeds.tsv`, `hard_negative_seeds.tsv`.
- Publisher texts are not committed. Keep them in `_workdir/curation_texts_archive/`.

## Deliverables

`data/controls/repeat-adhesins-other-species/` with the files of `COMMON-RULES.md`, plus
`curation_table.other_species.tsv` (same columns as `curation_table.tsv`; `evidence_basis` included),
`hard_negatives_other_species.tsv`, `search_log.tsv`, `candidates.tsv`, and a table of
rows and clusters per group.

## Acceptance checks

1. Every row has a PMID and a quote that supports the label **and** the claim of adhesion for that
   protein.
2. No row appears in both `curation_table.other_species.tsv` and `hard_negatives_other_species.tsv`.
3. The log lists every query, including those with zero hits.
4. `counts.md` has rows, sequences and clusters per group, and says which groups returned nothing.

## Do not

- Do not count a PubMed hit as an adhesin.
- Do not copy a mechanism label from an ortholog in another species.
- Do not infer repeats from a database annotation or from running a detector.
- Do not add a row for a protein that the paper only mentions without testing.
