# Task 05. Controls for allergen similarity (`iuis_allergen_similarity`, `iuis_allergen_homolog`)

*Read `COMMON-RULES.md` first. Written 2026-10-06.*

## Goal

Build positive and negative controls for the two allergen evidence calls, and count how many
**negatives** exist. The negatives decide whether the module can ever have a specificity.

## Why this matters

- `iuis_allergen_similarity`: identity at least 35% over at least 80 aa (one local alignment).
  `iuis_allergen_homolog`: identity at least 70% and coverage at least 80%, or an allergen-specific
  Pfam hit. Both are **evidence**, not mechanism categories. IgE cross-reactivity is a hypothesis at most.
- Status today: `unvalidated`. A leave-species-out test gives sensitivity only.
- No negatives exist. Absence from IUIS means "never tested", not "not an allergen".
- The reference set holds the 30 *A. fumigatus* allergens. Counts in *A. fumigatus* show
  self-recognition.
- WHO/IUIS lists no *Coccidioides* allergen. The only Onygenales entries are four *Trichophyton*
  proteases (contact route).

## Known today

`docs/reports/2026-10-04-fungal-allergen-scoping.md`; files in `analysis/allergen_scoping/`:
- WHO/IUIS: 120 fungal allergen molecules from 31 species (measured 2026-10-04). 121 isoallergen
  rows, 116 with a protein sequence, 111 usable (4 fragments, 1 free-text entry skipped).
- UniProt keyword KW-0020 (Allergen), Fungi: 111 reviewed entries; 15 have "homolog" in the name
  (not all are demonstrated allergens).
- Not yet used: AllergenOnline (FARRP; fungal count not known), COMPARE (HESI; the download is a
  form that **emails** the files; owner action), Fungal Allergen Database (Iwate; maintenance
  unclear), IEDB (already used for the antigen table; has IgE epitope data).

## Positive controls

- IUIS allergens with a sequence: each is a positive for the 70%/80% rule only against **other**
  species (leave-species-out). Build the all-against-all table inputs; do not run BLAST here unless
  the owner asks.
- Extra positives from AllergenOnline and COMPARE for fungi that are **not** in IUIS, with the
  evidence of IgE binding (PMID). `stratum`: `iuis`, `allergenonline`, `compare`, `uniprot_kw0020`.
  Mark KW-0020 rows with "homolog" in the name `evidence_level=similarity_only`.

## Negative controls (the main task)

A negative needs an **independent reason**:
- IEDB IgE assay records with a **negative** outcome for a fungal protein (decision C2). Count
  these **first**: how many fungal proteins, how many species, how many independent clusters.
  Report the count before building anything else.
- Published IgE-binding tests of a protein with a negative result in sera of allergic patients
  (PMID, number of sera, quote).
- Homologs of an allergen (at least 35% identity) that were tested and did not bind IgE.
Do not use "not in IUIS" as a reason.

If there are fewer than 20 clusters of tested negatives, the module stays without a specificity.
State this plainly in `counts.md`.

## Independence and leakage

- Cluster the positives **and** the negatives together (30% identity). A negative in the same
  cluster as a positive is a conflict: list it in `rejected.tsv` with both rows.
- The reference set of the module is the IUIS fungal set. A positive that is in the reference set
  has `leakage=in_reference` for any measure that compares it with itself. Leave-species-out removes
  the species of the query from the reference. Say how.

## Size target

Positives per species group; number of tested-negative clusters. Report separately for
Onygenales, Eurotiales (*Aspergillus*, *Penicillium*), Pleosporales (*Alternaria*,
*Cladosporium*), Malassezia, and others.

## Do not

- Do not train or test an allergenicity model.
- Do not call IgE cross-reactivity a result.
- Do not submit the COMPARE form. Tell the owner it needs an email address.
- Do not use publisher texts in the repository.
