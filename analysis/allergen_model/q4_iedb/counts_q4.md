# Q4: IEDB IgE-negative fungal records (measured 2026-10-10)

Source: IEDB Query API, `https://query-api.iedb.org/bcell_search`, run 2026-10-10.
Filters: `antibody_isotype = IgE`, `qualitative_measure = Negative` (328 records) and
`source_organism_iri_search` contains `NCBITaxon:4751` (Fungi). The positive query used the same
filters with `Positive` (178 records). Full pages returned (Content-Range 0-327/328 and 0-177/178).
Raw responses: `ige_Negative.json`, `ige_Positive.json`. Headers: `hdr_*.txt`.

All 328 negative records are linear peptides (`structure_type`). 327 are human hosts, 1 mouse.

## Counts

| Measure | Value |
|---|---|
| IgE-negative records (fungal, B-cell) | 328 |
| Distinct source proteins with a negative record | 16 |
| Source species with a negative record | 9 |
| Proteins with both IgE-negative and IgE-positive records (conflict) | 10 |
| Proteins with negative records only | 6 (from 5 species) |

## Negative-only proteins (6)

| UniProt | Name in IEDB | Species | Negative peptides |
|---|---|---|---|
| A0A084AFB8 | Sta c 3 | Stachybotrys chartarum | 23 |
| A0A177D5S5 | Polysaccharide synthase Cps1p (Fragment) | Alternaria alternata | 1 |
| A0ABY8EWX6 | Mala f 3 | Malassezia furfur | 1 |
| F2SLX7 | Tri r 2 | Trichophyton rubrum | 4 |
| P46075 | Asp f 5 | Aspergillus fumigatus | 1 |
| Q9Y750 | Mala f 4 | Malassezia furfur | 1 |

Five of the six are named WHO/IUIS allergens (Sta c 3, Mala f 3, Tri r 2, Asp f 5, Mala f 4). The sixth is
a fragment. An IgE-negative peptide on a known allergen is not a protein-level negative.

## Conflicts (10 proteins with negative and positive records)

Asp f 1, Asp f 2, Asp f 4, Asp f 6, Asp f 13, Pen c 13, Pen ch 13, Pen ch 18, Cur l 3, rhizopuspepsin.
Every one is a known or named allergen. A negative peptide does not make the protein a negative.

## Result against the spec

- Evidence-grade negative clusters: **zero** protein-level negatives are usable. The six negative-only
  proteins are allergens. The 10 conflicting proteins are allergens.
- The floor in spec section 3.3 is 20 independent negative clusters. Clusters are at most the number of
  proteins. There are 16 proteins with any negative record and 6 with negative records only. Both are below
  20, so the floor is not met under any counting rule. The cluster count itself needs MMseqs2 on HPCC, but
  the floor cannot be reached from these proteins, so the HPCC count does not change the answer.
- Consequence: the model has **no specificity claim** from IEDB. Sensitivity on positives is reported.
  Calls on weak Swiss-Prot negatives are reported and labelled optimistic (spec section 3.3).

## Limits

- IEDB negatives are epitope-level. This query covers IgE only. IgG and T-cell negatives were not counted.
- The query uses the source-organism lineage filter. It does not check for protein records that IEDB lists
  under a different organism name.
- Task 05's evidence rule still applies: a negative needs an independent reason. Peptide-level IgE negatives
  on named allergens do not qualify.
