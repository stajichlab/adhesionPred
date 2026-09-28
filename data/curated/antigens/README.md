# Antigen tables (DRAFT)

Purpose: **immunodiagnostics and serology** — finding proteins the host immune system
recognizes, which can be produced recombinantly and tested against patient sera. Focus on
*Coccidioides*, with other fungi included for transferable evidence.

## `antigens.tsv` — curated immune evidence
Built by `build_antigens.py` from the [IEDB](https://www.iedb.org) Query API, which curates
published T cell, B cell and MHC-elution assays. One row per protein with a resolvable
UniProt accession, recording assay counts by type.

| class | meaning |
|---|---|
| `antigen_tcell` / `antigen_bcell` / `antigen_both` | which assay types support it |

Evidence level is by assay count only (`A1` ≥5, `A2` 2–4, `A3` 1) and is **not** a quality
judgment. Coverage is thin for *Coccidioides* (5 proteins), which is exactly why the
candidate ranking below exists.

## `coccidioides_candidates.tsv` — ranked candidates for production
Built by `build_antigen_candidates.py`: all 1,069 *C. immitis* RS and *C. posadasii* C735
surface/secreted proteins, scored by transferable evidence:

| signal | points |
|---|---|
| homologous to a curated IEDB fungal antigen | +3 |
| already a curated *Coccidioides* antigen (positive control) | +2 |
| conserved between *C. immitis* and *C. posadasii* (a test must detect both) | +2 |
| GPI-anchored or cell-wall | +1 |
| signal peptide (reaches host fluids) | +1 |
| no close human homolog | +1 |
| close human homolog (cross-reactivity risk) | −2 |

**The ranking reproduces known biology it was not told about.** The top scores recover the
proline-rich antigen (Ag2/PRA) family, 1,3-β-glucanosyltransferases (Gel1), CFEM-domain
proteins and subtilisin-like proteases — all described *Coccidioides* antigens. The human
screen independently flags the classic cross-reactive decoys (BiP/HSP70, peptidyl-prolyl
isomerases, calnexin, protein disulfide isomerase), which are immunogenic but poor
diagnostic markers.

**What this is not.** The score ranks proteins for *wet-lab prioritization in diagnostics*.
It says nothing about whether a protein is protective; that requires immunological testing.
Sequence homology to an antigen also does not transfer epitopes: conservation of the region
actually recognized still has to be checked per candidate.
