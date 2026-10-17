# Lectin validation truth set (controls)

Built 2026-10-11 by claude-sonnet-5-5. Report: `docs/reports/2026-10-11-lectin-validation.md`. Log: `analysis/lectin_validation/LOG.md`.

This set is for a measurement of the lectin-like Pfam models. It is not a calibration set. It makes no call and no status (owner decision D1, option (c)).

## Files

| File | Content |
|---|---|
| `controls.tsv` | 733 rows. COMMON-RULES columns, plus `label_class`, `tier`, `h_flag`, `protein_name`, `length`, `strain_taxid`, `unilectin_class`, `pdb_ids`, `pdb_sugar_comps`, `uniprot_status`, `negative_evidence`, `phylum`, `selection_route`, `cluster30`, `seed_hits` |
| `sequences.faa` | Sequences, header = accession |
| `clusters.tsv` | accession, cluster representative (MMseqs2 30% identity, 50% coverage) |
| `rejected.tsv` | 82 candidates not added, with the reason |
| `counts.md` | Rows, sequences and clusters per class |
| `unverified_notes.md` | Claims not checked (copy of `../unverified_notes.md`) |

## Classes

Positive: `structure_backed`, `host_binding_h1a`, `chitin_masking_h1b`, `literature_binding`, `known_gap` (primary), `pa14_owned_elsewhere` (apart).
Negative with a reason: `er_housekeeping`, `enzyme_lectin_like_domain`, `n1_shared`.
Negative without a reason: `assumed`. Never mix it with the reasoned classes.

## Versions and commands

UniProtKB release 2026_03 (fetched 2026-10-11). UniLectin3D API and RCSB GraphQL (2026-10-11). InterPro seed alignments (versions equal to Pfam 38.2). Pfam 38.2 (`/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`). HMMER 3.4. MMseqs2 17-b804f. Rebuild: `analysis/lectin_validation/run_all.sh ROOT` (network steps are separate: `fetch_unilectin.py`, `fetch_uniprot.py`, `fetch_seeds.py`, `select_assumed.py`; their outputs are in `../data`).
The raw Swiss-Prot lists (`data/raw/swissprot_fungi.tsv.gz`, `.faa.gz`) are not in git. Fetch them with the UniProt stream query `taxonomy_id:4751 AND reviewed:true` (fields in `select_assumed.py`).

## Weakest rows

The 39 `structure_backed` rows (tier T2: the ligand was not checked per structure), Vd2LysM (G2X4U8), the 12 `er_housekeeping` rows with evidence by similarity, and all `enzyme_lectin_like_domain` rows. `assumed` rows have no reason.

## Location

COMMON-RULES asks for `data/controls/<task-name>/`. This set stays in `analysis/lectin_validation/controls/` because it is a measurement table and not a calibration set. Exception recorded in the report, section 9.
