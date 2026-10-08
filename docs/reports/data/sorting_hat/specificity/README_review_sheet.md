# Review sheet for the first Pfam sign-off

File: `review_sheet_first_signoff.tsv` (tab separated, 16 rows: CFEM 7, PA14 8, Hydrophobin 1).

For each row, fill the column `review_class` with one of:
`false_domain_hit` (Pfam finds the domain, but the protein is not a member of the family in the sense that matters for the call),
`uncharacterised_true_member` (a real member that is not on the draft list),
`receptor_like` (for example a CFEM domain in a membrane receptor),
or `unsure` (leave a note in `reviewer_note`).

`proposal_basis` is the table fact that left the row undecided. `link` opens the protein record. For the *C. albicans* rows the annotation is CGD's headline. The other rows of these three families (the ones that already have a proposed class) are in `pfam_proposals.tsv`; filter on the family column.

After review, tell the assistant the file name. The assistant records the classes, writes a `specificity_note` for each family, and prepares the `family_table.tsv` change for your sign-off. Nothing is set `active` until you say so.

## Ortholog columns (added 2026-10-07)

For the W72310 and A1163 rows, `af293_ortholog` is the best reciprocal BLAST hit in the Af293 UniProt proteome (identity and coverage), `af293_ortholog_annotation` is its UniProt name, and `af293_ortholog_hit_class` is the class proposed for that Af293 protein in the same family ("no hit in this family" if the ortholog has no hit). They help for "hypothetical protein" rows. A high identity with full coverage means the same protein in another strain. One row has low coverage (KAK9638745.1 against BGLK, coverage 0.605), so read it with care.

## GPI columns (added 2026-10-07)

`predgpi_call`, `predgpi_fpr` and `predgpi_omega_site` come from PredGPI (module `predgpi/202001`, the same wrapper as Phase B, `analysis/step1_compare/jobs/predgpi_scores.py`). The call is `highly_probable` (estimated false positive rate at most 0.0015), `probable` (at most 0.005), `weakly` (at most 0.01) or `none`. It is a prediction, not evidence of a GPI anchor. PredGPI failed with an internal error on YAL064C-A (TDA8) and YHR213W, so those two rows are not scored. `cimg00693_vs_calb.tsv` is a one-way BLAST of the RS protein XP_001246922.1 (CIMG_00693) against *C. albicans* (best hit PGA7, 28.9% identity over 121 aa, e-value 5.7e-7). It is not a reciprocal test.
