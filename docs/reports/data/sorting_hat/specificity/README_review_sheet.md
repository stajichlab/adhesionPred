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
