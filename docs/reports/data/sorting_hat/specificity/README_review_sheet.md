# Review sheet for the first Pfam sign-off

File: `review_sheet_first_signoff.tsv` (tab separated, 16 rows: CFEM 7, PA14 8, Hydrophobin 1).

For each row, fill the column `review_class` with one of:
`false_domain_hit` (Pfam finds the domain, but the protein is not a member of the family in the sense that matters for the call),
`uncharacterised_true_member` (a real member that is not on the draft list),
`receptor_like` (for example a CFEM domain in a membrane receptor),
or `unsure` (leave a note in `reviewer_note`).

`proposal_basis` is the table fact that left the row undecided. `link` opens the protein record. For the *C. albicans* rows the annotation is CGD's headline. The other rows of these three families (the ones that already have a proposed class) are in `pfam_proposals.tsv`; filter on the family column.

After review, tell the assistant the file name. The assistant records the classes, writes a `specificity_note` for each family, and prepares the `family_table.tsv` change for your sign-off. Nothing is set `active` until you say so.
