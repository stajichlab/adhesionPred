# HsbA family (PF12296): what it is and why the call is only relabelled (2026-10-26)

**Observation (from a review of the aggregate report).** PF12296 is the Pfam HsbA family. On 2026-10-08 it was grouped with the hydrophobins ("surface-active"). In the ten proteomes it hits 9 *A. fumigatus* proteins and none elsewhere, and the hand-check sample
shows that at least one hit is Q4WHJ8, "Cell wall mannoprotein 1" (Afmp1, the Mp1p-type galactomannoprotein family; the blind second opinion read the entry as a Ser/Thr-rich, C-terminally GPI-anchored wall protein
localised to hyphal and conidial walls). Read that way the family is a wall galactomannoprotein family, not a hydrophobin-like one. (Evidence level: one reviewed-looking entry and a model's literature reading; PMIDs not independently checked for the specific claim.)

**Owner decision 2026-10-26.** Reclassify the family as a wall family.

**What was done.** Only the description: `docs/class_descriptions.yaml` (and the generated `docs/CLASSES.md`) now describe `hsba_domain` as a wall-family hit, not a surface-active call. No behaviour changed.

**What was not done, and why.** Counting the family under `wall_family_domain` was tried and reverted:
* `wall_family_domain` would read two modules (`pfam_adhesion`, `pfam_hsba`). About 25 tests and the calibration paths assume it reads `pfam_adhesion` only (the `truth --module pfam_adhesion` command no longer fits a two-module call; the label, gap and golden tests change);
* moving PF12296 to the `pfam_adhesion` module leaves `pfam_hsba` without an active family, which makes the module `unavailable` and `hsba_domain` not assessable, and an unknown input in the `surface_attachment_candidate` OR makes the composite not assessable for every protein, so the call would have to be removed from the composites and the `other_*` lists as well;
* `data/sorting_hat/family_table.tsv` is frozen by a test (`test_family_table_and_family_measure_unchanged_from_origin_main`); changing its class text needs that test and the freeze note updated by the owner;
* any edit of `categories.yaml`, even a comment, changes its sha256 and makes the status files stale.

**Practical effect today.** None on the calls: `surface_attachment_candidate` already includes the HsbA hits, and `cell_wall_adhesion_candidate` reads the narrow call only (PF12296 has `narrow_ok=no`). To make the HsbA proteins count in `wall_family_domain` itself, a config change is needed:
bundle it with the next config change so the status files are re-measured once.
