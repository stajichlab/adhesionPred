# Task 09. Lectin class: write the spec (host-interaction tier)

*Read `COMMON-RULES.md` first. Written 2026-10-10 by the session that ran the discovery pass. The owner agreed on 2026-10-10 that lectins become a **new candidate class in the host-interaction tier, with a spec first**. The owner wants help choosing the right set of models, training/truth data, literature and curation.*

## Goal

Write a **spec** (not code) for a lectin class of `cellsurface_sorting_hat`: which models define it, which proteins are in the truth set, what literature and curation are needed, how it is measured, and what the report says. Then get an independent review (Opus) of the spec. Stop there. Do not implement and do not change `categories.yaml`, `family_table.tsv` or any status file.

Branch: work on a new local branch from `fungi-scan-and-specs` (or its successor). No merges, no deletions, no PRs, no pushes unless the owner says so.

## Context (facts, with where to check them)

- **Tiers (owner decision 2026-10-09):** broad (cell wall and surface), narrow (adhesion), host interaction and immune evasion. Lectins go in the host-interaction tier. Workstream order is B (discovery pass, done), C (broad-tier truth set), A (rename and split categories). The lectin spec is separate from C and A. See `docs/HANDOFF-2026-10-09-overnight.md`.
- **Current coverage (measured):** only PA14 (PF07691) and GLEYA (PF10528) are active families that overlap lectin folds. They hit 9 and 12 of 54 E1 adhesins. The other lectin-like Pfam models are almost all in the unexplained pool (no current call recovers them): Jacalin PF01419 307 of 308 pool S proteins in pool U; Ricin_B_lectin PF00652 270 of 270; Fungal_lectin_2 PF18647 251 of 251; RicinB_lectin_2 PF14200 146 of 146; Laminin_G_3 PF13385 74 of 74; GBD_Y3 PF22803 70 of 70; Fungal_lectin PF07938 31 of 32; LDL PF21691 30 of 34 (653 proteomes, `_workdir/discovery/pools/`, local only; regenerate with `analysis/secretome_discovery/pools.py`).
- **No lectin truth exists.** None of the lectin models above hit the 54 E1 adhesins or the 137 N1 hard negatives (`_workdir/family_tests/lectins.hmm`, `lectin_*.domtbl`; method in `analysis/family_tests/README.md`). There are 74 Pfam 38.2 models whose name or description mentions lectin, agglutinin, galectin or glycan binding (`_workdir/discovery/pfam_names.tsv`, local; made with `awk` over `Pfam-A.hmm`).
- **Discovery pass:** three lectin families are enriched in Onygenales pool U at q < 0.05: PF07938 (6 of 27 genera), PF18647 (15), PF21691 (6). See `docs/reports/2026-10-10-secretome-discovery-pass.md` and `docs/reports/data/sorting_hat/secretome_discovery/families_U_onygenales.tsv`. This is an association with the clade, not a function.

## Literature (checked 2026-10-10; extend)

- **Lebreton et al. 2021**, *J Fungi* 7(6):453, PMID 34200153, [doi:10.3390/jof7060453](https://doi.org/10.3390/jof7060453), PMC8227253 (full text read). A bioinformatic survey. It makes 109 motif HMMs from UniLectin3D classes (194 fungal structures, 129 with a ligand), screens 1,223 MycoCosm species (`hmmsearch`, E-value 10, similarity score at least 0.25), and reports 28,691 putative lectins in the MycoLec database. The predictions are not experimentally validated and the biology focuses on Agaricomycetes. Structurally characterised fungal classes it names: galectin-like (18%), PA14 yeast adhesin (15%), AAL-like (12%), fungal fruit body lectin (11%), beta-trefoil (7%). It names house-keeping lectins (calnexin, calreticulin, P-type lectins, ERGIC-53) that are ER proteins. It states that pathogenic filamentous fungi make soluble lectins involved in host recognition that trigger lung inflammation (the references are in the paper; read them).
- **Erjavec et al. 2012**, *Trends Biotechnol* 30(5):259-73, PMID 22341093, [doi:10.1016/j.tibtech.2012.01.004](https://doi.org/10.1016/j.tibtech.2012.01.004). Abstract only (the publisher page returned 403). A review of mushroom proteins (lectins, enzymes, protease inhibitors, hydrophobins) for biotechnology. Mostly Agaricales; check for fungal lectins relevant to Ascomycota.
- Not yet read, to be found and verified through PubMed (do not cite from memory): the *Aspergillus fumigatus* fucose-binding lectin (FleA) and its role in host recognition and inflammation; Epa and Flo lectin-domain adhesins (PA14); the *Coccidioides*, *Histoplasma* and *Blastomyces* literature on carbohydrate-binding surface proteins; Dectin-1 and other host lectin receptors are the host side and are out of scope except where they define the fungal ligand.

## What the spec must decide

1. **Which models define the class.** Start from the UniLectin3D classes that occur in surface or secreted fungal proteins and map each to Pfam models (or keep a UniLectin-derived HMM where no Pfam model exists). Exclude the ER house-keeping lectins. Say for each class: model, cutoff (the Pfam GA, or how a custom cutoff is chosen and recorded), why it is surface- or host-relevant, and evidence level. Decide what to do with ubiquitous classes (Jacalin is in 251 of 619 background genera, Ricin-B in 216): they would flood the call.
2. **Gate.** Require a signal peptide? A mature TM helix? GPI? Reuse the tier rules of the tool (`docs/CLASSES.md`).
3. **Truth set (the owner wants help with this).** Candidates: fungal lectins with a solved structure and a ligand in UniLectin3D (experimentally backed; check licence and access); fungal lectins with experimental host-binding or phenotype evidence from the literature; hard negatives (house-keeping ER lectins, intracellular lectin-fold proteins, enzymes that share a lectin-like fold such as the PA14 glucosidases already in `data/controls/hard-negatives-shared/`). Follow the tier rules of the hydrophobin work (`docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md`): T4 review judgement is reported only and never used as truth; cluster at 30% identity; 20 clusters per class for `estimated`; leakage labels. Say how many clusters are realistic and say plainly if the answer is "smoke only".
4. **Curation workflow.** What an agent can curate, what needs the owner, rubric and decision file as in `analysis/hydrophobin_truth/`.
5. **Relation to existing calls.** PA14 and GLEYA are already active and sit in the narrow (adhesion) tier. Which lectin classes are adhesion-like (narrow tier) and which are host-recognition lectins (host-interaction tier)? How does the new class avoid double counting with `wall_family_domain` and `surface_attachment_candidate`?
6. **Measurement and report.** Calls, status (`unvalidated`, `smoke`, `estimated`) and what the report shows. How the discovery-pass families (PF07938, PF18647, PF21691) are handled.
7. **Work plan** with tasks, checks and stop points, and the owner decisions (as a table with defaults). Use the format of `docs/superpowers/specs/2026-10-09-sowgp-like-calls-design.md` (facts table with file or PMID for every number).

## Do not

- Do not cite a paper you have not found in PubMed or read. Give PMID and DOI. If a source cannot be reached, say so.
- Do not use a protein as truth because the tool calls it, or because a Pfam model hits it (that is circular for a family-defined class).
- Do not edit `data/curated/adhesins/adhesins.tsv`, `docs/agent-tasks/05-allergen-controls.md` or `docs/agent-tasks/COMMON-RULES.md` (the last two have uncommitted changes from another session).
- Do not run SLURM jobs without checking the queue; the group CPU limit is shared (the discovery jobs waited about 7 hours).

## Output

- `docs/superpowers/specs/<date>-lectin-class-design.md` (revision 1), then `...-review-1.md` from an independent Opus review run, then revise until a review passes or a blocker needs the owner.
- A short list of owner questions, one at a time, with a recommended answer for each.
- Update `docs/ROADMAP-2026-10-08-DRAFT.md` with the new milestone only if the owner agrees.
