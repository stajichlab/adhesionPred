# Model card: custom Ecm33 profile HMM, `ecm33_all`, version 2026-10-23

Written 2026-10-24 by Claude Code (claude-sonnet-5-5). Every number below was read from a committed file (named in brackets). The model, `cutoff.json`, `SHA256SUMS` and the freeze were not changed.

## Words used (each defined once)

| Word | Plain meaning |
|---|---|
| Profile HMM | A statistical picture of a protein family made from an alignment of its members. It gives any protein a score (in bits) for how well it fits the family. |
| Ecm33 | A family of fungal cell-surface proteins, tied to the membrane by a lipid anchor (GPI anchor), that help build the cell wall. |
| Cutoff | The score at or above which a protein counts as a hit. |
| Report-only | The model's hits are listed in a separate section of the wall/fusion report. They change no call, tier, status or panel row. |
| Held-out | Proteins left out of training so they can be used to test. |
| Clade | A fungal group (for example Sordariomycetes). |
| UniRef50 | Proteins grouped so that members of a group are at least 50 percent identical; one protein stands for each group. |

## 1. Name and version

- Name: `ecm33_all` (HMMER profile NAME), folder `data/custom_models/ecm33_all_2026-10-23/`. Version = the folder name plus the file hashes in section 4. Profile length 500 positions, built from 481 sequences (`NSEQ 481`), made with HMMER 3.4 on 2026-10-23.
- It has no GA line (no curator-set threshold); the cutoff is in `cutoff.json`.

## 2. What it is for, and what it is not

It is a **report-only** model for the wall/fusion report: it lists proteins that look like members of the Ecm33 family, a family that Pfam's PF12454 mostly misses (Pfam GA reaches 28 of 187 PANTHER Ecm33 clusters). It is **not a call**: it does not decide wall, fusion or adhesin status, and it changes no tier, `family_table.tsv`, `family_measure.tsv` or `narrow_ok`. It measures family membership, not adhesion.

## 3. How it was built

[`docs/reports/2026-10-23-pr6-ecm33-build.md`, `docs/reports/data/pr6_ecm33_build_2026-10-23/`]

1. Members: phmmer (a one-protein search) with the Kiss Ecm33 Q5AGC4 (Candida albicans) against UniRef50 release 2025_03. A hit is a member if E <= 1e-5, the protein is fungal and has a fine clade. Result: **481 members in 19 clades** (the 681 raw hits included 195 non-fungal ones, never members; `members/summary.txt`).
2. Alignment: MAFFT 7.505 with L-INS-i, all 481 sequences.
3. Model: `hmmbuild --amino` (HMMER 3.4).
4. Pinned run: job 30287890 (task `all`), code commit `ba15e0c8`; `provenance.txt` lists the tool versions.

## 4. Cutoff rule and file hashes

- Cutoff **29.6 bits**, a **full-sequence score** (the score of the whole protein, the value `hmmsearch -T` tests). It was set **from false hits only** (no family member's score was used: `positives_used: false`) under the signed freeze `analysis/wall_fusion_truth/pr6_ecm33_freeze.json`, sha256 `a7553076d4963e4daedf3d7fdf0eba28b8e8349847b2644482fc2e854d7ce915` (the value in the build record and in `cutoff.json`). The rule: the lowest cutoff at which the four kept proteomes (A. fumigatus, C. albicans, C. immitis, S. cerevisiae; 32,491 proteins) have at most 5 unlabelled extra calls per 10,000 proteins (per proteome and in all, 16.2455) and no shared hard-negative group is called above a rate of 0.05. Result: 4 extra calls (all A. fumigatus), hard-negative rates 0.0 in all 9 groups. Margins: `docs/reports/2026-10-24-ecm33-margin-report.md`.

SHA256SUMS (verbatim; check with `shasum -a 256 -c SHA256SUMS` inside the folder):

```
9a6e2c812540b0b1892f081611bf4044ec674bc524d68f5c6fc3d16001c50937  model.hmm
1f51e55bcdeb78a1820a540f1762a661f39dc167836a2c1c502ce4ddb65ee486  cutoff.json
416f3c64dc208bd539064973d301b654be705a24dcd1533d0c8520c6770c44b0  provenance.txt
```

## 5. How to run it

Through the report (`analysis/wall_fusion_report/`, table `wall_fusion_custom_models.tsv`). The command it runs, as built in `hmmsearch_command()` of `analysis/wall_fusion_report/wall_fusion_report.py`:

```
hmmsearch -T 29.6 --noali --cpu N --tblout hits.tbl data/custom_models/ecm33_all_2026-10-23/model.hmm proteins.faa
```

`-T 29.6` keeps only proteins whose full-sequence score is at least 29.6 bits; there is no `--cut_ga`. Its hits are never merged with Pfam GA counts.

## 6. What the evaluations showed

Build record (11 held-out clade units, two independence modes each; unit models, not the `all` model):

- The model hits **469 of 469** held-out members (106 + 55 + 8 + 37 + 38 + 33 + 12 + 5 + 107 + 57 + 11), lowest score 86.9, unit cutoffs 19.8 to 45.9. Both true Ecm33s among the shared negatives are hit (`all` model: P38248 387.1, Q4WNS8 444.2).
- Pfam GA PF12454 hits far fewer (for example 0 of 107 in Saccharomycotina, 24 of 106 in Dothideomycetes); the iteration comparator (jackhmmer seeded by the training alignment) hits every member too, so **u = 0 in every unit and no unit is testable**.

Independent PANTHER evaluation [`docs/reports/2026-10-24-pr6-ecm33-independent-eval.md`]: 187 UniRef50 clusters of PANTHER Ecm33 families (725 proteins), 66 separated from training at 30 percent identity.

- Verdict: **not testable (u < 3): u = 1** (the one headline cluster missed by every method).
- Headline 15 clusters: custom **14/15** (95 percent interval 0.681 to 0.998); **seed-only jackhmmer 14/15, a tie**; one-pass phmmer 5/15; Pfam GA 1/15.
- All 187: custom 185, seed-only jackhmmer 185, Pfam GA PF12454 **28**, relaxed Pfam 48, phmmer 124.

## 7. Known limits

- The build test is **circular**: the held-out members come from the same Q5AGC4 search that supplied the training members.
- The independent test is **PANTHER agreement only**: PANTHER is itself a profile HMM, so it measures agreement between two models, not truth and not adhesion.
- The model **ties a single-protein jackhmmer** (14/15 and 185/187). What it adds is a fixed, hash-pinned, fast artefact, not more recall.
- Ecm33 is a **conserved family**, so the independent set is mostly close relatives of training: only 15 headline clusters are separated; 94 of 109 strict clusters share a 30 percent cluster with training. There is **no curated independent set**: all 13 Swiss-Prot family entries are in training clusters (`docs/reports/2026-10-24-ecm33-curated-ortholog-count.md`).
- **Untested outside fungi.**
- The cutoff margin is thin on the false-hit side (the next false hit scores 29.5; `2026-10-24-ecm33-margin-report.md`).
- 164 of 481 training members are in no PANTHER-called cluster (`purity.json`), a disagreement between two models that is unresolved.

## 8. What "additional development" could mean

Options, not decisions (`docs/reports/2026-10-24-pr6-decision-note.md`): a held-out set that is not the product of a search from Q5AGC4 (for example hand-curated Ecm33 entries from clades outside the three that Swiss-Prot covers); examining the 94 PANTHER clusters that share a 30 percent cluster with training; a new held-out set plus a new signed freeze if the cutoff is ever to move; applying the same test design to other PANTHER families.

## 9. Owner decision

2026-10-24: the owner decided to KEEP this model, report-only, and to record that it may need additional development (`README.md` in this folder). 2026-10-23: integrated as report-only and unvalidated (build record, section 6).
