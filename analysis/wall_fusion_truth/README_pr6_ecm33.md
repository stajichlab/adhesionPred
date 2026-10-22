# PR 6 steps 6.1 and 6.2 for ecm33: build and leave-clade-out (DRAFT, nothing run on HPCC)

Owner decisions of 2026-10-22 that bind this work (spec `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md`, sections 4 to 7b): **ecm33 only** is built
(KS0086, Q5AGC4, wall slice); hyd1_like is a known gap under the strict reading; aga1 and ece1 wait for a drift check. One model per
homologous family from the **conservative one-pass phmmer set** (never the jackhmmer set). Wall = leave-clade-out per held-out clade at the fine level.
Training rule: at least 5 sequences in at least 2 clades outside the held-out unit (provisional). Cutoff by the hydrophobin rule:
at most 5 unlabelled extra calls per 10,000 proteins over the four kept proteomes (32,491 proteins: 16.2455, about 16 in all) and a fixed hard-negative
rate, both written at a **freeze before any held-out score exists**. Comparators on the same units: Pfam GA, relaxed-GA Pfam, a jackhmmer iteration.
Custom models are report-only. No call, status, tier, `narrow_ok`, `family_measure.tsv` or `family_table.tsv` row is written by anything here.

**Decided by the owner on 2026-10-22** (now fixed in the freeze proposal, no longer open): (1) **P38248** is left out of the hard-negative rate and reported on its own line; the shared negative set is not changed for anyone else. (2) Independence: **`cluster30` is the HEADLINE** and `none` (only the held-out members are removed) is a **second, descriptive column**; both run for every held-out unit and are reported side by side. (3) **`min_u` = 3**. Everything else in the freeze file is still a PROPOSAL and keeps its label (listed below).

## Files

| File | Purpose |
|---|---|
| `pr6_ecm33_members.py` | stdlib; the members and the held-out units from the committed round 2 table (`gap_homologs.tsv.gz`). Ids, clades, scores only. |
| `pr6_ecm33_extract.sbatch` | HPCC job (DRAFT): runs `pr6_ecm33_members.py`, streams UniRef50 once with `seqkit grep`, stages `ecm33_members.faa` + sha256. |
| `pr6_ecm33_build.py` | stdlib; `manifest` (training set of a unit under one independence mode), `check-alignment` (manifest and a pre-registered occupancy check), `cutoff` (false-hit side only; refuses without a frozen freeze). |
| `pr6_ecm33_build.sbatch` | HPCC array job (DRAFT): per (unit, independence mode) MAFFT L-INS-i, `hmmbuild`, false-hit searches, cutoffs, then the held-out searches and the comparators. |
| `pr6_ecm33_report.py` | stdlib; the held-out report (headline and descriptive modes side by side, P38248 on its own line); refuses without the frozen freeze and a matching sha256 in every cutoff file. |
| `pr6_ecm33_freeze_PROPOSED.json` | the PROPOSED freeze (`"frozen": false`): the decided items and the unchanged proposals. The owner freezes it (steps in the run order). |
| `tests/wall_fusion_truth/test_pr6_ecm33.py` | fixtures, plus the real round 2 table for the member counts, plus static checks of both jobs. |

## What the members are, and are not

Q5AGC4 gives 681 phmmer hits in UniRef50 (E 1e-5): 195 non-fungal (never members), 5 fungal without a fine clade (listed in `excluded.tsv`, **not a clade**, not
counted), and **481 members in 19 fine clades** (480 UniRef50 representatives + the query itself, read from the Kiss control set). Held-out units
(a fine clade with at least 5 members, training rule met): **11**; too small to test: 8.

| Held-out unit | members | | too small (still training members) | members |
|---|---|---|---|---|
| Saccharomycotina | 107 (106 + the query) | | Lecanoromycetes | 4 |
| Dothideomycetes | 106 | | Geoglossomycetes | 2 |
| Sordariomycetes | 57 | | Agaricomycetes, Basidiobolomycota, Chytridiomycota, Mortierellomycota, Xylonomycetes, Zoopagomycota | 1 each |
| Eurotiomycetes | 55 | | | |
| Leotiomycetes | 38 | | | |
| Kickxellomycota | 37 | | | |
| Mucoromycota | 33 | | | |
| Orbiliomycetes | 12 | | | |
| Taphrinomycotina | 11 | | | |
| Glomeromycota | 8 | | | |
| Pezizomycetes | 5 | | | |

**A member is a homolog of Q5AGC4 by sequence search, not an independently curated Ecm33.** The held-out "truth" of a unit is defined by the same
search that supplies the training members, and the members of a clade are related, not independent clusters. Recall on a unit therefore measures whether a
model built from the other clades recovers the family in an unseen clade. It is **not** evidence about adhesion and **not** call sensitivity, and **every
cell is `descriptive`**: counts and a recall fraction, no interval, no `estimated`. (The Kiss clusters and the floor of 5 clusters of decision 4 do
not apply to UniRef50 members; the only Kiss member is Q5AGC4 itself.)

## Run order (each HPCC step needs the owner's go-ahead; exact commands, NOT run)

0. **The freeze commit (the OWNER does this; nothing here freezes anything).**
   ```
   cd analysis/wall_fusion_truth
   cp pr6_ecm33_freeze_PROPOSED.json pr6_ecm33_freeze.json
   # edit pr6_ecm33_freeze.json: set "frozen": true and "owner_signature" (name and date); change any unchanged proposal you do not accept
   git add pr6_ecm33_freeze.json && git commit -m "freeze: ecm33 build and leave-clade-out (owner)"
   ```
   The build job stops if `pr6_ecm33_freeze.json` is missing or not `"frozen": true`; the report stops unless every cutoff file records the sha256 of that file.
1. Update the HPCC clone to the merge commit of the PR that holds these scripts and the freeze.
2. **Extract** (one job; the file of the round 2 record, 28,781,193,344 bytes, is checked by size, not hashed again):
   ```
   sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<sha>,OUT=$W/ecm33_extract \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_ecm33_extract.sbatch
   ```
   Resources are an UNMEASURED guess (batch, 4 cores, 8 GB, 2 h): one pass of `seqkit` over 28.8 GB should take minutes to a few tens of minutes.
3. **Build** (array: **2 tasks per held-out unit (the headline mode `cluster30`, then the descriptive mode `none`) + 1 task for the `all` model = 23 tasks, 0 to 22**; at most 4 at a time):
   ```
   sbatch --array=0-22%4 --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<sha>,EXTRACT=$W/ecm33_extract,OUT=$W/ecm33_build,PFAM_RELEASE=38.2 \
     -o $W/logs/%x.%A_%a.log -e $W/logs/%x.%A_%a.log $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_ecm33_build.sbatch
   ```
   (the size is `2 x units + 1`; with another unit list compute `N=$(( 2 * $(awk -F'\t' 'NR>1 && $7=="unit"' $W/ecm33_extract/members/units.tsv | wc -l) ))` and use `--array=0-$N%4`).
   **Cost.** A single-mode build would have 12 tasks (11 units + `all`); with both modes every unit is built twice, so up to **23 tasks**, each of the same kind (clustering, MAFFT L-INS-i, `hmmbuild`, searches). The run time is **still unmeasured**: at 4 tasks at a time, 23 tasks take about
   6 waves instead of 3, so the wall time roughly doubles if every task runs to the 12 h request. A cut: the descriptive task of a unit is **not rebuilt** when the headline rule removed no neighbour of that unit (the two training sets are then identical; the task writes
   `SAME_AS_HEADLINE` and exits in seconds), so the real number of builds is between 12 and 23. Results are written to `<out>/<unit>/<mode>/`, `all` to `<out>/all/none/`.
   Inside each task the order is fixed: MMseqs2 clustering of the members (30 percent identity, 50 percent coverage) -> manifest (headline manifest first, for the count of excluded neighbours; then the manifest of the task's mode) -> MAFFT -> alignment checks (the held-out unit
   and, in `cluster30`, its 30 percent neighbours are absent from the alignment: the same check for both modes) -> `hmmbuild` -> searches over the four kept proteomes and the 161 shared negatives -> the two cutoffs -> **only then** the held-out searches. L-INS-i on up to 480 sequences is slow: the 12 h request is a guess.
4. **Report** (local, after the results are copied back):
   `python3 -I analysis/wall_fusion_truth/pr6_ecm33_report.py --results <ecm33_build> --freeze analysis/wall_fusion_truth/pr6_ecm33_freeze.json --units <members/units.tsv> --out <dir>`: `recall_by_unit.tsv` has two rows per unit, labelled headline (`cluster30`) and descriptive (`none`); `side_by_side.tsv` puts the two next to each other with the number of excluded neighbours.
5. The results go into their own record PR (tables and provenance, no sequences).

## The freeze: what is decided and what is still a proposal

**Decided by the owner (2026-10-22), written in the freeze file:**
1. **P38248** (*S. cerevisiae* Ecm33, a shared hard negative of stratum `gpi_wall_enzyme` that is a true family member) is left out of the hard-negative rate (`negatives_ignored: ["P38248"]`) and reported on its own line (`ignored_negatives_line.tsv` (one line per accession of `negatives_ignored` and mode): its custom score, the cutoff and whether it is hit, per unit and mode). The shared negative set is not changed for anyone else.
2. **Independence**: HEADLINE `cluster30` (spec 4.4), SECOND DESCRIPTIVE column `none` (`independence_mode`, `independence_second_column`). The manifest check that the held-out unit is absent from the alignment runs for both modes. Read the count of excluded neighbours per unit before judging any result: for a conserved family the 30 percent clusters can span several clades, and a unit whose headline training set falls under the minimum is `not_buildable` there while its descriptive column may still be built.
3. **`min_u` = 3**: u is the number of held-out members that every comparator misses; a unit is testable when u >= 3 and the custom model must recover at least ceil(u/2) of them (spec 4.6; the hydrophobin run had u = 0).

**Still PROPOSALS (unchanged, labelled in the file under `proposals_unchanged`):** cutoff floor 0.0; at most 5 unlabelled extra calls per 10,000 proteins (5 x 32,491 / 10,000 = 16.2455, about 16 in all); hard-negative rate 0.05; at least 5 members for a unit; the training rule (5 sequences in 2 clades);
the alignment check (at least 100 columns with a residue in at least half of the sequences, median coverage of those columns at least 0.6; a stop condition, the hydrophobin cysteine check does not apply).

## Pre-freeze check of the hard negatives (do this BEFORE the freeze commit)

**Why.** `negatives_ignored` of the freeze keeps true Ecm33 family members out of the hard-negative rate (a correct model should hit them). The freeze lists only
`P38248` (*S. cerevisiae* Ecm33). The shared negatives also hold `Q4WNS8` (*A. fumigatus* ecm33, stratum `gpi_wall_enzyme`), and a member under another name cannot be
excluded by reading names. If one stays in the rate, the cutoff is pushed above the score of a true Ecm33. The owner chose (2026-10-22) to list ALL members before
freezing, with the freeze's own definition of a family member (`labelled_rule`): `phmmer` of Q5AGC4 finds it at E <= 1e-5.

**What the job does** (`pr6_ecm33_negatives_precheck.sbatch`, DRAFT, tiny: partition short, 2 cores, 2 GB, 30 min). It checks the pin and the staged FASTA hash, takes Q5AGC4
from `kiss_sequences.faa`, runs `phmmer` against the 161 shared negatives at E 1e-5 (the definition) and, separately and descriptively, at E 1e-3, and writes the tables of
`pr6_ecm33_negatives_precheck.py`. It does not read or change the freeze.

```
# 1. (owner go-ahead) update the HPCC clone to the merge commit of the PR that holds this script
# 2. (owner go-ahead) submit; no new staging: kiss_sequences.faa is already in $W/in
B=/bigdata/stajichlab/jstajich/projects; W=$B/adhesionPred/_workdir/kiss_cluster_split
sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/ecm33_negcheck \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_ecm33_negatives_precheck.sbatch
```

**Output** in `$OUT`: `negatives_family_table.tsv` (accession, gene, species, stratum, E-value, score, length, whether the freeze lists it), `weak_E1e-3.tsv` (descriptive only,
never in the list), `proposed_negatives_ignored.json` (the exact one-line edit), `summary.txt` (counts, per stratum, NEW accessions), `provenance.txt` and the two `.tbl` files.
A `LOUD:` line appears if P38248 or Q4WNS8 is not found (both are expected), if a hit is not in `controls.tsv`, or if an accession the freeze lists is not found.

**What the owner does with it.** Read `summary.txt`. If the list is what you expect, put the proposed line in the freeze file before signing it (replace the one line
`"negatives_ignored": [...]`). The report then writes one line per listed accession (`ignored_negatives_line.tsv`), and the build leaves all of them out of the rate. The list
is a proposal made from the data: a protein with a weak hit (the descriptive table) is NOT added, and you decide any borderline case.

## Checks and expected outputs

- `extract`: `members/` (members.tsv, excluded.tsv, units.tsv, member_ids.txt, summary.txt), `ecm33_members.faa` (+ sha256), `provenance.txt`. The job fails if
  UniRef50 is not the file of the round 2 record, if any member id is not found, or if a header is duplicated.
- `build`, per (unit, mode): `manifest.json` (held-out ids, training ids with clade and 30 percent cluster, the excluded neighbours), `train.aln.fa`, `model.hmm`, `cutoff.json`
  and `relaxed_cutoff.json` (each with the freeze sha256 and `positives_used: false`), the held-out hit tables, `provenance.txt`. A unit that is not buildable keeps only its manifest.
- `report`: `recall_by_unit.tsv` (per unit and mode, role headline or descriptive: custom, Pfam GA, relaxed Pfam and iteration recall; u; whether the custom model recovers ceil(u/2) of u), `side_by_side.tsv`, `extra_calls.tsv`
  (extra calls against the freeze limit), `hard_negative_calls.tsv` (rates against the frozen rate, P38248 not in them), `ignored_negatives_line.tsv` (one line per accession of `negatives_ignored` and mode), `summary.txt`.

## Open items I could not confirm (nothing was run)

- **MAFFT**: the module `mafft/7.505` is listed on HPCC; the options are the precedent's (`--localpair --maxiterate 1000`); the run time on 480 sequences is unmeasured.
- **seqkit streaming 28.8 GB**: module `seqkit/2.4.0` is listed; the 2 h request is a guess.
- **`esl-reformat`** (Easel, bundled with HMMER) for the Stockholm input of the jackhmmer comparator, and whether `jackhmmer` takes that alignment as its query as written.
- **`hmmfetch -f` by the name `Ecm33`** (the name of PF12454 in `wall_fusion_models.tsv`); the job checks the fetched model's accession.
- **Proteome FASTA paths** (`$DATA_ROOT/_workdir/sorting_hat/<name>.faa`): not re-checked on HPCC; the job stops if a count differs from `inputs.tsv`.
- **MMseqs2** `MMseqs2/17-b804f` (worked in the Task 5.3 job).
- The extra-call count depends on the labelled set (phmmer of Q5AGC4 against each proteome, E 1e-5): a proteome protein that is a real family member but not found by that search counts as an extra call.

## Next owner go-aheads

The freeze commit; the clone update; the extract job; the build array; copying the results back; the record PR. aga1 and ece1 are not part of this build; hyd1_like is a known gap
(its `buildable` verdict in the round 2 table rests on the "(no clade)" bucket of the summary script, which this work does not count as a clade).
