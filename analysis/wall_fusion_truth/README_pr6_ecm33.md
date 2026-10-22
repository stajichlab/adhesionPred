# PR 6 steps 6.1 and 6.2 for ecm33: build and leave-clade-out (DRAFT, nothing run on HPCC)

Owner decisions of 2026-10-22 that bind this work (spec `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md`, sections 4 to 7b): **ecm33 only** is built
(KS0086, Q5AGC4, wall slice); hyd1_like is a known gap under the strict reading; aga1 and ece1 wait for a drift check. One model per
homologous family from the **conservative one-pass phmmer set** (never the jackhmmer set). Wall = leave-clade-out per held-out clade at the fine level.
Training rule: at least 5 sequences in at least 2 clades outside the held-out unit (provisional). Cutoff by the hydrophobin rule:
at most 5 unlabelled extra calls per 10,000 proteins over the four kept proteomes (32,491 proteins: 16.2455, about 16 in all) and a fixed hard-negative
rate, both written at a **freeze before any held-out score exists**. Comparators on the same units: Pfam GA, relaxed-GA Pfam, a jackhmmer iteration.
Custom models are report-only. No call, status, tier, `narrow_ok`, `family_measure.tsv` or `family_table.tsv` row is written by anything here.

## Files

| File | Purpose |
|---|---|
| `pr6_ecm33_members.py` | stdlib; the members and the held-out units from the committed round 2 table (`gap_homologs.tsv.gz`). Ids, clades, scores only. |
| `pr6_ecm33_extract.sbatch` | HPCC job (DRAFT): runs `pr6_ecm33_members.py`, streams UniRef50 once with `seqkit grep`, stages `ecm33_members.faa` + sha256. |
| `pr6_ecm33_build.py` | stdlib; `manifest` (training set of a unit after the independence rule), `check-alignment` (manifest and a pre-registered occupancy check), `cutoff` (false-hit side only; refuses without a frozen freeze). |
| `pr6_ecm33_build.sbatch` | HPCC array job (DRAFT): per unit MAFFT L-INS-i, `hmmbuild`, false-hit searches, cutoffs, then the held-out searches and the comparators. |
| `pr6_ecm33_report.py` | stdlib; the held-out report; refuses without the frozen freeze and a matching sha256 in every cutoff file. |
| `pr6_ecm33_freeze_PROPOSED.json` | the PROPOSED freeze (`"frozen": false`): the numbers and the four owner decisions it needs. |
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

0. **The freeze commit (owner).** Copy `pr6_ecm33_freeze_PROPOSED.json` to `pr6_ecm33_freeze.json`, decide the four items below, set `"frozen": true`, commit. Nothing
   is built before it: the build job stops if the file is missing or not frozen.
1. Update the HPCC clone to the merge commit of the PR that holds these scripts and the freeze.
2. **Extract** (one job; the file of the round 2 record, 28,781,193,344 bytes, is checked by size, not hashed again):
   ```
   sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<sha>,OUT=$W/ecm33_extract \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_ecm33_extract.sbatch
   ```
   Resources are an UNMEASURED guess (batch, 4 cores, 8 GB, 2 h): one pass of `seqkit` over 28.8 GB should take minutes to a few tens of minutes.
3. **Build** (array: 11 units + `all` = tasks 0 to 11; at most 4 at a time):
   ```
   sbatch --array=0-11%4 --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<sha>,EXTRACT=$W/ecm33_extract,OUT=$W/ecm33_build,PFAM_RELEASE=38.2 \
     -o $W/logs/%x.%A_%a.log -e $W/logs/%x.%A_%a.log $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_ecm33_build.sbatch
   ```
   Inside each task the order is fixed: MMseqs2 clustering of the members (30 percent identity, 50 percent coverage) -> manifest -> MAFFT -> alignment checks ->
   `hmmbuild` -> searches over the four kept proteomes and the 161 shared negatives -> the two cutoffs -> **only then** the held-out searches. L-INS-i on up to
   480 sequences is slow: the 12 h request is a guess.
4. **Report** (local, after the results are copied back):
   `python3 -I analysis/wall_fusion_truth/pr6_ecm33_report.py --results <ecm33_build> --freeze analysis/wall_fusion_truth/pr6_ecm33_freeze.json --units <members/units.tsv> --out <dir>`.
5. The results go into their own record PR (tables and provenance, no sequences).

## What the freeze needs from the owner (four decisions, each with a proposal in the PROPOSED file)

1. **P38248.** The *S. cerevisiae* Ecm33 is one of the 161 shared hard negatives (stratum `gpi_wall_enzyme`) **and** a true member of the family: a good custom
   model should hit it. Proposal: leave it out of the hard-negative rate and report it on its own line. Alternative: keep it, which lifts the cutoff above the
   score of a real Ecm33. (The round 2 Ecm33 check already showed PF12454 scores 0.6 on it, far under GA 23.2.)
2. **Independence mode.** Spec 4.4 excludes every training member that clusters with a held-out member at 30 percent identity and 50 percent coverage
   (`cluster30`). For a conserved family such clusters can span several clades, so the rule can empty a unit's training set (the manifest then says `not_buildable`
   and the unit is reported, not built). Alternative `none` removes only the held-out members. Proposal: `cluster30`, as in the spec, and read the
   count of excluded neighbours per unit before judging any result.
3. **min_u.** u is the number of held-out members that every comparator misses; a unit is testable when u >= min_u and the custom model must recover at least
   ceil(u/2) of them (spec 4.6; the hydrophobin run had u = 0). Proposal: 3.
4. **The numbers.** Cutoff floor 0.0 and hard-negative rate 0.05 (the hydrophobin precedent), extra calls 5 per 10,000, the alignment check (at least 100 columns
   with a residue in at least half of the sequences, median coverage of those columns at least 0.6; a stop condition, the hydrophobin cysteine check does not apply),
   the unit size (at least 5 members). All are proposals.

## Checks and expected outputs

- `extract`: `members/` (members.tsv, excluded.tsv, units.tsv, member_ids.txt, summary.txt), `ecm33_members.faa` (+ sha256), `provenance.txt`. The job fails if
  UniRef50 is not the file of the round 2 record, if any member id is not found, or if a header is duplicated.
- `build`, per unit: `manifest.json` (held-out ids, training ids with clade and 30 percent cluster, the excluded neighbours), `train.aln.fa`, `model.hmm`, `cutoff.json`
  and `relaxed_cutoff.json` (each with the freeze sha256 and `positives_used: false`), the held-out hit tables, `provenance.txt`. A unit that is not buildable keeps only its manifest.
- `report`: `recall_by_unit.tsv` (custom, Pfam GA, relaxed Pfam and iteration recall; u; whether the custom model recovers ceil(u/2) of u), `extra_calls.tsv`
  (extra calls against the freeze limit), `hard_negative_calls.tsv` (rates against the frozen rate), `summary.txt`.

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
