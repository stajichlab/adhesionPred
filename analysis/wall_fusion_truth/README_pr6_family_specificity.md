# PR 6: specificity tests of the candidate Pfam families (DRAFT, not run on HPCC)

Owner decisions of 2026-10-22: **1.** the candidate list is the union of a Pfam hit that covers at least 0.5 of the sequence and a
best hit with an E-value of 1e-30 or better (both limits inclusive): 26 clusters, 30 families
(`docs/reports/data/pr6_family_candidates_2026-10-22/`). **2.** Every candidate family gets a specificity test BEFORE any panel
addition: nothing enters the active panel without a `family_measure.tsv` row and the owner's sign-off, and `narrow_ok` stays `no` until
the owner signs it.

**This work changes no tier, no panel, no status, no `narrow_ok`, and writes neither `data/sorting_hat/family_measure.tsv` nor
`family_table.tsv`.** The job writes report-only tables. The converter writes a PROPOSED measurement file elsewhere (or prints a dry run)
and refuses the two committed paths (tests: `tests/wall_fusion_truth/test_pr6_family_specificity.py`).

## Files

| File | Purpose |
|---|---|
| `pr6_family_specificity.py` | stdlib. `candidates` (the rule, offline), `summarize` (the tables from the `hmmsearch` output), `convert` (dry-run proposal of `family_measure.tsv` rows). |
| `pr6_family_specificity.sbatch` | HPCC job (DRAFT, NOT SUBMITTED): pin check, inputs check, `hmmfetch` of the 30 models, `hmmsearch --cut_ga` over six FASTA files, `summarize`. |
| `docs/reports/data/pr6_family_candidates_2026-10-22/` | the committed candidate tables. The job makes them again and stops if one differs. |

## What is searched, and what is counted

- Models: the 30 candidate families, fetched by NAME from Pfam 38.2 (`hmmfetch -f` looks up by name). A candidate without a GA line stops the
  job with a message (every Pfam-A family has one; this is a guard, not a silent skip).
- **Own members (D19).** The family's own evidenced members are all members of the candidate clusters in which the family has a qualifying hit
  (`family_members.tsv`: 73 rows for the 30 families). The own-member count is the denominator. **The own-member recall is circular**: the family was
  selected because it hits these members. It is a count, never a sensitivity, and the `leakage` field says `tuned_on_truth`. In the proposed
  measurement rows `n_pos_unknown` equals `n_pos`: nobody has checked whether the members are in the family's Pfam seed.
- Other Kiss positives that hit the family are counted apart (per slice), outside the denominator: a descriptive cross-reaction.
- **Negatives:** the shared hard-negative set (`data/controls/hard-negatives-shared/`, 137 N1 and 24 N2 rows; 91 clusters in all), scored under the plain
  accession from the committed FASTA. The counts follow the 2026-10-09 family test: N1 rows and N1 clusters are the denominators, N2 hits a separate
  count. The `#neg` copies exist only in the Kiss clustering FASTA and are not involved; the six accessions that are both a Kiss control and a shared
  negative are a positive in one table and a negative in the other (D20).
- **Kept proteomes (decision 9b, one per species):** Afum_Af293_UniProt (9,647), Calb_SC5314 (6,212), Cimm_RS (9,910), Scer_S288C (6,722); 32,491
  proteins. Hits per proteome and per 10,000 proteins, and every hit with its FASTA header text (`proteome_hits.tsv`). Descriptive: the proteomes
  are not labelled. The excluded proteomes (Afum_Af293_Fungi5k, Bder_ER3 with internal-ID headers; Scer_S288C_nodubious, Afum_A1163,
  Afum_W72310, a second proteome of a species) are neither read nor named in the job.
- No interval, no sensitivity and no specificity is computed: the tier comes from `cellsurface_sorting_hat.family_tier`, never from this script.

## Run order (all commands NOT run; each HPCC step waits for the owner's explicit go-ahead)

0. **Locally, offline (safe, already done for the committed tables):**
   `python3 -I analysis/wall_fusion_truth/pr6_family_specificity.py candidates --out /new/empty/dir` and compare with
   `docs/reports/data/pr6_family_candidates_2026-10-22/` (the tests do this).
1. **HPCC, update the clone** (a write): `git -C $B/adhesionPred-kiss53 pull --ff-only origin main`; it must print the merge commit of this PR and the
   tree must be clean (`status --short` prints nothing).
2. **Submit** (a write; no new staging is expected: `kiss_sequences.faa` with its hash is already at `$W/in`, the rest is in the clone and the
   proteomes are read in place):
   ```
   B=/bigdata/stajichlab/jstajich/projects
   W=$B/adhesionPred/_workdir/kiss_cluster_split
   sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<merge commit>,IN=$W/in,OUT=$W/family_specificity_out,PFAM_RELEASE=38.2 \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log \
     $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_family_specificity.sbatch
   ```
   Partition `short`, 4 cores, 8 GB, 1 hour; writes only to a new `$OUT`.
3. **Read the outputs** (read only): `scp` the small tables to the local machine; the domain tables stay on HPCC.
4. **Record PR** (the results as a run record and a short report, with no recommendation beyond the tables). Only after that, and only with the
   owner's sign-off per family: `convert` (dry run first):
   `PYTHONPATH=src python3 analysis/wall_fusion_truth/pr6_family_specificity.py convert --specificity <family_specificity.tsv> --source <report path> --date <date>`
   prints the tier that the derivation gives each proposed row (`smoke` or `untested`; `descriptive_shared_sets` is capped at `smoke`). Writing rows
   into `data/sorting_hat/family_measure.tsv` and re-deriving `specificity_tier` is a separate owner-approved PR (D19, issue #69), not part of this one.

## Expected outputs (`$OUT`)

`family_specificity.tsv` (one row per family: own members and hits, other-positive hits per slice, N1 and N2 hits and clusters, hits per proteome and
per 10,000, kept total), `positives_hits.tsv`, `negatives_hits.tsv`, `proteome_hits.tsv`, `provenance.tsv`, `provenance.txt`, and the six `domtbl_*.txt`.
The tables hold accessions, counts and FASTA header text; no sequence.

## Checks the job makes

Pin of the code, the step 6.0 record, the candidate tables, the shared negatives and `inputs.tsv`; candidate tables made again equal the committed ones;
132 positives found in the staged FASTA; 161 negatives named as `controls.tsv`; each kept proteome exists and holds the `n_proteins` of `inputs.tsv`;
Pfam release name `*Pfam38.2` and the sha256 against the earlier jobs (a difference is reported, not fatal); 30 models fetched with the expected
accessions and names; every `hmmsearch` ends with `# [ok]`; `summarize` refuses hits on a sequence that is not in its FASTA and refuses a non-empty output.

## Not confirmed (the author of this draft did not touch HPCC)

- The four proteome files at `$PROTEOME_DIR/<name>.faa` (default `.../_workdir/sorting_hat/`, the `OLD` directory of `lib.sh`): these
  files were read earlier in the PR 6 work (headers and protein counts equal `inputs.tsv`), but this draft did not re-check them. The job stops if one is missing or its count differs.
- `hmmfetch -f` by NAME for the 30 models (the same method as the Task 5.4 job, which fetched 39 models; no SSI index is next to `Pfam-A.hmm`). The job
  stops with `hmmfetch found N of 30` if a name is missing.
- Clone contents: the clone must hold the merges of the PRs that carry the step 6.0 record and this PR.

## Observations to carry into the owner's decision (from the committed tables, not from a run)

- `PF11765` (Hyphal_reg_CWP) is already a row of `family_table.tsv` (an inactive family); `convert` flags rows that already have a measurement.
- Several candidates are sub-domains of one protein: five Mok13 families (KS0115), PMT and PMT_4TMC (KS0016, KS0035), Glyco_hydro_63 and 63N (KS0114),
  Glyco_transf_20 and Trehalose_PPase (KS0119). Their hits are not independent.
- A Pfam hit is not specificity; this job is the test that says whether a family stays quiet on negatives and proteomes.
