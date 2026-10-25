# Aggregate report with Pfam, TMHMM and repeat02 modules: 10 proteomes (2026-10-26)

Second aggregate report, after adding the `tm`, `repeat02` and `pfam_*` modules to the 10 proteome runs of
`docs/reports/data/aga1_aggregate_2026-10-26/`. Only 3 calls stay NA (`cocci_specificity_rank_top15`, `iuis_allergen_similarity`,
`iuis_allergen_homolog`): the allergen BLAST and antigen lookup have not been run on these proteomes.

How it was made (inputs are outside the repository; no sequences are stored):

1. Tool output: the HPCC proteome array (`analysis/aga1_class/proteome_array.sbatch`, TOOL = pfam_family, tmhmm; Pfam 38.2 `--cut_ga` over the models of `data/sorting_hat/family_table.tsv`, TMHMM 2.0c `-short`) was copied read-only into the local results folders.
2. `python3 -I analysis/sorting_hat_run/aga1_domain_modules.py --runs-root <runs> --taxdump nodes.dmp --taxa analysis/sorting_hat_run/aga1_acceptance_taxa.tsv --results <results> ...` merges the chunks per species, runs `02_repeat_profile.py` on the proteome, writes the modules and runs the engine into `<species>/out_domains`.
3. `cellsurface_sorting_hat_aggregate --manifest manifest.tsv --taxdump nodes.dmp --names names.dmp --rank class --out <new folder>` (manifest: the ten `out_domains` folders, labels as in the first report).

`categories.yaml` sha256 is unchanged (`d12d160697632ff5f92d673205c691a54714d2e6b13eec816ee10b6a515fe50b`), so the smoke status files of the first report still apply (4 of 200 proteome-call cells have a measured status).

| file | sha256 |
|---|---|
| counts_by_group.tsv | 9e74b38f6ec0cf90308236772a9e74257525477dead4fbf5d5b38941da6076dd |
| counts_by_proteome.tsv | 4a8bb5541e290780e66ed7108ac9141e83fc68dd4f379a8ab4d270caf8f1e7b6 |
| proteomes.tsv | 3e7bcc1b345aa6a245d73766a6d38def822991891c77e4a1e92741c1a749a686 |
| report.md | 55897782169d039a0cb1ccd4ec5c8e7a8006a8334fb9e849f136415200111e99 |

## Reading the numbers

* `pfam_allergen` has `run_state: unavailable` in every proteome (no active allergen family in the table); `hydrophobin_domain` is called only in the Eurotiomycetes (9 of 19,359 proteins), none elsewhere.
* The family-table calls and `tandem_repeat_protein` have status `unvalidated` in these proteomes (marked `*`). Counts are counts of calls, not prevalences.
* `repeat02` was run locally with the project detector at its default thresholds (coverage 0.25, copies 2.5 in the wrapper).
