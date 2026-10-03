# Handoff: migrating adhesionPred to run natively on HPCC

> **2026-10-02:** see `docs/HANDOFF-2026-10-02.md` for the current state and worklist. This
> document is kept for its environment notes, standing rules and known traps. Its migration task
> is partly done (`config/site.yaml` and `analysis/_common/paths.*` exist).

**Written 2026-09-28** for the agent taking over on UCR HPCC. Read this first, then
`docs/TOOL-ARCHITECTURE.md`.

> **Why this exists.** The prior session developed from a laptop, pushing scripts to HPCC over
> SSH and running them in a scratch directory (`/bigdata/stajichlab/jstajich/projects/adhesionPred_review`).
> That worked but left hardcoded paths and a dependence on a scratch dir that should not
> survive. **Your primary task is the migration.** A second, scientific task is queued behind it.

---

## 1. Current state

| | |
|---|---|
| repo | `github.com/stajichlab/adhesionPred` |
| `main` | PR #18 merged 2026-09-28 (25 commits) — all analysis and reports are on `main` |
| open branch | `docs-spherule-refs` (one commit, RNA-seq citations) — needs merging, trivial |
| working tree | clean |

**Read before changing anything:**
- `docs/TOOL-ARCHITECTURE.md` — the framework. Different predictors for different intentions; explains why there is no single "adhesin predictor".
- `docs/reports/2026-09-27-coccidioides-antigen-findings.md` — antigen results, shareable with collaborators.
- `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md` — class-2a curation, the most recent work.
- `docs/model-review/STATUS.md` — what the tools can and cannot do, per clade.

---

## 2. PRIMARY TASK — the migration

### 2.1 What is wrong now

**16 script files hardcode `/bigdata/stajichlab/...` paths.** Two also assume a scratch
directory that is not part of the repo.

Files needing edits:
```
analysis/adhesion_properties/01_build_protein_universe.py
analysis/adhesion_properties/db.py
analysis/adhesion_properties/run_task7_seq_props.sbatch
analysis/chytrid_batrach/01_cbm18_vs_calls.py
analysis/chytrid_batrach/02_cbm18_architecture.py
analysis/chytrid_batrach/03_domain_enrichment.py
analysis/chytrid_batrach/04_vwd_profile.py
analysis/chytrid_batrach/05_known_virulence_families.py
analysis/cocci_antigens/01_build_inputs.sh          <- also assumes scratch dir
analysis/cocci_antigens/02_score_antigens.py
analysis/cocci_antigens/03_run.sh                   <- also assumes scratch dir
analysis/cocci_antigens/04_annotate_candidates.py
analysis/cocci_antigens/05_secretion_filter.py
analysis/cocci_antigens/06_tier_candidates.py
analysis/cocci_repeats/01_signalp.sh
analysis/kingdom_survey/01_build_species_table.py
```

### 2.2 The distinct paths to centralise

| key | path |
|---|---|
| `fungi5k_input` | `/bigdata/stajichlab/shared/projects/Fungi_5k/input` |
| `fungi5k_duckdb` | `/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb` (278 GB, read-only) |
| `fungi5k_samples` | `/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv` |
| `cocci_pangenome` | `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome` |
| `cocci_longread` | `/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc` |
| `spherule_rnaseq` | `/bigdata/stajichlab/jstajich/projects/Coccidioides_UCSD_SpheruleMycelium/reports` |
| `rhodotorula_biofilm` | `/bigdata/stajichlab/shared/projects/Rhodotorula/Biofilm/00_input` |
| `workdir` | **decide** — outputs and intermediates; must NOT be `adhesionPred_review` |

### 2.3 Suggested approach

1. `config/site.yaml` (or `site.sh` if you prefer shell-native) holding the table above.
2. A tiny loader — `analysis/_common/paths.py` and a shell equivalent — that reads it and
   exposes the keys. Keep it dependency-light; `yaml` is fine, or use a flat `KEY=value` file
   readable from both Python and bash if you would rather avoid the dependency.
3. Replace hardcoded strings file by file. Run each script after editing — most are fast.
4. Delete the `adhesionPred_review` dependence in `cocci_antigens/01_build_inputs.sh` and
   `03_run.sh`; they should work from a repo clone with `$WORKDIR` from config.
5. Settle **one output rule**. Suggested: large intermediates → `$WORKDIR`; small summary
   tables → tracked in the repo next to the report that cites them. Current state is
   inconsistent and that is worth fixing while you are in here.

### 2.4 Do NOT

- **Do not reorganise `analysis/`.** Eight project dirs with numbered scripts and a report
  each is working. Reorganising for its own sake will only break provenance links in the docs.
- **Do not rewrite the reports** in `docs/reports/`. They are shareable deliverables with
  collaborators; update them additively if results change.
- **Do not regenerate** `data/curated/*.tsv` casually — they query live APIs and row counts
  drift. Rebuild only when intentionally refreshing, and note it in the commit.

---

## 3. Environment notes — hard-won, will save you time

| issue | resolution |
|---|---|
| **SLURM scripts need `#!/bin/bash -l`** | without `-l` the module system is not initialised and `module load` silently fails, then the binary is "command not found". This cost a wasted job. |
| **SignalP 6** | `module load signalp/6-gpu`; binary is **`signalp6`** (not `signalp`). Runs on `short_gpu` with `--gres=gpu:1 --constraint=gpu_latest`, ~285 sequences/s. A container also exists in `$NXF_SINGULARITY_CACHE`. |
| **`short` partition caps at 2 h** | a 4 h request is rejected outright with "Requested time limit is invalid". Use `short_gpu` (also 2 h) or `epyc`/`gpu` for longer. |
| **MMseqs2 crashes on some `short` nodes** | "Illegal instruction (core dumped)" on node c18 — the module is built for newer CPU instructions. **Use `-p epyc` for MMseqs2.** |
| **`AssocGrpMemLimit` pending** | the lab's group memory allocation gets saturated by other jobs. Request less memory; most of these jobs need 8–12 GB, not 32. |
| **login-node memory cap** | loading PLM weights on the login node gets OOM-killed (exit 137). Do model loading inside `srun`, and download-only on the login node. |
| **GPU partitions** | `short_gpu` has ada6000/a100/h100/blackwell. `exfab` (gpu12, 2× ada6000) is high-priority if you have access; the shared `gpu` partition had ~10 h fairshare waits. |
| **DuckDB memory** | set `PRAGMA memory_limit='3GB'` explicitly — it auto-detects host RAM and ignores the cgroup cap, then gets OOM-killed. |
| **Python env** | a working venv with torch/esm/transformers/duckdb/scipy/sklearn is at `adhesionPred_review/analysis/model_review/pilot/.venv`. Rebuild it properly in the new location — `esm` pulls a mismatched `torchvision` from PyPI, so install torch+torchvision together from the cu126 index, and `httpx` separately (the `esm` package imports it but does not declare it). |

---

## 4. QUEUED SCIENTIFIC TASK — the 2a retrain test

Do this after the migration. It is the experiment that says whether the curation strategy works.

**Background.** Class 2a (repeat/avidity surface proteins) works in Saccharomycotina but misses
SOWgp. `analysis/model_review/repeat_structure_transfer.py` currently recovers **3/5** of the
cross-clade repeat proteins. The hypothesis: the training positives are long and Ser/Thr-rich
(FLO11 1367 aa, ALS1 1260 aa) while SOWgp is short and Pro/Cys-rich (324–422 aa), so length and
composition are confounded with clade.

**The test.** `analysis/cocci_repeats/class2a_candidates.tsv` holds 41 curated candidates,
**20 of them Pro/Cys-rich** — exactly the architecture the training set lacks. Add them as
training positives and re-run the transfer test.

- **Success looks like:** recovery improves above 3/5, and specifically SOWgp58/SOWgp66 (currently 0.334 / 0.450) cross 0.5.
- **A null result is equally publishable** and should be reported, not buried — it would mean the gap is not about training composition.
- **Caution:** the 41 candidates are computational predictions with no functional validation. Adding them as positives risks training on noise. Consider holding out SOWgp itself, and report performance with and without the new positives.

Ranked after that (detail in the reports): read-depth at the SOWgp locus to separate real
absence from collapsed arrays; identify the period-17 Ser/Thr family (10 members, 5 copies in
Cpos1038 vs 1 in Silveira, uncharacterised); add NetGPI to the surface filter, which is
currently SignalP-only and therefore missing GPI-anchored proteins.

**Parked, awaiting data:** *C. posadasii* host-outcome phenotypes, to be tested against the
1,580 accessory / 286 copy-variable orthogroups. And class 2b (small invasins) — scoped as
structure + Foldseek homology search, **not** classification, because there are only 2–3 known
examples and no classifier can be trained on that.

---

## 5. Standing rules in this project

These were adopted after specific failures. Please keep them.

1. **Name a tool after what it detects, not what you want it to find.** The shipped classifier
   is a tandem-repeat surface protein detector, not an adhesin predictor.
2. **State clade scope.** Adhesin families are phylogenetically disjoint.
3. **Hold every predictor to an acceptance test on known positives, and let it fail loudly.**
   The *Coccidioides* scorer prints `NOT CALIBRATED` at 3/4 and that warning is deliberate.
4. **Keep score axes separate.** Merging antigenicity with specificity hid that the clinical CF
   antigen is both immunogenic and cross-reactive.
5. **Absence is not established by a missing call.** Collapsed repeat arrays, gene-model failure
   and an arbitrary threshold all look identical. Say "not assessable".
6. **Say when the right tool is not ML.** Stage 1 is a SignalP problem; hydrophobins are an HMM
   problem; antigen specificity is comparative genomics.

## 6. Known traps that already cost time

- UniProt's *C. posadasii* reference proteome is strain **C735 ΔSOWgp** — a SOWgp **deletion** strain.
- Pfam **PF10528 is GLEYA, not Flo11** (Flo11 is PF10182).
- **CBM18 is `Chitin_bind_1` (PF00187)** and is absent from the CAZy overview table — querying CAZy for "CBM18" silently returns nothing.
- Fungi_5k SignalP coverage is thin (~4% of *C. immitis* RS proteins), so "no signal peptide" often means unannotated, not negative. **SOWgp itself has no Fungi_5k match.**
- Pfam HMM `NAME` is not always the label you would use — keying on the wrong one silently returns zero.
- `hash()` in Python is randomised per process; do not use it for CV fold assignment.
