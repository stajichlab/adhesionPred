# Step 1 Phase B: Features and Embeddings Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the tracked SignalP 6 and PredGPI feature table for all truth species and *C. immitis* RS (D2), the resumable ESM-2 8M and 35M embedding job with the C-terminal window variants (D3), and the throughput pilot J0, so that a controller can submit J0, J1 and J2 after an independent review.

**Architecture:** A stdlib-only step (05) collects every input set, dedupes the sequences by `seq_sha256` and gives each unique sequence a fixed matrix row. J0 measures the real embedding throughput. A stdlib step (06) derives the chunk size and the J2 job count from that measurement. J1 (one GPU job) runs SignalP 6 and PredGPI on the unique sequences; J2 (GPU, one or more array tasks) embeds deterministic, length-sorted chunks through `surface_glyco.embeddings.get_esm_embeddings` and stores each chunk with a hash, so a rerun skips finished chunks. A stdlib step (07) merges the tool outputs into `features.tsv.gz` with the truth labels and the embedding row numbers; `assemble_embeddings.py` checks the chunks and builds one float32 matrix per model and window.

**Tech Stack:** Python 3.12 standard library (`/usr/bin/python3.12`) for every module in `analysis/step1_compare/`; the conda env `adhesionPred` (Python 3.14.2, torch 2.10.0+cu128, fair-esm 2.0.0, numpy 2.4.2) for `analysis/step1_compare/jobs/*.py`; module `signalp/6-gpu` (SignalP 6.0i GPU build, its own env with torch 2.8.0+cu128); module `predgpi/202001` (Python 3.9.18 with numpy); SLURM partition `exfab` (node gpu12, 2 x ada6000); bash; pytest; ruff 0.3.5 through pre-commit.

**Spec:** `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` (sections 5, 6, 7, 8, decision Q10 in 10, deliverables D2 and D3 in 11). Phase A plan and code: `docs/superpowers/plans/2026-10-01-step1-truth-set.md`, `analysis/step1_compare/` (README.md, COLUMNS.md).

**Scope.** This plan covers D2, D3 and J0 only. It does not cover training, the evaluation harness, gates or rule fitting. It submits no job. The section "Run (controller, after review)" lists the exact `sbatch` commands and acceptance checks.

## Facts measured for this plan (2026-10-01)

All numbers below were computed on 2026-10-01 from the files named. Numbers marked "assumption" are not measured.

| Fact | Value | Source |
|---|---|---|
| Phase A production run | finished 02:32:23 PDT, every script rc=0, all four run logs say `all_sources: true` | `$STEP1_WORKDIR/logs/phaseA_production.log`, `extract_log.json`, `sequence_run.json`, `d8_run.json`, `keyword_tier_run.json` |
| `truth_sequences.tsv.gz` | 50,844 rows (51,106 truth genes minus 262 unmatched); 41,293 unique sequences; 22,247,415 residues in unique sequences; 3,920 unique sequences longer than 1,022 aa; longest 8,515 aa; shortest 4 aa; only non-standard residue: X (418 occurrences); 0 stored hashes differ from a recomputed hash | prototype count script |
| Phase B universe (05 on the production inputs, all 11 sets in `sequence_sets.tsv`) | 123,567 members; 69,941 unique sequences; 34,524,441 residues; 5,541 unique sequences longer than 1,022 aa; no character outside the ESM-2 alphabet | prototype run of `05_prepare_sequences.py` |
| Residues to embed per model | 37,860,229 (32,197,327 N-terminal window residues + 5,541 x 1,022 C-terminal window residues) | prototype run of `06_plan_embedding.py --rate 29500` |
| Proteome FASTA records (unique sequences) | S288C 6,722 (6,638); *C. albicans* 6,250 (6,192); *S. pombe* PomBase 5,126 (5,078); Af293 9,647; *A. nidulans* 10,561; H99 7,427; JEC21 6,740 (6,739); *U. maydis* 6,805; *C. immitis* RS 9,910 (9,875); keyword cache 3,573 (3,509) | prototype count script |
| ESM-2 weights | `esm2_t6_8M_UR50D.pt` and `esm2_t12_35M_UR50D.pt` are in `~/.cache/torch/hub/checkpoints/` | `ls` |
| conda env torch | 2.10.0+cu128; `torch.cuda.is_available()` is False on the CPU node c01. Not run on a GPU node yet | `python -c` on c01 |
| CPU embedding rate (c01, 2 cores, 20 real truth sequences, batch 8) | 8M: 364.3 residues/s; 35M: 133.4 residues/s. At these rates 37.86 M residues take 28.9 h (8M) and 78.8 h (35M). J2 needs a GPU | `jobs/throughput_pilot.py --device cpu --n 20` |
| GPU embedding rate for 8M and 35M | not measured. Review 9.1 gives ESM-2 150M at 60 proteins/s and 29,500 residues/s on one RTX 6000 Ada, but with bf16, HuggingFace transformers and token-budget batching, not the fp32 fair-esm code path of `get_esm_embeddings`. Used below only as an assumption for the worked example | `docs/model-review/2026-09-27-review-and-framework-plan.md` 9.1 |
| SignalP 6 GPU build | `module load signalp/6-gpu` works only in a login shell (`bash -l`), because the module runs `conda activate`. On the CPU node c01 it stops with `RuntimeError: The model checkpoint is configured to run on GPU, but no CUDA device was found`. The CPU build has no weights (verified 2026-09-30, `analysis/cocci_repeats/01_signalp.sh`) | run on c01 |
| SignalP 6 output format | `prediction_results.txt`: TAB separated, two `#` lines, header `# ID<TAB>Prediction<TAB>OTHER<TAB>SP(Sec/SPI)<TAB>CS Position`; SP rows end with `CS pos: 24-25. Pr: 0.5487`, OTHER rows with an empty field. `output.gff3`: one `signal_peptide` line per SP protein, start 1, end = first number of `CS pos`. The ID column is the whole FASTA header line, so J1 uses the 64-character `seq_sha256` as the only header token | `analysis/cocci_repeats/signalp/CimmitisRS_FungiDB/` (job 29280458) |
| SignalP 6 GPU rate | 82,326 proteins in 6 min 16 s wall time (7 proteomes, model load included; 217 to 244 proteins/s in the progress bar) on gpu12 with 8 cores and one ada6000 | `logs/signalp.29280458.log`, `sacct -j 29280458` |
| PredGPI CLI | `predgpi.py -f FASTA -o OUT -m {json,gff3}`; writes one line per protein: `GPI-anchor` (omega site, score 1.0 if FPR <= 0.0015, 0.70 if <= 0.005, 0.55 if <= 0.01) or `Chain`. No continuous score. Sequences of 40 aa or less are never GPI. Duplicate FASTA IDs collapse (a dict) | `$PREDGPI_HOME/predgpi.py`, `predgpilib/utils.py` |
| PredGPI rate | 1,000 S288C proteins (432,782 residues) in 32.5 s on one core of c01 (31 proteins/s) | `time predgpi.py` |
| PredGPI wrapper | `jobs/predgpi_scores.py` gives the same class, omega site and score as the CLI for all 1,000 proteins (0 mismatches) and adds `fpr` and `svm`; same run time (32.3 s) | prototype comparison |
| Partitions | `exfab`: node gpu12, `gres/gpu:ada6000=2`, MaxTime 30 days, DefaultTime 7 days, account `stajichlab` works (job 29190082 runs there). `short_gpu`: MaxTime 2 h. NetGPI: no module | `scontrol show partition`, `sacctmgr`, `module avail` |

## Global Constraints

- Analysis modules in `analysis/step1_compare/` (not `jobs/`): Python 3.12 standard library only (`PY=/usr/bin/python3.12`); `test_every_module_imports_only_stdlib_or_local` must stay green.
- `analysis/step1_compare/jobs/*.py`: the conda env Python `ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python` with `PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/jobs`. `jobs/predgpi_scores.py` runs under the PredGPI module Python 3.9 and must stay Python 3.9 compatible.
- Pooling is never re-implemented: every embedding comes from `surface_glyco.embeddings.get_esm_embeddings` (layer 6, residue mean, padding excluded, `MAX_RESIDUES = 1022`).
- ruff 0.3.5 through pre-commit, line length 100. Run `pre-commit run --all-files` before every commit; re-stage if a hook reformats.
- STOP contract (Phase A README): a script that cannot continue prints `STOP: <reason>` to stderr, exits with code 2 and writes no output. Outputs are written with `runinfo.atomic_write_all` (temp name, then `os.replace`). Run logs record input hashes, `all_sources`, git commit, Python version and arguments, and no time stamp.
- Work directory from `PROJ_ROOT` and `STEP1_WORKDIR` (`paths.workdir()`); Phase B outputs go to `$STEP1_WORKDIR/phaseb/`. Data never go into the repository.
- No script uses `BASH_SOURCE` (not even in a comment: `test_no_shell_script_uses_bash_source` scans the text).
- SLURM scripts: `#!/bin/bash -l`, `set -euo pipefail`, `${SCRATCH:?}` node-local work, results copied to `/bigdata` before exit, `PROJ_ROOT` and `STEP1_WORKDIR` required from the environment (`: "${PROJ_ROOT:?...}"`), explicit `#SBATCH --time`.
- Job size: about 1 to 1.5 h of real run time per job, derived from measured throughput; never split work into many jobs of minutes.
- Compress large text outputs (`gzip -n` or the byte-stable `truth_table.write_tsv` for `.gz`); scripts read gzip input transparently (`gaf.open_text`, `truth_table.read_tsv`).
- Prose in code, docs and logs: Simplified Technical English; no claim without a measured number or a test.

## Review Focus

1. **A tool silently skips sequences** (SignalP drops a record; PredGPI collapses two records with one ID). Expected: 07 stops and names the counts; with `--allow-missing-calls` the gap is in `feature_coverage.tsv`. Pinned by `test_missing_predgpi_call_stops_unless_allowed` (Task 8) and by the duplicate-id check in `jobs/predgpi_scores.py` (Task 7).
2. **Stale outputs after 05 is re-run** (the unique set changes, old J1 parts or J2 chunks remain). Expected: 07 stops on ids that are not unique sequences; J2 recomputes a chunk whose JSON names other members; edited chunk members are refused. Pinned by `test_stale_tool_output_stops` (Task 8), `test_read_plan_detects_edited_members` (Task 3), `test_corrupted_chunk_is_recomputed_with_the_same_hash` (Task 4).
3. **A killed or failed job leaves half-written files under final names.** Expected: no done marker for an unfinished part or chunk; the rerun resumes and gives the same arrays. Pinned by `test_kill_and_resume_gives_the_same_arrays` (Task 4) and `test_j1_signalp_failure_leaves_no_done_marker` (Task 9).
4. **A character that `sanitize_sequence` would delete** (`-`, `.`) shifts the window boundary, so the C-terminal window is not the last 1,022 model residues. Expected: 05 stops. Pinned by `test_prepare_stops_on_a_non_esm_character` (Task 2).
5. **A GPU job without a visible GPU falls back to CPU** (28.9 h and 78.8 h at the measured CPU rates). Expected: J0, the diff harness and J2 stop with exit 2. Pinned by `test_cuda_request_without_gpu_stops` (Task 4) and `test_gpu_cpu_diff_stops_without_a_gpu` (Task 6).

---

## Design

### Data model (later phases join on these keys)

- **Member**: one protein of one input set, key `(set_id, source_id, gene_id)`. For the `truth` set, `(source_id, gene_id)` is the key of `truth_set.tsv.gz` and `truth_set_triaged.tsv.gz`.
- **Unique sequence**: key `seq_sha256` (`seqhash.seq_sha256`, the Phase A hash). Sorted by `seq_sha256`; `row` = position (0-based). Each unique sequence is embedded and scored once.
- **Embedding rows**: `emb/<model>.nterm.npy` has one row per unique sequence (`row`). `emb/<model>.cterm.npy` has one row per sequence longer than 1,022 aa (`cterm_row`, numbered in `seq_sha256` order). For shorter sequences the C-terminal window is the whole sequence, so M8-C and M35-C use the `nterm` row; `assemble_embeddings.window_matrix` builds the full C-terminal matrix. This saves 2 x 64.6 M floats of duplicate rows.
- **Features**: `features_unique.tsv.gz` (one row per unique sequence) and `features.tsv.gz` (one row per member, with `label`, `subset`, `stratum`, `d8_class`, `homology_only`, `role` for truth members, the features, `emb_row` and `emb_cterm_row`).

### Input sets (`sequence_sets.tsv`)

D2 asks for features "for all truth species and *C. immitis* RS", and spec section 6 asks for a rule-versus-ML agreement matrix on whole S288C, *C. albicans*, *S. pombe* and *C. immitis* RS proteomes. So the universe is: `truth_sequences.tsv.gz`; all records of `keyword_sequences.fasta.gz` (T-c rows and literature seeds, needed for V-kw training); the eight proteome FASTA files that 00 downloaded; and the *C. immitis* RS FASTA that `analysis/cocci_repeats` uses (path from `config/site.yaml` key `cocci_pangenome`). The whole proteomes add 28,648 unique sequences to the 41,293 truth sequences (69,941 total). This is a design decision of this plan; the owner can drop a set by deleting its row.

### Windows (Q10)

`seqwindow.nterm_window` = first 1,022 residues (what `get_esm_embeddings` keeps); `seqwindow.cterm_window` = last 1,022 residues, only for sequences longer than 1,022 aa. 05 refuses any character outside the ESM-2 alphabet, so `sanitize_sequence` is the identity on every stored sequence and the windows are exactly the model input.

### Chunks and job sizing (J2)

A chunk is the unit of resumability: entries `(window length, seq_sha256)` sorted, filled greedily up to `chunk_residues`. Membership depends only on the unique set and `chunk_residues`, and `members_sha256` records it. A job is a group of chunks; chunk k goes to array task k mod n_jobs. The size parameters come from the J0 result (`chunk_plan.py`):

```
chunk_residues = max(10,000, floor(min_m r_m * 300 s / 10,000) * 10,000)
T_total        = sum over models m of (R / r_m + load_m)      R = 37,860,229 residues
n_jobs         = max(1, ceil(T_total * 1.25 / 4,500 s))
time_minutes   = ceil((T_total / n_jobs * 1.5 + 600) / 60)
```

r_m = residues per second of model m at its best batch size in `j0/throughput.json`; load_m = its load time. 300 s bounds the work lost when a task is killed; 4,500 s (1.25 h) is the middle of the 1 to 1.5 h rule; 1.25 is a safety factor.

**Worked example (assumption, until J0 runs):** r = 29,500 residues/s for both models (the review 9.1 rate of ESM-2 150M on another code path), load 10 s each. Then chunk_residues = 8,850,000; the real unique set gives 5 chunks (nterm 4, cterm 1; prototype dry plan); T_total = 2 x 1,283.4 + 20 = 2,586.8 s; n_jobs = ceil(3,233.5 / 4,500) = 1; time_minutes = 75. If J0 measured 2,000 residues/s for one model alone, the same formula gives 6 jobs (`test_plan_jobs_splits_slow_rates_into_jobs_under_the_target`).

### J1 sizing (one job)

Measured rates: SignalP 6 on gpu12 about 219 proteins/s over a whole job (82,326 in 376 s); PredGPI 31 proteins/s per core. For 69,941 sequences: SignalP about 320 s; PredGPI 2,256 core-seconds, about 282 s per part with 8 parallel parts (assumption: gpu12 cores are as fast as c01 cores). The whole J1 work is about 10 minutes, so one job is the minimum count, and it cannot be sized up to 1 h. It runs on `exfab` with an explicit `--time=1:00:00`. `exfab` and not `short_gpu` (2 h limit, also enough): `analysis/cocci_repeats/01_signalp.sh` records that `short_gpu` scheduled about 4 h out on 2026-09-30 while `exfab` had no pending job. SignalP needs a GPU (CPU build has no weights; GPU build refuses CPU). The FASTA is cut into 8 parts so a rerun resumes per part.

### Storage

Embeddings are binary float32 `.npy` files, not text, so the text compression rule does not apply; compression of float32 embeddings was not measured. Sizes (computed): 8M nterm 89.5 MB, cterm 7.1 MB; 35M nterm 134.3 MB, cterm 10.6 MB; 241.5 MB for the four matrices, and the same again for the chunk files. All text tables are gzip.

## File structure

| File | Responsibility | Task |
|---|---|---|
| `analysis/step1_compare/seqwindow.py` | N- and C-terminal windows, ESM-2 alphabet check | 1 |
| `analysis/step1_compare/sequence_sets.tsv` | the 11 input sets | 2 |
| `analysis/step1_compare/seqsets.py` | members, dedupe, unique rows | 2 |
| `analysis/step1_compare/05_prepare_sequences.py` | Phase B step 1 script | 2 |
| `analysis/step1_compare/runinfo.py` (modify) | two more entries in `PRODUCERS` | 2, 8 |
| `analysis/step1_compare/chunk_plan.py` | chunks, job sizing formula, plan reader | 3 |
| `analysis/step1_compare/06_plan_embedding.py` | chunk plan and job plan script | 3 |
| `analysis/step1_compare/jobs/embed_store.py` | chunk files: save, verify, load | 4 |
| `analysis/step1_compare/jobs/embed_chunks.py` | J2 runner | 4 |
| `analysis/step1_compare/jobs/assemble_embeddings.py` | chunk check, matrices, chunk manifest | 5 |
| `analysis/step1_compare/jobs/throughput_pilot.py` | J0 pilot | 6 |
| `analysis/step1_compare/jobs/gpu_cpu_diff.py` | spec 7 GPU-CPU harness | 6 |
| `analysis/step1_compare/feature_parsers.py` | SignalP and PredGPI parsers, Ser+Thr | 7 |
| `analysis/step1_compare/jobs/predgpi_scores.py` | PredGPI with continuous scores | 7 |
| `analysis/step1_compare/07_build_features.py` | feature table, coverage report | 8 |
| `analysis/step1_compare/jobs/j0_pilot.sh`, `j1_features.sh`, `j2_embed.sh` | SLURM jobs | 9 |
| `analysis/step1_compare/COLUMNS.md`, `README.md` (modify) | output columns, Phase B section | 10 |
| `tests/step1_compare/conftest.py` (modify) | put `jobs/` on `sys.path` | 4 |
| `tests/step1_compare/phaseb_fixture.py` | tiny embedding plan for tests | 4 |
| `tests/step1_compare/fixtures/phaseb/` | tool output fixtures, SignalP stub | 7, 9 |
| `tests/step1_compare/test_*.py` (9 new files) | tests | 1 to 10 |

## How to run the tests

```bash
cd /bigdata/stajichlab/jstajich/projects/adhesionPred-feat    # the worktree of branch step1-features
PY=/usr/bin/python3.12
ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
$PY -m pytest tests/step1_compare -q                          # stdlib run: torch test files skip
module load predgpi/202001                                    # enables the PredGPI wrapper test
PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q       # full run
```

The J1 smoke tests (Task 9) need the HPCC module system (`BASH_FUNC_module%%` in the environment) and `predgpi/202001`; elsewhere they skip.

---

### Task 1: ESM-2 windows

**Files:**
- Create: `analysis/step1_compare/seqwindow.py`
- Test: `tests/step1_compare/test_seqwindow.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `seqwindow.MAX_RESIDUES = 1022`; `seqwindow.ESM_RESIDUES: frozenset[str]`; `seqwindow.WINDOWS = ("nterm", "cterm")`; `bad_characters(sequence: str) -> set[str]`; `needs_cterm(sequence: str) -> bool`; `nterm_window(sequence: str) -> str`; `cterm_window(sequence: str) -> str`; `window(sequence: str, name: str) -> str` (ValueError for an unknown name).

- [ ] **Step 1: Write the failing test**

```python
import pytest
import seqwindow


def _protein(n: int) -> str:
    # A sequence whose residue at position i (0-based) is unambiguous for slicing checks.
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(alphabet[i % 20] for i in range(n))


def test_cterm_window_takes_last_1022():
    seq = "M" * 500 + _protein(1022)
    got = seqwindow.cterm_window(seq)
    assert len(got) == 1022
    assert got == seq[500:]
    assert got[-1] == seq[-1] and got[0] == seq[500]


def test_nterm_window_takes_first_1022():
    seq = _protein(1500)
    assert seqwindow.nterm_window(seq) == seq[:1022]


@pytest.mark.parametrize("n", [1, 40, 1021, 1022])
def test_short_protein_windows_are_the_whole_sequence(n):
    seq = _protein(n)
    assert seqwindow.nterm_window(seq) == seq
    assert seqwindow.cterm_window(seq) == seq
    assert not seqwindow.needs_cterm(seq)


def test_needs_cterm_only_above_the_limit():
    assert not seqwindow.needs_cterm("A" * 1022)
    assert seqwindow.needs_cterm("A" * 1023)


def test_window_dispatch_and_unknown_name():
    seq = _protein(1100)
    assert seqwindow.window(seq, "nterm") == seq[:1022]
    assert seqwindow.window(seq, "cterm") == seq[78:]
    with pytest.raises(ValueError, match="unknown window"):
        seqwindow.window(seq, "middle")


def test_bad_characters_lists_what_esm_would_change():
    assert seqwindow.bad_characters("MKXBUZO") == set()
    assert seqwindow.bad_characters("MK-J.*") == {"-", "J", ".", "*"}


def test_constants_match_the_package():
    # The package imports torch; skip where it is not installed (py3.12 stdlib runs).
    pytest.importorskip("torch")
    pytest.importorskip("esm")
    from surface_glyco import card, embeddings

    assert seqwindow.MAX_RESIDUES == card.MAX_RESIDUES
    assert set(embeddings.ESM_RESIDUES) == set(seqwindow.ESM_RESIDUES)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_seqwindow.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'seqwindow'`.

- [ ] **Step 3: Write the implementation**

```python
"""ESM-2 input windows for Phase B (spec 5, decision Q10). Standard library only.

ESM-2 takes at most MAX_RESIDUES residues (1,024 tokens with BOS and EOS). The N-terminal
window is the first MAX_RESIDUES residues; `surface_glyco.embeddings.get_esm_embeddings` cuts
there too. The C-terminal window is the last MAX_RESIDUES residues. It exists only for proteins
longer than MAX_RESIDUES; for shorter proteins both windows are the whole sequence.

Every sequence must hold only ESM_RESIDUES. Then `sanitize_sequence` in the package changes
nothing, and a window cut here is exactly the input the model sees.
"""

MAX_RESIDUES = 1022  # same value as surface_glyco.card.MAX_RESIDUES
ESM_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYXBUZO")  # surface_glyco.embeddings.ESM_RESIDUES
WINDOWS = ("nterm", "cterm")


def bad_characters(sequence: str) -> set[str]:
    """Characters of `sequence` that ESM-2 does not accept as they are."""
    return set(sequence) - ESM_RESIDUES


def needs_cterm(sequence: str) -> bool:
    return len(sequence) > MAX_RESIDUES


def nterm_window(sequence: str) -> str:
    return sequence[:MAX_RESIDUES]


def cterm_window(sequence: str) -> str:
    return sequence[-MAX_RESIDUES:]


def window(sequence: str, name: str) -> str:
    if name == "nterm":
        return nterm_window(sequence)
    if name == "cterm":
        return cterm_window(sequence)
    raise ValueError(f"unknown window {name!r}; valid: {', '.join(WINDOWS)}")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare/test_seqwindow.py -q`
Expected: `9 passed, 1 skipped` (the package constants test needs torch).
Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_seqwindow.py -q`
Expected: `10 passed`.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/seqwindow.py tests/step1_compare/test_seqwindow.py
git commit -m "step1_compare: N- and C-terminal ESM-2 windows (Q10)"
```

### Task 2: Sequence universe and dedupe (05)

**Files:**
- Create: `analysis/step1_compare/sequence_sets.tsv`, `analysis/step1_compare/seqsets.py`, `analysis/step1_compare/05_prepare_sequences.py`
- Modify: `analysis/step1_compare/runinfo.py` (`PRODUCERS`)
- Test: `tests/step1_compare/test_phaseb_prepare.py`

**Interfaces:**
- Consumes: `seqwindow.bad_characters`, `seqwindow.needs_cterm` (Task 1); Phase A `seqhash.clean`, `seqhash.seq_sha256`, `sequences.read_fasta`, `sequences.fasta_key`, `truth_table.read_tsv`, `truth_table.write_tsv`, `manifest.sha256_file`, `runinfo.require_full`, `runinfo.says_all_sources`, `runinfo.atomic_write_all`, `paths.site_value`.
- Produces: `seqsets.SET_COLUMNS`, `seqsets.MEMBER_COLUMNS = ("set_id", "source_id", "gene_id", "seq_sha256", "length")`, `seqsets.UNIQUE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", "sequence")`; `seqsets.SequenceSetError(ValueError)`; `read_sets(path) -> list[dict]`; `resolve(row, work, downloads, site_value) -> Path`; `truth_members(rows, set_id="truth") -> (members, seqs)`; `fasta_members(set_id, path, kind) -> (members, seqs, empty)`; `unique_rows(seqs: dict[str, str]) -> list[dict]`. Script `05_prepare_sequences.py` with `OUTPUT_NAMES`, `phaseb_dir(work) -> Path`, `run(set_rows, work, downloads, site_value, provenance=None)`, `main(argv=None) -> int`. Outputs in `$STEP1_WORKDIR/phaseb/`: `sequence_members.tsv.gz`, `unique_sequences.tsv.gz`, `unique_sequences.fasta.gz`, `prepare_run.json`.

- [ ] **Step 1: Write the failing test**

```python
import gzip
import json

import pytest
import seqhash
import seqsets
import truth_table
from conftest import load_script

LONG = "MKV" + "ST" * 600  # 1,203 aa: has a C-terminal window
SHARED = "MKTAYIAKQRQISFVKSHFSRQ"


def _truth_rows():
    rows = []
    for source_id, gene_id, seq, label in (
        ("Scer_SGD", "S1", SHARED, "P-ext"),
        ("Scer_SGD", "S2", LONG, "N-sec"),
        ("Calb_CGD", "C1", SHARED, "N-int"),  # same sequence as S1 in another source
    ):
        rows.append(
            {
                "source_id": source_id,
                "gene_id": gene_id,
                "label": label,
                "fasta_id": gene_id,
                "length": str(len(seq)),
                "seq_sha256": seqhash.seq_sha256(seq),
                "sequence": seq,
            }
        )
    return rows


def _phaseb_inputs(tmp_path, truth_rows=None):
    work = tmp_path / "work"
    downloads = work / "downloads"
    downloads.mkdir(parents=True)
    truth_table.write_tsv(
        work / "truth_sequences.tsv.gz", list(_truth_rows()[0]), truth_rows or _truth_rows()
    )
    with gzip.open(work / "keyword_sequences.fasta.gz", "wt") as handle:
        handle.write(f">sp|P11111|KW1_YEAST desc\n{SHARED}\n>tr|Q22222|KW2_CANAL\nMSSPLLA*\n")
    (downloads / "prot.fasta").write_text(
        f">YAL001C desc\n{SHARED.lower()}\n>YAL002W\nMKLLV\n>YAL003W\n\n"
    )
    site = tmp_path / "site"
    site.mkdir()
    (site / "rs.fasta").write_text(">CIMG_1-t1-p1 | gene=CIMG_1\nMPPPP\n")
    for name in ("sequence_run.json", "keyword_tier_run.json"):
        (work / name).write_text(json.dumps({"all_sources": True}))
    sets = [
        {"set_id": "truth", "kind": "truth", "location": "truth_sequences.tsv.gz", "note": ""},
        {"set_id": "kw", "kind": "keyword", "location": "keyword_sequences.fasta.gz", "note": ""},
        {"set_id": "Scer_proteome", "kind": "download", "location": "prot.fasta", "note": ""},
        {"set_id": "Cimm", "kind": "site", "location": "cocci:rs.fasta", "note": ""},
    ]
    sets_path = tmp_path / "sets.tsv"
    truth_table.write_tsv(sets_path, seqsets.SET_COLUMNS, sets)
    return work, downloads, site, sets_path


def test_prepare_dedupes_by_hash_across_sets(tmp_path, monkeypatch):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 0
    out = work / "phaseb"
    members = truth_table.read_tsv(out / "sequence_members.tsv.gz")
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    # 3 truth + 2 keyword + 2 proteome (one empty record skipped) + 1 site = 8 members
    assert len(members) == 8
    # SHARED (5 members), LONG, MSSPLLA, MKLLV, MPPPP
    assert len(unique) == 5
    assert [r["seq_sha256"] for r in unique] == sorted(r["seq_sha256"] for r in unique)
    assert [r["row"] for r in unique] == ["0", "1", "2", "3", "4"]
    shared = seqhash.seq_sha256(SHARED)
    assert sum(m["seq_sha256"] == shared for m in members) == 4
    kw = {m["gene_id"] for m in members if m["set_id"] == "kw"}
    assert kw == {"P11111", "Q22222"}
    cimm = [m for m in members if m["set_id"] == "Cimm"]
    assert cimm[0]["gene_id"] == "CIMG_1-t1-p1" and cimm[0]["source_id"] == "Cimm"
    log = json.loads((out / "prepare_run.json").read_text())
    assert log["all_sources"] is True
    assert log["inputs"]["Scer_proteome"]["empty_records"] == 1
    assert log["over_max_residues"] == 1


def test_cterm_row_numbers_only_long_sequences(tmp_path):
    seqs = {seqhash.seq_sha256(s): s for s in (SHARED, LONG, "A" * 1023, "A" * 1022)}
    rows = seqsets.unique_rows(seqs)
    assert len(rows) == 4
    long_rows = [r for r in rows if int(r["length"]) > 1022]
    assert sorted(r["cterm_row"] for r in long_rows) == ["0", "1"]
    assert all(r["cterm_row"] == "" for r in rows if int(r["length"]) <= 1022)


def test_prepare_stops_on_a_non_esm_character(tmp_path, monkeypatch, capsys):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    (downloads / "prot.fasta").write_text(">YAL009W\nMK-LV\n")
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "not ESM-2 residues" in capsys.readouterr().err
    assert not (work / "phaseb").exists()


def test_prepare_stops_on_a_stale_truth_hash(tmp_path, monkeypatch, capsys):
    rows = _truth_rows()
    rows[0]["seq_sha256"] = "0" * 64
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path, rows)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "re-run 02_attach_sequences.py" in capsys.readouterr().err


def test_prepare_refuses_a_partial_truth_set(tmp_path, monkeypatch, capsys):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    (work / "sequence_run.json").write_text(json.dumps({"all_sources": False}))
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "02_attach_sequences.py" in capsys.readouterr().err
    assert prepare.main(argv + ["--allow-partial-truth-set"]) == 0


def test_prepare_stops_on_a_missing_input(tmp_path, monkeypatch, capsys):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    (downloads / "prot.fasta").unlink()
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "not found" in capsys.readouterr().err


def test_same_gene_with_two_sequences_stops(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">G1\nMKV\n>G1\nMKL\n")
    with pytest.raises(seqsets.SequenceSetError, match="two different sequences"):
        seqsets.fasta_members("x", path, "download")


def test_unique_fasta_headers_are_hashes(tmp_path, monkeypatch):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 0
    text = gzip.open(work / "phaseb" / "unique_sequences.fasta.gz", "rt").read()
    headers = [line[1:] for line in text.splitlines() if line.startswith(">")]
    unique = truth_table.read_tsv(work / "phaseb" / "unique_sequences.tsv.gz")
    assert headers == [r["seq_sha256"] for r in unique]
    first = (work / "phaseb" / "unique_sequences.fasta.gz").read_bytes()
    assert prepare.main(argv) == 0
    assert (work / "phaseb" / "unique_sequences.fasta.gz").read_bytes() == first  # byte-stable


def test_real_sets_file_is_valid():
    import paths

    rows = seqsets.read_sets(paths.STEP1_DIR / "sequence_sets.tsv")
    assert [r["set_id"] for r in rows][:2] == ["truth", "uniprot_kw"]
    assert sum(r["kind"] == "download" for r in rows) == 8
    assert any(r["set_id"] == "Cimm_RS_proteome" for r in rows)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_prepare.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'seqsets'`.

- [ ] **Step 3: Write the sets file**

Write it with this command, so that the TAB characters are exact:

```bash
/usr/bin/python3.12 - <<'EOF'
from pathlib import Path

rows = [
    ('set_id', 'kind', 'location', 'note'),
    ('truth', 'truth', 'truth_sequences.tsv.gz', '02_attach_sequences.py output; one row per matched truth gene'),
    ('uniprot_kw', 'keyword', 'keyword_sequences.fasta.gz', '04_build_keyword_tier.py cache; surface.tsv accessions and literature seeds'),
    ('Scer_proteome', 'download', 'orf_trans_all.fasta.gz', 'whole proteome for the agreement matrix (spec 6)'),
    ('Calb_proteome', 'download', 'C_albicans_SC5314_A22_current_default_protein.fasta.gz', 'whole proteome for the agreement matrix (spec 6)'),
    ('Spom_proteome', 'download', 'peptide.fa.gz', 'whole proteome for the agreement matrix (spec 6)'),
    ('Afum_proteome', 'download', 'UP000002530.fasta.gz', 'whole proteome of a truth species'),
    ('Anid_proteome', 'download', 'UP000000560.fasta.gz', 'whole proteome of a truth species'),
    ('Cneo_H99_proteome', 'download', 'UP000010091.fasta.gz', 'whole proteome of a truth species'),
    ('Cneo_JEC21_proteome', 'download', 'UP000002149.fasta.gz', 'whole proteome of a truth species'),
    ('Umay_proteome', 'download', 'UP000000561.fasta.gz', 'whole proteome of a truth species'),
    ('Cimm_RS_proteome', 'site', 'cocci_pangenome:input_run2/CimmitisRS_FungiDB.fasta', 'C. immitis RS (spec 5, 6, D2); FungiDB file used by analysis/cocci_repeats'),
]
path = Path("analysis/step1_compare/sequence_sets.tsv")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text("".join("\t".join(r) + "\n" for r in rows))
EOF
```

- [ ] **Step 4: Write `seqsets.py`**

```python
"""Phase B sequence universe: members of every input set and the unique sequences.

Standard library only. A member is one protein of one input set (one row per set_id,
source_id, gene_id). Members with the same cleaned sequence share one seq_sha256 and are
embedded and scored once. Unique sequences are sorted by seq_sha256; `row` is the position in
that order and is the row of the N-terminal embedding matrix. `cterm_row` numbers the
sequences longer than MAX_RESIDUES in the same order and is the row of the C-terminal matrix.
"""

from collections.abc import Iterable
from pathlib import Path

import seqhash
import sequences
import seqwindow
import truth_table

SET_COLUMNS = ("set_id", "kind", "location", "note")
KINDS = ("truth", "keyword", "download", "site")
MEMBER_COLUMNS = ("set_id", "source_id", "gene_id", "seq_sha256", "length")
UNIQUE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", "sequence")


class SequenceSetError(ValueError):
    """An input set is malformed, missing, or holds a sequence ESM-2 cannot take as it is."""


def read_sets(path: str | Path) -> list[dict[str, str]]:
    rows = truth_table.read_tsv(path)
    seen = set()
    for row in rows:
        missing = [c for c in SET_COLUMNS if c not in row]
        if missing:
            raise SequenceSetError(f"{path}: row {row} lacks {missing}")
        if row["kind"] not in KINDS:
            raise SequenceSetError(f"{path}: set {row['set_id']} has unknown kind {row['kind']}")
        if row["set_id"] in seen:
            raise SequenceSetError(f"{path}: set_id {row['set_id']} occurs twice")
        seen.add(row["set_id"])
    return rows


def resolve(row: dict[str, str], work: Path, downloads: Path, site_value) -> Path:
    """Path of a set's input. `site` locations are `<site.yaml key>:<relative path>`."""
    kind, location = row["kind"], row["location"]
    if kind in ("truth", "keyword"):
        return Path(work) / location
    if kind == "download":
        return Path(downloads) / location
    key, _, rel = location.partition(":")
    if not rel:
        raise SequenceSetError(f"set {row['set_id']}: site location needs 'key:path'")
    return Path(site_value(key)) / rel


def _check(seq: str, where: str) -> None:
    bad = seqwindow.bad_characters(seq)
    if bad:
        raise SequenceSetError(f"{where}: characters {sorted(bad)} are not ESM-2 residues")


def truth_members(rows: Iterable[dict[str, str]], set_id: str = "truth"):
    """Members from truth_sequences.tsv.gz rows. The stored seq_sha256 is checked."""
    members, seqs = [], {}
    for r in rows:
        seq = r["sequence"]
        digest = seqhash.seq_sha256(seq)
        if digest != r["seq_sha256"] or seq != seqhash.clean(seq):
            raise SequenceSetError(
                f"truth_sequences: {r['source_id']}:{r['gene_id']} has a stored seq_sha256 or "
                "sequence that does not match the cleaned sequence; re-run 02_attach_sequences.py"
            )
        _check(seq, f"{set_id} {r['source_id']}:{r['gene_id']}")
        members.append(_member(set_id, r["source_id"], r["gene_id"], digest, seq))
        seqs[digest] = seq
    return members, seqs


def fasta_members(set_id: str, path: str | Path, kind: str):
    """Members from a FASTA file. Return (members, sequences by hash, empty records).

    gene_id is the UniProt accession for kind `keyword`, else the first header token.
    A gene_id that occurs twice with the same cleaned sequence is kept once; with a different
    sequence it raises SequenceSetError."""
    members, seqs, by_gene = [], {}, {}
    empty = 0
    for header, raw in sequences.read_fasta(path):
        seq = seqhash.clean(raw)
        if not seq:
            empty += 1
            continue
        if kind == "keyword":
            gene_id = sequences.fasta_key(header, "uniprot")
            if gene_id is None:
                raise SequenceSetError(f"{set_id}: header {header[:40]!r} is not UniProt FASTA")
        else:
            gene_id = header.split()[0]
        digest = seqhash.seq_sha256(seq)
        if gene_id in by_gene:
            if by_gene[gene_id] != digest:
                raise SequenceSetError(f"{set_id}: {gene_id} has two different sequences")
            continue
        _check(seq, f"{set_id} {gene_id}")
        by_gene[gene_id] = digest
        members.append(_member(set_id, set_id, gene_id, digest, seq))
        seqs[digest] = seq
    return members, seqs, empty


def _member(set_id, source_id, gene_id, digest, seq) -> dict[str, str]:
    return {
        "set_id": set_id,
        "source_id": source_id,
        "gene_id": gene_id,
        "seq_sha256": digest,
        "length": str(len(seq)),
    }


def unique_rows(seqs: dict[str, str]) -> list[dict[str, str]]:
    """One row per unique sequence, sorted by seq_sha256, with row and cterm_row."""
    rows, cterm = [], 0
    for i, digest in enumerate(sorted(seqs)):
        seq = seqs[digest]
        crow = ""
        if seqwindow.needs_cterm(seq):
            crow = str(cterm)
            cterm += 1
        rows.append(
            {
                "row": str(i),
                "seq_sha256": digest,
                "length": str(len(seq)),
                "cterm_row": crow,
                "sequence": seq,
            }
        )
    return rows
```

- [ ] **Step 5: Write `05_prepare_sequences.py`**

```python
#!/usr/bin/env python3
"""Phase B step 1: collect every input set and dedupe the sequences by seq_sha256.

Reads the sets in sequence_sets.tsv: truth_sequences.tsv.gz (02), keyword_sequences.fasta.gz
(04), the proteome FASTA files in $STEP1_WORKDIR/downloads, and the C. immitis RS FASTA (path
from config/site.yaml). Writes to $STEP1_WORKDIR/phaseb/:

  sequence_members.tsv.gz   one row per (set_id, source_id, gene_id)
  unique_sequences.tsv.gz   one row per unique sequence (row, seq_sha256, length, cterm_row,
                            sequence), sorted by seq_sha256
  unique_sequences.fasta.gz the same sequences as FASTA; the header is the seq_sha256
  prepare_run.json          input hashes, counts, all_sources, git commit, arguments

STOP (exit 2, no output written): a missing input; a run log (sequence_run.json,
keyword_tier_run.json) that does not say all_sources: true, unless --allow-partial-truth-set;
a stored seq_sha256 that does not match its sequence; a character that ESM-2 does not accept;
one gene_id with two different sequences in one set. Columns are listed in COLUMNS.md.
"""

import argparse
import gzip
import sys
from pathlib import Path

import manifest
import paths
import runinfo
import seqsets
import truth_table

OUTPUT_NAMES = (
    "sequence_members.tsv.gz",
    "unique_sequences.tsv.gz",
    "unique_sequences.fasta.gz",
    "prepare_run.json",
)


def phaseb_dir(work: Path) -> Path:
    return Path(work) / "phaseb"


def write_fasta_gz(path: Path, rows) -> None:
    """FASTA with the seq_sha256 as header, gzip with a fixed mtime (byte-stable)."""
    data = "".join(f">{r['seq_sha256']}\n{r['sequence']}\n" for r in rows).encode("ascii")
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(data)


def run(set_rows, work: Path, downloads: Path, site_value, provenance: dict | None = None):
    members, seqs, inputs = [], {}, {}
    for row in set_rows:
        path = seqsets.resolve(row, work, downloads, site_value)
        if not path.exists():
            raise seqsets.SequenceSetError(f"set {row['set_id']}: {path} not found")
        if row["kind"] == "truth":
            got, got_seqs = seqsets.truth_members(truth_table.read_tsv(path), row["set_id"])
            empty = 0
        else:
            got, got_seqs, empty = seqsets.fasta_members(row["set_id"], path, row["kind"])
        if not got:
            raise seqsets.SequenceSetError(f"set {row['set_id']}: {path} gave no sequence")
        members += got
        seqs.update(got_seqs)
        inputs[row["set_id"]] = {
            "file": str(path),
            "sha256": manifest.sha256_file(path),
            "members": len(got),
            "empty_records": empty,
        }
    unique = seqsets.unique_rows(seqs)
    lengths = [int(r["length"]) for r in unique]
    log = {
        "inputs": inputs,
        "members": len(members),
        "unique_sequences": len(unique),
        "unique_residues": sum(lengths),
        "over_max_residues": sum(r["cterm_row"] != "" for r in unique),
        "all_sources": None,
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
        "arguments": [],
    }
    log.update(provenance or {})
    out = phaseb_dir(work)
    runinfo.atomic_write_all(
        out,
        {
            "sequence_members.tsv.gz": lambda p: truth_table.write_tsv(
                p, seqsets.MEMBER_COLUMNS, members
            ),
            "unique_sequences.tsv.gz": lambda p: truth_table.write_tsv(
                p, seqsets.UNIQUE_COLUMNS, unique
            ),
            "unique_sequences.fasta.gz": lambda p: write_fasta_gz(p, unique),
            "prepare_run.json": lambda p: runinfo.write_json(p, log),
        },
    )
    return members, unique, log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sets", default=str(paths.STEP1_DIR / "sequence_sets.tsv"))
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument(
        "--allow-partial-truth-set",
        action="store_true",
        help="use inputs whose run logs do not say all_sources: true",
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    downloads = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    try:
        set_rows = seqsets.read_sets(args.sets)
        kinds = {r["kind"] for r in set_rows}
        allow = args.allow_partial_truth_set
        logs = []
        if "truth" in kinds:
            logs.append("sequence_run.json")
        if "keyword" in kinds:
            logs.append("keyword_tier_run.json")
        for name in logs:
            runinfo.require_full(work / name, name, allow)
        provenance = {
            "all_sources": all(runinfo.says_all_sources(work / n) for n in logs),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        _, unique, log = run(set_rows, work, downloads, paths.site_value, provenance)
    except (seqsets.SequenceSetError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(
        f"members={log['members']} unique={log['unique_sequences']} "
        f"residues={log['unique_residues']} over_1022={log['over_max_residues']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 6: Name the producer of `keyword_tier_run.json` in `runinfo.PRODUCERS`**

Add only the `keyword_tier_run.json` line now (the `d8_run.json` line is added in Task 8). The final state after Task 8 is:

```diff
--- a/analysis/step1_compare/runinfo.py
+++ b/analysis/step1_compare/runinfo.py
@@ -16,6 +16,8 @@
 PRODUCERS = {
     "extract_log.json": "01_extract_go_truth.py",
     "sequence_run.json": "02_attach_sequences.py",
+    "keyword_tier_run.json": "04_build_keyword_tier.py",
+    "d8_run.json": "03_triage_pm.py",
 }
 
 
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_prepare.py tests/step1_compare/test_paths.py -q`
Expected: all pass (`test_phaseb_prepare.py`: 9 passed). `test_every_module_imports_only_stdlib_or_local` still passes.

- [ ] **Step 8: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/sequence_sets.tsv analysis/step1_compare/seqsets.py \
  analysis/step1_compare/05_prepare_sequences.py analysis/step1_compare/runinfo.py \
  tests/step1_compare/test_phaseb_prepare.py
git commit -m "step1_compare: 05 collects the Phase B input sets and dedupes by seq_sha256"
```

### Task 3: Chunk plan and job sizing (06)

**Files:**
- Create: `analysis/step1_compare/chunk_plan.py`, `analysis/step1_compare/06_plan_embedding.py`
- Test: `tests/step1_compare/test_chunk_plan.py`

**Interfaces:**
- Consumes: `seqwindow.window`, `seqwindow.needs_cterm` (Task 1); `seqsets.unique_rows`, `seqsets.UNIQUE_COLUMNS` (Task 2, tests); `truth_table`, `manifest`, `runinfo`, `paths`.
- Produces: `chunk_plan.CHUNK_SECONDS = 300`, `MIN_CHUNK = 10_000`, `SAFETY = 1.25`, `TARGET_SECONDS = 4500`, `PLAN_COLUMNS`, `MEMBER_COLUMNS`; dataclass `Chunk(chunk_id: str, window: str, rows: list[int], hashes: list[str], residues: int)` with property `members_sha256`; `members_digest(hashes) -> str`; `window_entries(unique_rows, window) -> list[tuple[int, str, int]]`; `build_chunks(unique_rows, window, chunk_residues) -> list[Chunk]`; `chunk_residues_from_rate(min_rate: float) -> int`; `plan_jobs(rates: dict[str, float], load_s: dict[str, float], residues: int) -> dict` (keys `total_seconds`, `n_jobs`, `seconds_per_job`, `time_minutes`); `chunks_for_job(chunk_ids, job_index, job_count) -> list[str]`; `read_plan(out_dir) -> list[Chunk]`. Script `06_plan_embedding.py` with `MODELS`, `OUTPUT_NAMES`, `PlanError`, `rates_from_j0(j0, models)`, `run(...)`, `main(argv=None)`. Outputs: `chunk_plan.tsv`, `chunk_members.tsv.gz`, `job_plan.json` (keys include `batch_size[model]`, `n_jobs`, `time_minutes`, `rate_source`).

- [ ] **Step 1: Write the failing test**

```python
import json

import chunk_plan
import pytest
import seqhash
import seqsets
import truth_table
from conftest import load_script


def _unique(lengths):
    seqs = {}
    for i, n in enumerate(lengths):
        seq = ("MKV" + "ACDEFGHIKLMNPQRSTVWY"[i % 20] * n)[:n]
        seqs[seqhash.seq_sha256(seq)] = seq
    return seqsets.unique_rows(seqs)


def test_chunks_cover_every_sequence_once_and_respect_the_budget():
    unique = _unique([50, 60, 70, 300, 900, 1100, 1500, 40, 41, 2000])
    chunks = chunk_plan.build_chunks(unique, "nterm", 1200)
    rows = [r for c in chunks for r in c.rows]
    assert sorted(rows) == list(range(len(unique)))
    for c in chunks:
        assert c.residues <= 1200 or len(c.rows) == 1
    ids = [c.chunk_id for c in chunks]
    assert ids == [f"nterm_{i:04d}" for i in range(len(chunks))]


def test_cterm_chunks_hold_only_long_sequences():
    unique = _unique([50, 1100, 1500, 2000, 1022])
    chunks = chunk_plan.build_chunks(unique, "cterm", 5000)
    rows = sorted(r for c in chunks for r in c.rows)
    long_rows = sorted(int(u["row"]) for u in unique if int(u["length"]) > 1022)
    assert rows == long_rows
    assert all(c.residues == 1022 * len(c.rows) for c in chunks)


def test_chunk_membership_is_deterministic_and_length_sorted():
    unique = _unique([500, 20, 300, 800, 45, 1000, 60])
    a = chunk_plan.build_chunks(unique, "nterm", 900)
    b = chunk_plan.build_chunks(list(reversed(unique)), "nterm", 900)
    assert [c.members_sha256 for c in a] == [c.members_sha256 for c in b]
    lengths = [int(unique[r]["length"]) for c in a for r in c.rows]
    assert lengths == sorted(lengths)


def test_chunk_residues_from_rate_formula():
    assert chunk_plan.chunk_residues_from_rate(29_500) == 8_850_000
    assert chunk_plan.chunk_residues_from_rate(1.0) == chunk_plan.MIN_CHUNK
    with pytest.raises(ValueError):
        chunk_plan.chunk_residues_from_rate(0)


def test_plan_jobs_worked_example():
    # Assumed rate (review 9.1, ESM-2 150M) until J0 measures 8M and 35M.
    plan = chunk_plan.plan_jobs(
        {"esm2_t6_8M_UR50D": 29_500.0, "esm2_t12_35M_UR50D": 29_500.0},
        {"esm2_t6_8M_UR50D": 10.0, "esm2_t12_35M_UR50D": 10.0},
        37_860_229,
    )
    assert plan["total_seconds"] == pytest.approx(2586.8, abs=0.1)
    assert plan["n_jobs"] == 1
    assert plan["time_minutes"] == 75


def test_plan_jobs_splits_slow_rates_into_jobs_under_the_target():
    plan = chunk_plan.plan_jobs({"m": 2_000.0}, {"m": 0.0}, 37_860_229)
    assert plan["n_jobs"] == 6  # 18,930 s * 1.25 / 4,500 s = 5.26 -> 6
    assert plan["seconds_per_job"] * chunk_plan.SAFETY <= chunk_plan.TARGET_SECONDS


def test_chunks_for_job_round_robin():
    ids = [f"c{i}" for i in range(7)]
    jobs = [chunk_plan.chunks_for_job(ids, i, 3) for i in range(3)]
    assert jobs == [["c0", "c3", "c6"], ["c1", "c4"], ["c2", "c5"]]
    assert sorted(sum(jobs, [])) == sorted(ids)
    with pytest.raises(ValueError):
        chunk_plan.chunks_for_job(ids, 3, 3)


def _plan_inputs(tmp_path, lengths):
    out = tmp_path / "phaseb"
    out.mkdir(parents=True)
    unique = _unique(lengths)
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    return out, unique


def test_plan_script_reads_j0_and_writes_a_consistent_plan(tmp_path):
    out, unique = _plan_inputs(tmp_path, [50, 400, 1100, 1300, 900])
    j0 = {
        "best": {
            "esm2_t6_8M_UR50D": {"batch_size": 32, "residues_per_s": 5.0},
            "esm2_t12_35M_UR50D": {"batch_size": 16, "residues_per_s": 4.0},
        },
        "model_load_s": {"esm2_t6_8M_UR50D": 1.0, "esm2_t12_35M_UR50D": 2.0},
    }
    (out / "j0").mkdir()
    (out / "j0" / "throughput.json").write_text(json.dumps(j0))
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--chunk-residues", "1000"]) == 0
    plan = json.loads((out / "job_plan.json").read_text())
    assert plan["rate_source"] == "J0" and plan["batch_size"]["esm2_t12_35M_UR50D"] == 16
    assert plan["residues_per_model"] == 50 + 400 + 1022 + 1022 + 900 + 1022 + 1022
    chunks = chunk_plan.read_plan(out)
    assert sum(len(c.rows) for c in chunks if c.window == "nterm") == len(unique)
    assert sum(len(c.rows) for c in chunks if c.window == "cterm") == 2


def test_plan_script_stops_without_j0(tmp_path, capsys):
    _plan_inputs(tmp_path, [50, 60])
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path)]) == 2
    assert "run J0 first" in capsys.readouterr().err
    assert not (tmp_path / "phaseb" / "chunk_plan.tsv").exists()
    assert plan_script.main(["--work-dir", str(tmp_path), "--rate", "100"]) == 0
    plan = json.loads((tmp_path / "phaseb" / "job_plan.json").read_text())
    assert plan["rate_source"] == "assumed"


def test_read_plan_detects_edited_members(tmp_path):
    out, _ = _plan_inputs(tmp_path, [50, 60, 70])
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--rate", "100"]) == 0
    rows = truth_table.read_tsv(out / "chunk_members.tsv.gz")
    rows[0]["seq_sha256"] = "f" * 64
    truth_table.write_tsv(out / "chunk_members.tsv.gz", chunk_plan.MEMBER_COLUMNS, rows)
    with pytest.raises(ValueError, match="members differ"):
        chunk_plan.read_plan(out)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_chunk_plan.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'chunk_plan'`.

- [ ] **Step 3: Write `chunk_plan.py`**

```python
"""Embedding chunks and job sizing for Phase B (J2). Standard library only.

A chunk is the unit of resumability: a fixed, ordered list of unique sequences for one window.
Entries are sorted by (window length, seq_sha256), so a chunk holds sequences of similar length
(little padding) and its membership depends only on the unique set and `chunk_residues`.

A job is a group of chunks. Job count and chunk size come from the measured J0 throughput:

  chunk_residues = max(MIN_CHUNK, floor(min_m(r_m) * CHUNK_SECONDS / 10,000) * 10,000)
  T_total        = sum_m (R / r_m + load_m)          R = residues of all windows to embed
  n_jobs         = max(1, ceil(T_total * SAFETY / TARGET_SECONDS))
  time_minutes   = ceil((T_total / n_jobs * 1.5 + 600) / 60)

r_m is residues per second of model m at its best batch size, and load_m its load time, both
from the J0 throughput JSON. Chunk k of the plan goes to job k mod n_jobs.
"""

import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path

import seqwindow
import truth_table

CHUNK_SECONDS = 300  # at most about 5 min of work is lost when a job is killed
MIN_CHUNK = 10_000
SAFETY = 1.25
TARGET_SECONDS = 4500  # 1.25 h: the middle of the 1 to 1.5 h job size rule
PLAN_COLUMNS = ("chunk_id", "window", "n_seqs", "residues", "members_sha256")
MEMBER_COLUMNS = ("chunk_id", "position", "row", "seq_sha256")


@dataclass(slots=True)
class Chunk:
    chunk_id: str
    window: str
    rows: list[int] = field(default_factory=list)
    hashes: list[str] = field(default_factory=list)
    residues: int = 0

    @property
    def members_sha256(self) -> str:
        return members_digest(self.hashes)


def members_digest(hashes) -> str:
    return hashlib.sha256("\n".join(hashes).encode("ascii")).hexdigest()


def window_entries(unique_rows, window: str) -> list[tuple[int, str, int]]:
    """(window length, seq_sha256, row) for every sequence that has this window."""
    out = []
    for r in unique_rows:
        seq = r["sequence"]
        if window == "cterm" and not seqwindow.needs_cterm(seq):
            continue
        out.append((len(seqwindow.window(seq, window)), r["seq_sha256"], int(r["row"])))
    return sorted(out)


def build_chunks(unique_rows, window: str, chunk_residues: int) -> list[Chunk]:
    if chunk_residues < 1:
        raise ValueError("chunk_residues must be positive")
    chunks: list[Chunk] = []
    current = None
    for length, digest, row in window_entries(unique_rows, window):
        if current is None or (current.rows and current.residues + length > chunk_residues):
            current = Chunk(f"{window}_{len(chunks):04d}", window)
            chunks.append(current)
        current.rows.append(row)
        current.hashes.append(digest)
        current.residues += length
    return chunks


def chunk_residues_from_rate(min_rate: float) -> int:
    if min_rate <= 0:
        raise ValueError("throughput must be positive")
    return max(MIN_CHUNK, int(min_rate * CHUNK_SECONDS // 10_000) * 10_000)


def plan_jobs(rates: dict[str, float], load_s: dict[str, float], residues: int) -> dict:
    """Job count and wall time for embedding `residues` residues with every model in `rates`."""
    if not rates:
        raise ValueError("no model rates")
    total = sum(residues / rates[m] + load_s.get(m, 0.0) for m in rates)
    n_jobs = max(1, math.ceil(total * SAFETY / TARGET_SECONDS))
    per_job = total / n_jobs
    return {
        "total_seconds": round(total, 1),
        "n_jobs": n_jobs,
        "seconds_per_job": round(per_job, 1),
        "time_minutes": math.ceil((per_job * 1.5 + 600) / 60),
    }


def chunks_for_job(chunk_ids: list[str], job_index: int, job_count: int) -> list[str]:
    if not 0 <= job_index < job_count:
        raise ValueError(f"job index {job_index} outside 0..{job_count - 1}")
    return [c for k, c in enumerate(chunk_ids) if k % job_count == job_index]


def read_plan(out_dir) -> list[Chunk]:
    """Chunks from chunk_plan.tsv and chunk_members.tsv.gz, checked against each other."""
    out_dir = Path(out_dir)
    plan = truth_table.read_tsv(out_dir / "chunk_plan.tsv")
    chunks = {r["chunk_id"]: Chunk(r["chunk_id"], r["window"]) for r in plan}
    for m in truth_table.read_tsv(out_dir / "chunk_members.tsv.gz"):
        chunk = chunks.get(m["chunk_id"])
        if chunk is None or int(m["position"]) != len(chunk.rows):
            raise ValueError(f"chunk_members.tsv.gz: unexpected row {m}")
        chunk.rows.append(int(m["row"]))
        chunk.hashes.append(m["seq_sha256"])
    for r in plan:
        chunk = chunks[r["chunk_id"]]
        if len(chunk.rows) != int(r["n_seqs"]) or chunk.members_sha256 != r["members_sha256"]:
            raise ValueError(f"chunk {r['chunk_id']}: members differ from chunk_plan.tsv")
        chunk.residues = int(r["residues"])
    return [chunks[r["chunk_id"]] for r in plan]
```

- [ ] **Step 4: Write `06_plan_embedding.py`**

```python
#!/usr/bin/env python3
"""Phase B step 2: cut the unique sequences into embedding chunks and size the J2 jobs.

Reads $STEP1_WORKDIR/phaseb/unique_sequences.tsv.gz and the J0 throughput JSON
($STEP1_WORKDIR/phaseb/j0/throughput.json). Writes to $STEP1_WORKDIR/phaseb/:

  chunk_plan.tsv          one row per chunk (chunk_id, window, n_seqs, residues, members_sha256)
  chunk_members.tsv.gz    one row per chunk member (chunk_id, position, row, seq_sha256)
  job_plan.json           the formula inputs (rates, load times, residues), chunk_residues,
                          n_jobs, time_minutes, batch size per model

The formula is in chunk_plan.py. --rate gives one assumed residues/s rate for every model
instead of J0; the plan then records `rate_source: assumed`. Use it only for a dry plan.
STOP (exit 2, no output): no throughput JSON and no --rate; a model without a successful J0 run.
"""

import argparse
import json
import sys
from pathlib import Path

import chunk_plan
import manifest
import paths
import runinfo
import truth_table

OUTPUT_NAMES = ("chunk_plan.tsv", "chunk_members.tsv.gz", "job_plan.json")
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


class PlanError(ValueError):
    """The throughput input is missing or incomplete."""


def rates_from_j0(j0: dict, models) -> tuple[dict, dict, dict]:
    rates, load_s, batch = {}, {}, {}
    for m in models:
        best = j0.get("best", {}).get(m)
        if not best or best.get("residues_per_s", 0) <= 0:
            raise PlanError(f"J0 throughput has no successful run for {m}")
        rates[m] = float(best["residues_per_s"])
        batch[m] = int(best["batch_size"])
        load_s[m] = float(j0.get("model_load_s", {}).get(m, 0.0))
    return rates, load_s, batch


def run(unique, rates, load_s, batch, out_dir: Path, chunk_residues=None, provenance=None):
    size = chunk_residues or chunk_plan.chunk_residues_from_rate(min(rates.values()))
    chunks = []
    for window in ("nterm", "cterm"):
        chunks += chunk_plan.build_chunks(unique, window, size)
    residues = sum(c.residues for c in chunks)
    plan = {
        "models": list(rates),
        "rates_residues_per_s": rates,
        "model_load_s": load_s,
        "batch_size": batch,
        "chunk_residues": size,
        "chunks": len(chunks),
        "residues_per_model": residues,
        "unique_sequences": len(unique),
        **chunk_plan.plan_jobs(rates, load_s, residues),
    }
    plan.update(provenance or {})
    plan_rows = [
        {
            "chunk_id": c.chunk_id,
            "window": c.window,
            "n_seqs": str(len(c.rows)),
            "residues": str(c.residues),
            "members_sha256": c.members_sha256,
        }
        for c in chunks
    ]
    member_rows = [
        {"chunk_id": c.chunk_id, "position": str(i), "row": str(r), "seq_sha256": h}
        for c in chunks
        for i, (r, h) in enumerate(zip(c.rows, c.hashes, strict=True))
    ]
    runinfo.atomic_write_all(
        out_dir,
        {
            "chunk_plan.tsv": lambda p: truth_table.write_tsv(
                p, chunk_plan.PLAN_COLUMNS, plan_rows
            ),
            "chunk_members.tsv.gz": lambda p: truth_table.write_tsv(
                p, chunk_plan.MEMBER_COLUMNS, member_rows
            ),
            "job_plan.json": lambda p: runinfo.write_json(p, plan),
        },
    )
    return chunks, plan


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--throughput", default=None, help="default: phaseb/j0/throughput.json")
    parser.add_argument("--rate", type=float, default=None, help="assumed residues/s (dry plan)")
    parser.add_argument("--chunk-residues", type=int, default=None)
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    out = work / "phaseb"
    j0_path = Path(args.throughput) if args.throughput else out / "j0" / "throughput.json"
    try:
        unique_path = out / "unique_sequences.tsv.gz"
        unique = truth_table.read_tsv(unique_path)
        if args.rate is not None:
            rates = {m: args.rate for m in args.models}
            load_s, batch = {m: 0.0 for m in args.models}, {m: 16 for m in args.models}
            source = {"rate_source": "assumed", "throughput_sha256": ""}
        else:
            if not j0_path.exists():
                raise PlanError(f"{j0_path} not found; run J0 first or give --rate for a dry plan")
            rates, load_s, batch = rates_from_j0(json.loads(j0_path.read_text()), args.models)
            source = {"rate_source": "J0", "throughput_sha256": manifest.sha256_file(j0_path)}
        provenance = {
            **source,
            "unique_sequences_sha256": manifest.sha256_file(unique_path),
            "git_commit": runinfo.git_commit(),
            "python": runinfo.python_version(),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        chunks, plan = run(unique, rates, load_s, batch, out, args.chunk_residues, provenance)
    except (PlanError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(
        f"chunks={plan['chunks']} chunk_residues={plan['chunk_residues']} "
        f"n_jobs={plan['n_jobs']} time_minutes={plan['time_minutes']} "
        f"rate_source={plan['rate_source']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare/test_chunk_plan.py tests/step1_compare/test_paths.py -q`
Expected: all pass (`test_chunk_plan.py`: 10 passed).

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/chunk_plan.py analysis/step1_compare/06_plan_embedding.py \
  tests/step1_compare/test_chunk_plan.py
git commit -m "step1_compare: 06 cuts embedding chunks and sizes J2 from the J0 throughput"
```

### Task 4: Chunk store and J2 runner

**Files:**
- Create: `analysis/step1_compare/jobs/embed_store.py`, `analysis/step1_compare/jobs/embed_chunks.py`, `tests/step1_compare/phaseb_fixture.py`
- Modify: `tests/step1_compare/conftest.py` (one line)
- Test: `tests/step1_compare/test_phaseb_embed.py`

**Interfaces:**
- Consumes: `chunk_plan.read_plan`, `chunk_plan.chunks_for_job`, `Chunk.members_sha256` (Task 3); `seqwindow.window` (Task 1); `06_plan_embedding.main` (tests, Task 3); `surface_glyco.embeddings.get_esm_embeddings(sequences, model_name, batch_size, device, repr_layer, return_indices=True) -> (array, ids, kept)`.
- Produces: `embed_store.DTYPE = np.float32`; `array_sha256(arr) -> str`; `chunk_paths(publish, model, chunk_id) -> (npy, json)`; `save_chunk(publish, scratch, model, chunk_id, arr, meta) -> dict` (refuses NaN/inf and non-2-D arrays; writes to scratch, then copies npy and then json atomically to publish); `load_chunk(publish, model, chunk_id) -> (arr, meta)` (ValueError on shape, dtype or hash mismatch); `chunk_is_done(publish, model, chunk_id, members_sha256, n) -> bool`. `embed_chunks.MODELS`, `REPR_LAYER = 6`, `EmbedError(RuntimeError)`, `embed_window_sequences(seqs, model, batch_size, device, layer) -> np.ndarray`, `run_job(work, models, job_index, job_count, device, scratch, batch_size=None, stop_after=None) -> dict` (keys `done`, `skipped`, `stopped`), `main(argv=None)` (exit 0, 2 STOP, 3 test hook). Chunk files: `$STEP1_WORKDIR/phaseb/emb/<model>/<chunk_id>.npy` and `.json`. Test helper `phaseb_fixture.make_plan(work, chunk_residues=400) -> (out, unique)`, `phaseb_fixture.run_cpu(work, scratch, **kw) -> dict`, `MODEL`, `LONG`, `SEQS`.

- [ ] **Step 1: Put `jobs/` on the test path**

```diff
--- a/tests/step1_compare/conftest.py
+++ b/tests/step1_compare/conftest.py
@@ -9,6 +9,7 @@
 STEP1_DIR = Path(__file__).resolve().parents[2] / "analysis" / "step1_compare"
 TESTS_DIR = Path(__file__).resolve().parent
 sys.path.insert(0, str(STEP1_DIR))
+sys.path.insert(0, str(STEP1_DIR / "jobs"))
 sys.path.insert(0, str(TESTS_DIR))
 
 
```

- [ ] **Step 2: Write the test helper**

```python
"""Tiny Phase B embedding plan for the J2 runner and assembly tests (CPU, ESM-2 8M)."""

import seqhash
import seqsets
import truth_table
from conftest import load_script

MODEL = "esm2_t6_8M_UR50D"
LONG = "MKLSTA" + "STPSSTSA" * 140  # 1,126 aa: the only sequence with a C-terminal window
SEQS = [
    "MKTLLVAGLLSSAAFA",
    "MSTTSSTTSTPSSTSA" * 4,
    "MKVLAAGIVALLLAAGCSSS" * 3,
    "MQRSLLLAVAALATPAFAAS" * 6,
    "MAEEKKAVEEVKSAGEW" * 2,
    LONG,
]


def make_plan(work, chunk_residues=400):
    """Write unique_sequences.tsv.gz and a dry chunk plan (assumed rate) under work/phaseb."""
    out = work / "phaseb"
    out.mkdir(parents=True)
    unique = seqsets.unique_rows({seqhash.seq_sha256(s): s for s in SEQS})
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    plan = load_script("06_plan_embedding")
    argv = ["--work-dir", str(work), "--rate", "100", "--models", MODEL]
    assert plan.main(argv + ["--chunk-residues", str(chunk_residues)]) == 0
    return out, unique


def run_cpu(work, scratch, **kw):
    import embed_chunks
    import torch

    cpu = torch.device("cpu")
    return embed_chunks.run_job(work, [MODEL], 0, 1, cpu, scratch, batch_size=2, **kw)
```

- [ ] **Step 3: Write the failing test**

```python
"""J2 runner on CPU with ESM-2 8M (needs torch, fair-esm and the 8M weights)."""

import os
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")
torch = pytest.importorskip("torch")
pytest.importorskip("esm")

import chunk_plan  # noqa: E402
import embed_chunks  # noqa: E402
import embed_store  # noqa: E402
from conftest import STEP1_DIR  # noqa: E402
from phaseb_fixture import LONG, MODEL, make_plan, run_cpu  # noqa: E402

CPU = torch.device("cpu")


def test_job_embeds_every_chunk_member_once(tmp_path):
    work = tmp_path / "w"
    out, unique = make_plan(work)
    result = run_cpu(work, tmp_path / "scratch")
    chunks = chunk_plan.read_plan(out)
    assert result == {"done": len(chunks), "skipped": 0, "stopped": False}
    for c in chunks:
        arr, meta = embed_store.load_chunk(out / "emb", MODEL, c.chunk_id)
        assert arr.shape == (len(c.rows), 320) and arr.dtype == np.float32
        assert np.isfinite(arr).all()
        assert meta["members_sha256"] == c.members_sha256 and meta["repr_layer"] == 6


def test_rerun_skips_finished_chunks_after_hash_check(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    first = run_cpu(work, tmp_path / "scratch")
    again = run_cpu(work, tmp_path / "scratch")
    assert again == {"done": 0, "skipped": first["done"], "stopped": False}


def test_corrupted_chunk_is_recomputed_with_the_same_hash(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk_id = chunk_plan.read_plan(out)[0].chunk_id
    npy, js = embed_store.chunk_paths(out / "emb", MODEL, chunk_id)
    good = embed_store.load_chunk(out / "emb", MODEL, chunk_id)[1]["array_sha256"]
    np.save(npy, np.zeros_like(np.load(npy)))
    assert not embed_store.chunk_is_done(out / "emb", MODEL, chunk_id, "x", 1)
    again = run_cpu(work, tmp_path / "scratch")
    assert again["done"] == 1
    # determinism on CPU: the recomputed chunk is bit-identical to the first run
    assert embed_store.load_chunk(out / "emb", MODEL, chunk_id)[1]["array_sha256"] == good


def test_kill_and_resume_gives_the_same_arrays(tmp_path):
    work_a, work_b = tmp_path / "a", tmp_path / "b"
    out_a, _ = make_plan(work_a)
    out_b, _ = make_plan(work_b)
    run_cpu(work_b, tmp_path / "scratch")  # uninterrupted reference
    script = STEP1_DIR / "jobs" / "embed_chunks.py"
    env = {
        **os.environ,
        "PYTHONPATH": os.pathsep.join(
            [str(STEP1_DIR.parents[1] / "src"), str(STEP1_DIR), str(STEP1_DIR / "jobs")]
        ),
    }
    base = [sys.executable, str(script), "--work-dir", str(work_a), "--models", MODEL]
    base += ["--device", "cpu", "--batch-size", "2", "--scratch-dir", str(tmp_path / "s")]
    killed = subprocess.run(base + ["--stop-after-chunks", "1"], env=env, capture_output=True)
    assert killed.returncode == 3, killed.stderr.decode()[-2000:]
    resumed = subprocess.run(base, env=env, capture_output=True, text=True)
    assert resumed.returncode == 0, resumed.stderr[-2000:]
    assert '"skipped": 1' in resumed.stdout
    for c in chunk_plan.read_plan(out_a):
        a = embed_store.load_chunk(out_a / "emb", MODEL, c.chunk_id)[1]["array_sha256"]
        b = embed_store.load_chunk(out_b / "emb", MODEL, c.chunk_id)[1]["array_sha256"]
        assert a == b, c.chunk_id


def test_chunk_rows_do_not_depend_on_batch_size(tmp_path):
    seqs = ["MKTLLVAGLLSSAAFA", "MSTTSSTTSTPSSTSA" * 4, LONG[-1022:]]
    one = embed_chunks.embed_window_sequences(seqs, MODEL, 1, CPU, 6)
    three = embed_chunks.embed_window_sequences(seqs, MODEL, 3, CPU, 6)
    np.testing.assert_allclose(one, three, atol=1e-5)


def test_cterm_chunk_row_is_the_embedding_of_the_last_1022_residues(tmp_path):
    from surface_glyco.embeddings import get_esm_embeddings

    work = tmp_path / "w"
    out, unique = make_plan(work, chunk_residues=5000)
    run_cpu(work, tmp_path / "scratch")
    cterm = [c for c in chunk_plan.read_plan(out) if c.window == "cterm"]
    assert len(cterm) == 1 and len(cterm[0].rows) == 1
    arr, _ = embed_store.load_chunk(out / "emb", MODEL, cterm[0].chunk_id)
    direct, _ = get_esm_embeddings([{"id": "x", "sequence": LONG[-1022:]}], MODEL, 1, CPU)
    np.testing.assert_allclose(arr[0], direct[0], atol=1e-5)
    nterm_direct, _ = get_esm_embeddings([{"id": "x", "sequence": LONG}], MODEL, 1, CPU)
    assert not np.allclose(arr[0], nterm_direct[0], atol=1e-3)


def test_a_skipped_sequence_stops_the_chunk(monkeypatch):
    import surface_glyco.embeddings as emb

    def drop_last(records, **kw):
        n = len(records) - 1
        return np.zeros((n, 320), dtype=np.float32), [r["id"] for r in records[:n]], list(range(n))

    monkeypatch.setattr(emb, "get_esm_embeddings", drop_last)
    with pytest.raises(embed_chunks.EmbedError, match="not embedded"):
        embed_chunks.embed_window_sequences(["MKV", "MKL"], MODEL, 2, CPU, 6)


def test_nan_result_is_refused(tmp_path):
    arr = np.full((2, 4), np.nan, dtype=np.float32)
    with pytest.raises(ValueError, match="NaN or inf"):
        embed_store.save_chunk(tmp_path / "p", tmp_path / "s", MODEL, "nterm_0000", arr, {})
    assert not (tmp_path / "p").exists()


@pytest.mark.skipif(torch.cuda.is_available(), reason="needs a host without a GPU")
def test_cuda_request_without_gpu_stops(tmp_path, capsys):
    work = tmp_path / "w"
    make_plan(work)
    argv = ["--work-dir", str(work), "--device", "cuda", "--scratch-dir", str(tmp_path / "s")]
    assert embed_chunks.main(argv) == 2
    assert "cuda" in capsys.readouterr().err
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phaseb_embed.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'embed_chunks'`.

- [ ] **Step 5: Write `jobs/embed_store.py`**

```python
"""Chunk files of the Phase B embeddings: save, verify, load. Needs numpy.

Each chunk is two files in <publish>/<model>/: `<chunk_id>.npy` (float32, one row per chunk
member, in chunk order) and `<chunk_id>.json` (the done marker). The JSON is written last, so
a chunk without JSON is not done. Both files are first written to a temporary name and then
moved with os.replace, so a killed job never leaves a half-written chunk under its final name.
"""

import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np

DTYPE = np.float32


def array_sha256(arr: np.ndarray) -> str:
    """SHA-256 of the raw array bytes (C order). Shape and dtype are stored beside it."""
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def chunk_paths(publish: Path, model: str, chunk_id: str) -> tuple[Path, Path]:
    base = Path(publish) / model
    return base / f"{chunk_id}.npy", base / f"{chunk_id}.json"


def _atomic_copy(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".tmp.{dest.name}")
    try:
        shutil.copyfile(src, tmp)
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)


def save_chunk(publish: Path, scratch: Path, model: str, chunk_id: str, arr, meta: dict) -> dict:
    """Write the chunk to `scratch`, then copy it to `publish`. Return the stored metadata."""
    arr = np.ascontiguousarray(arr, dtype=DTYPE)
    if arr.ndim != 2:
        raise ValueError(f"{chunk_id}: expected a 2-D array, got shape {arr.shape}")
    if not np.isfinite(arr).all():
        raise ValueError(f"{chunk_id}: array has NaN or inf")
    meta = {
        **meta,
        "chunk_id": chunk_id,
        "model": model,
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "array_sha256": array_sha256(arr),
    }
    s_npy, s_json = chunk_paths(scratch, model, chunk_id)
    s_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(s_npy, arr)
    s_json.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    p_npy, p_json = chunk_paths(publish, model, chunk_id)
    _atomic_copy(s_npy, p_npy)
    _atomic_copy(s_json, p_json)
    return meta


def load_chunk(publish: Path, model: str, chunk_id: str):
    """Return (array, metadata) of a finished chunk. Raise ValueError if it fails a check."""
    npy, js = chunk_paths(publish, model, chunk_id)
    meta = json.loads(js.read_text())
    arr = np.load(npy, allow_pickle=False)
    if list(arr.shape) != meta["shape"] or str(arr.dtype) != meta["dtype"]:
        raise ValueError(f"{model}/{chunk_id}: shape or dtype differs from its JSON")
    if array_sha256(arr) != meta["array_sha256"]:
        raise ValueError(f"{model}/{chunk_id}: array SHA-256 differs from its JSON")
    return arr, meta


def chunk_is_done(publish: Path, model: str, chunk_id: str, members_sha256: str, n: int) -> bool:
    """True if the chunk exists, matches the plan (members, row count) and its hash."""
    npy, js = chunk_paths(publish, model, chunk_id)
    if not (npy.exists() and js.exists()):
        return False
    try:
        arr, meta = load_chunk(publish, model, chunk_id)
    except (OSError, ValueError, KeyError):
        return False
    return meta.get("members_sha256") == members_sha256 and arr.shape[0] == n
```

- [ ] **Step 6: Write `jobs/embed_chunks.py`**

```python
#!/usr/bin/env python
"""J2: embed the chunks of one job with ESM-2 (layer 6, residue-mean pooling). Resumable.

Run with the adhesionPred conda env and
PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/jobs.
Pooling is done by
surface_glyco.embeddings.get_esm_embeddings; this script only cuts windows, groups chunks,
checks the result and stores it. Reads $STEP1_WORKDIR/phaseb/unique_sequences.tsv.gz,
chunk_plan.tsv, chunk_members.tsv.gz and job_plan.json (batch size per model).

For each model and each chunk of this job (chunk k goes to job k mod --job-count):
- skip the chunk if phaseb/emb/<model>/<chunk_id>.npy and .json exist, the JSON names the
  same members_sha256, and the array SHA-256 matches the JSON;
- else embed it, write it under --scratch-dir, and copy it to phaseb/emb/<model>/.

STOP (exit 2): the plan files disagree; --device cuda without a GPU; a sequence that
get_esm_embeddings skipped; NaN or inf in a result. Exit 3: --stop-after-chunks was reached
(test hook that simulates a killed job).
"""

import argparse
import json
import sys
import time
from pathlib import Path

import chunk_plan
import embed_store
import numpy as np
import paths
import seqwindow
import truth_table

MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")
REPR_LAYER = 6


class EmbedError(RuntimeError):
    """A chunk cannot be embedded completely and correctly."""


def embed_window_sequences(seqs: list[str], model: str, batch_size: int, device, layer: int):
    """float32 array (len(seqs), dim) from get_esm_embeddings, in input order."""
    from surface_glyco.embeddings import get_esm_embeddings

    records = [{"id": str(i), "sequence": s} for i, s in enumerate(seqs)]
    arr, _, kept = get_esm_embeddings(
        records,
        model_name=model,
        batch_size=batch_size,
        device=device,
        repr_layer=layer,
        return_indices=True,
    )
    if list(kept) != list(range(len(seqs))):
        missing = sorted(set(range(len(seqs))) - set(kept))
        raise EmbedError(f"{len(missing)} sequences were not embedded (positions {missing[:5]})")
    return np.asarray(arr, dtype=np.float32)


def run_job(
    work: Path,
    models,
    job_index: int,
    job_count: int,
    device,
    scratch: Path,
    batch_size: int | None = None,
    stop_after: int | None = None,
) -> dict:
    out = Path(work) / "phaseb"
    publish = out / "emb"
    chunks = {c.chunk_id: c for c in chunk_plan.read_plan(out)}
    mine = chunk_plan.chunks_for_job(list(chunks), job_index, job_count)
    job_plan = json.loads((out / "job_plan.json").read_text())
    seq_by_row = {
        int(r["row"]): r["sequence"] for r in truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    }
    done = skipped = 0
    for model in models:
        size = batch_size or int(job_plan["batch_size"][model])
        for chunk_id in mine:
            chunk = chunks[chunk_id]
            if embed_store.chunk_is_done(
                publish, model, chunk_id, chunk.members_sha256, len(chunk.rows)
            ):
                skipped += 1
                print(f"skip {model}/{chunk_id} (done, hash verified)", flush=True)
                continue
            seqs = [seqwindow.window(seq_by_row[r], chunk.window) for r in chunk.rows]
            start = time.time()
            arr = embed_window_sequences(seqs, model, size, device, REPR_LAYER)
            meta = {
                "window": chunk.window,
                "members_sha256": chunk.members_sha256,
                "repr_layer": REPR_LAYER,
                "batch_size": size,
                "device": str(device),
                "seconds": round(time.time() - start, 2),
            }
            embed_store.save_chunk(publish, scratch, model, chunk_id, arr, meta)
            done += 1
            print(f"done {model}/{chunk_id} n={len(seqs)} s={meta['seconds']}", flush=True)
            if stop_after is not None and done >= stop_after:
                return {"done": done, "skipped": skipped, "stopped": True}
    return {"done": done, "skipped": skipped, "stopped": False}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    parser.add_argument("--job-index", type=int, default=0)
    parser.add_argument("--job-count", type=int, default=1)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--scratch-dir", required=True, help="node-local, e.g. $SCRATCH/emb")
    parser.add_argument("--batch-size", type=int, default=None, help="default: job_plan.json")
    parser.add_argument("--stop-after-chunks", type=int, default=None, help="test hook")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        import torch

        if args.device == "cuda" and not torch.cuda.is_available():
            raise EmbedError("--device cuda but torch.cuda.is_available() is False")
        result = run_job(
            work,
            args.models,
            args.job_index,
            args.job_count,
            torch.device(args.device),
            Path(args.scratch_dir),
            args.batch_size,
            args.stop_after_chunks,
        )
    except (EmbedError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result), flush=True)
    return 3 if result["stopped"] else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phaseb_embed.py -q`
Expected: `9 passed` (the prototype ran this file and the assembly tests, 11 tests, in 80 s on 2 CPU cores; ESM-2 8M on CPU). The determinism test asserts a bit-identical recomputed chunk on CPU. The kill test runs the CLI in a subprocess, stops it after one chunk (exit 3), resumes it, and compares every chunk hash with an uninterrupted run.
Run: `$PY -m pytest tests/step1_compare -q`
Expected: no failure; `test_phaseb_embed.py` is skipped (no torch), `test_every_module_imports_only_stdlib_or_local` passes (it scans the top folder only, not `jobs/`).

- [ ] **Step 8: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/jobs/embed_store.py analysis/step1_compare/jobs/embed_chunks.py \
  tests/step1_compare/conftest.py tests/step1_compare/phaseb_fixture.py \
  tests/step1_compare/test_phaseb_embed.py
git commit -m "step1_compare: resumable J2 embedding runner with hashed chunk files"
```

### Task 5: Assembly and per-chunk manifest

**Files:**
- Create: `analysis/step1_compare/jobs/assemble_embeddings.py`
- Test: `tests/step1_compare/test_phaseb_assemble.py`

**Interfaces:**
- Consumes: `embed_store.load_chunk`, `embed_store.array_sha256`, `embed_store.DTYPE` (Task 4); `chunk_plan.read_plan` (Task 3); `phaseb_fixture.make_plan`, `run_cpu` (tests, Task 4).
- Produces: `assemble_embeddings.MODELS`, `MANIFEST_COLUMNS`, `AssembleError(ValueError)`, `window_matrix(nterm, cterm, cterm_rows: list[str]) -> np.ndarray`, `assemble_model(publish, model, chunks, unique, manifest_rows) -> dict[str, np.ndarray]`, `run(work, models) -> dict`, `main(argv=None)`. Outputs in `phaseb/emb/`: `<model>.nterm.npy`, `<model>.cterm.npy`, `chunk_manifest.tsv`, `embedding_run.json`.

- [ ] **Step 1: Write the failing test**

```python
"""Assembly of the J2 chunks into matrices (CPU, ESM-2 8M)."""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("torch")
pytest.importorskip("esm")

import assemble_embeddings  # noqa: E402
import chunk_plan  # noqa: E402
import embed_store  # noqa: E402
import truth_table  # noqa: E402
from phaseb_fixture import MODEL, make_plan, run_cpu  # noqa: E402


def test_assemble_builds_row_ordered_matrices(tmp_path):
    work = tmp_path / "w"
    out, unique = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    assert assemble_embeddings.main(["--work-dir", str(work), "--models", MODEL]) == 0
    nterm = np.load(out / "emb" / f"{MODEL}.nterm.npy")
    cterm = np.load(out / "emb" / f"{MODEL}.cterm.npy")
    assert nterm.shape == (len(unique), 320) and nterm.dtype == np.float32
    assert cterm.shape == (1, 320)
    for c in chunk_plan.read_plan(out):
        arr, _ = embed_store.load_chunk(out / "emb", MODEL, c.chunk_id)
        for pos, row in enumerate(c.rows):
            target = nterm[row] if c.window == "nterm" else cterm[int(unique[row]["cterm_row"])]
            np.testing.assert_array_equal(target, arr[pos])
    manifest = truth_table.read_tsv(out / "emb" / "chunk_manifest.tsv")
    assert [m["chunk_id"] for m in manifest] == [c.chunk_id for c in chunk_plan.read_plan(out)]
    for m in manifest:
        meta = embed_store.load_chunk(out / "emb", MODEL, m["chunk_id"])[1]
        assert m["array_sha256"] == meta["array_sha256"]
    full_c = assemble_embeddings.window_matrix(nterm, cterm, [u["cterm_row"] for u in unique])
    long_row = next(int(u["row"]) for u in unique if u["cterm_row"] != "")
    np.testing.assert_array_equal(full_c[long_row], cterm[0])
    short_rows = [int(u["row"]) for u in unique if u["cterm_row"] == ""]
    np.testing.assert_array_equal(full_c[short_rows], nterm[short_rows])


def test_assemble_stops_on_a_missing_chunk(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk_id = chunk_plan.read_plan(out)[-1].chunk_id
    embed_store.chunk_paths(out / "emb", MODEL, chunk_id)[1].unlink()
    assert assemble_embeddings.main(["--work-dir", str(work), "--models", MODEL]) == 2
    assert chunk_id in capsys.readouterr().err
    assert not (out / "emb" / f"{MODEL}.nterm.npy").exists()
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phaseb_assemble.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'assemble_embeddings'`.

- [ ] **Step 3: Write `jobs/assemble_embeddings.py`**

```python
#!/usr/bin/env python
"""After J2: check every chunk and build one matrix per model and window. Needs numpy.

Reads phaseb/unique_sequences.tsv.gz, the chunk plan and phaseb/emb/<model>/<chunk_id>.*.
Writes to phaseb/emb/:

  <model>.nterm.npy        float32, shape (unique sequences, dim); row = `row` column
  <model>.cterm.npy        float32, shape (sequences > 1,022 aa, dim); row = `cterm_row`
  chunk_manifest.tsv       one row per model and chunk: members and array SHA-256
  embedding_run.json       shapes, matrix SHA-256 values, chunk count, plan inputs

A sequence of 1,022 aa or less has no C-terminal row: its C-terminal window is the whole
sequence, so the M8-C and M35-C candidates use its `nterm` row (see `window_matrix`).
STOP (exit 2, no output): a chunk is missing, fails its hash, or holds other members than the
plan; a matrix row is filled twice or not at all; NaN or inf; a row count that is not the
unique sequence count.
"""

import argparse
import sys
from pathlib import Path

import chunk_plan
import embed_store
import manifest
import numpy as np
import paths
import runinfo
import truth_table

MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")
MANIFEST_COLUMNS = ("model", "chunk_id", "window", "n_seqs", "members_sha256", "array_sha256")


class AssembleError(ValueError):
    """The chunks do not form complete matrices."""


def window_matrix(nterm: np.ndarray, cterm: np.ndarray, cterm_rows: list[str]) -> np.ndarray:
    """Full C-terminal matrix (one row per unique sequence) for the M8-C / M35-C candidates."""
    out = nterm.copy()
    for row, crow in enumerate(cterm_rows):
        if crow != "":
            out[row] = cterm[int(crow)]
    return out


def assemble_model(publish: Path, model: str, chunks, unique, manifest_rows: list) -> dict:
    n_rows = {"nterm": len(unique), "cterm": sum(r["cterm_row"] != "" for r in unique)}
    target = {int(r["row"]): int(r["cterm_row"]) if r["cterm_row"] != "" else None for r in unique}
    mats: dict[str, np.ndarray] = {}
    filled = {w: np.zeros(n, dtype=bool) for w, n in n_rows.items()}
    for chunk in chunks:
        try:
            arr, meta = embed_store.load_chunk(publish, model, chunk.chunk_id)
        except (OSError, ValueError, KeyError) as exc:
            raise AssembleError(f"{model}/{chunk.chunk_id}: {exc}") from exc
        if meta.get("members_sha256") != chunk.members_sha256 or arr.shape[0] != len(chunk.rows):
            raise AssembleError(f"{model}/{chunk.chunk_id}: members differ from the plan")
        manifest_rows.append(
            {
                "model": model,
                "chunk_id": chunk.chunk_id,
                "window": chunk.window,
                "n_seqs": str(len(chunk.rows)),
                "members_sha256": chunk.members_sha256,
                "array_sha256": meta["array_sha256"],
            }
        )
        w = chunk.window
        if w not in mats:
            mats[w] = np.zeros((n_rows[w], arr.shape[1]), dtype=embed_store.DTYPE)
        idx = [r if w == "nterm" else target[r] for r in chunk.rows]
        if any(i is None for i in idx) or filled[w][idx].any():
            raise AssembleError(f"{model}/{chunk.chunk_id}: a row is filled twice or has no slot")
        mats[w][idx] = arr
        filled[w][idx] = True
    for w, n in n_rows.items():
        if n and (w not in mats or not filled[w].all()):
            raise AssembleError(f"{model}: {int((~filled[w]).sum())} {w} rows have no embedding")
        if w not in mats:
            mats[w] = np.zeros((0, mats["nterm"].shape[1]), dtype=embed_store.DTYPE)
        if not np.isfinite(mats[w]).all():
            raise AssembleError(f"{model}: {w} matrix has NaN or inf")
    return mats


def run(work: Path, models) -> dict:
    out = Path(work) / "phaseb"
    publish = out / "emb"
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    chunks = chunk_plan.read_plan(out)
    log = {"models": {}, "unique_sequences": len(unique), "chunks": len(chunks)}
    writers, manifest_rows = {}, []
    for model in models:
        mats = assemble_model(publish, model, chunks, unique, manifest_rows)
        log["models"][model] = {
            w: {
                "shape": list(m.shape),
                "dtype": str(m.dtype),
                "array_sha256": embed_store.array_sha256(m),
            }
            for w, m in mats.items()
        }
        for w, m in mats.items():
            writers[f"{model}.{w}.npy"] = lambda p, m=m: _save(p, m)
    log.update(
        {
            "unique_sequences_sha256": _sha(out / "unique_sequences.tsv.gz"),
            "chunk_plan_sha256": _sha(out / "chunk_plan.tsv"),
            "git_commit": runinfo.git_commit(),
            "python": runinfo.python_version(),
        }
    )
    writers["chunk_manifest.tsv"] = lambda p: truth_table.write_tsv(
        p, MANIFEST_COLUMNS, manifest_rows
    )
    writers["embedding_run.json"] = lambda p: runinfo.write_json(p, log)
    runinfo.atomic_write_all(publish, writers)
    return log


def _save(path: Path, arr: np.ndarray) -> None:
    with open(path, "wb") as handle:  # np.save on a handle keeps the temp name unchanged
        np.save(handle, arr)


def _sha(path: Path) -> str:
    return manifest.sha256_file(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        log = run(work, args.models)
    except (AssembleError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for model, mats in log["models"].items():
        print(model, {w: m["shape"] for w, m in mats.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phaseb_assemble.py -q`
Expected: `2 passed`.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/jobs/assemble_embeddings.py tests/step1_compare/test_phaseb_assemble.py
git commit -m "step1_compare: assemble J2 chunks into checked matrices with a chunk manifest"
```

### Task 6: J0 throughput pilot and GPU-CPU difference harness

**Files:**
- Create: `analysis/step1_compare/jobs/throughput_pilot.py`, `analysis/step1_compare/jobs/gpu_cpu_diff.py`
- Test: `tests/step1_compare/test_phaseb_pilot.py`

**Interfaces:**
- Consumes: `seqwindow.nterm_window`, `cterm_window`, `needs_cterm` (Task 1); `chunk_plan.members_digest` (Task 3); `06_plan_embedding.main` (test, Task 3); `seqsets.unique_rows` (test, Task 2); `surface_glyco.embeddings.get_cached_model`, `get_esm_embeddings`.
- Produces: `throughput_pilot.SCHEMA = "step1-phaseb-j0-throughput/1"`, `PilotError(RuntimeError)`, `sample_sequences(truth_rows, n, seed) -> list[tuple[str, str]]`, `time_model(model, sample, batch_sizes, device) -> (runs, load_s)`, `best_runs(runs) -> dict`, `run(truth_path, out_path, models, batch_sizes, n, seed, device) -> dict`, `main(argv=None)`; JSON keys used by 06: `best[model].batch_size`, `best[model].residues_per_s`, `model_load_s[model]`. `gpu_cpu_diff.SCHEMA = "step1-phaseb-gpu-cpu-diff/1"`, `windows_for(truth_rows, n, n_long, seed) -> list[str]`, `embed(seqs, model, device, batch_size)`, `compare(a, b) -> dict` (keys `max_abs_diff`, `mean_abs_diff`, `max_rel_diff`, `min_cosine`), `run(...)`, `main(argv=None)`. Outputs: `phaseb/j0/throughput.json`, `phaseb/j0/gpu_cpu_diff.json`.

- [ ] **Step 1: Write the failing test**

```python
"""J0 pilot and the GPU-CPU difference harness, run on CPU with 20 sequences (ESM-2 8M)."""

import json

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("esm")

import gpu_cpu_diff  # noqa: E402
import seqhash  # noqa: E402
import throughput_pilot  # noqa: E402
import truth_table  # noqa: E402

MODEL = "esm2_t6_8M_UR50D"
COLUMNS = ("source_id", "gene_id", "label", "fasta_id", "length", "seq_sha256", "sequence")


def _truth(tmp_path, n=30):
    rows = []
    for i in range(n):
        seq = ("MKTLLVAGLLSSAAFA" + "ACDEFGHIKLMNPQRSTVWY"[i % 20] * (5 + 7 * i))[: 40 + 9 * i]
        if i == 0:
            seq = "MKLSTA" + "STPSSTSA" * 140  # one protein longer than 1,022 aa
        rows.append(
            {
                "source_id": "Scer_SGD",
                "gene_id": f"G{i}",
                "label": "N-int",
                "fasta_id": f"G{i}",
                "length": str(len(seq)),
                "seq_sha256": seqhash.seq_sha256(seq),
                "sequence": seq,
            }
        )
    rows.append({**rows[1], "source_id": "Calb_CGD", "gene_id": "C1"})  # duplicate sequence
    path = tmp_path / "truth_sequences.tsv.gz"
    truth_table.write_tsv(path, COLUMNS, rows)
    return path, rows


def test_sample_is_unique_fixed_by_seed_and_capped(tmp_path):
    _, rows = _truth(tmp_path)
    a = throughput_pilot.sample_sequences(rows, 20, 7)
    b = throughput_pilot.sample_sequences(list(reversed(rows)), 20, 7)
    assert a == b
    assert len({k for k, _ in a}) == 20
    assert max(len(s) for _, s in a) <= 1022
    assert throughput_pilot.sample_sequences(rows, 20, 8) != a
    with pytest.raises(throughput_pilot.PilotError):
        throughput_pilot.sample_sequences(rows, 31, 7)  # 30 unique sequences only


def test_pilot_writes_the_throughput_schema_on_cpu(tmp_path):
    truth_path, _ = _truth(tmp_path)
    out = tmp_path / "phaseb" / "j0" / "throughput.json"
    argv = ["--work-dir", str(tmp_path), "--device", "cpu", "--models", MODEL]
    assert throughput_pilot.main(argv + ["--n", "20", "--batch-sizes", "2,4"]) == 0
    rec = json.loads(out.read_text())
    assert rec["schema"] == throughput_pilot.SCHEMA
    assert rec["device"] == "cpu" and rec["gpu_name"] == "cpu"
    assert rec["n_proteins"] == 20 and rec["seed"] == 20261001
    assert len(rec["sample_sha256"]) == 64 and len(rec["truth_sequences_sha256"]) == 64
    assert [r["batch_size"] for r in rec["runs"]] == [2, 4]
    for r in rec["runs"]:
        assert r["status"] == "ok" and r["dim"] == 320
        assert r["proteins_per_s"] > 0 and r["residues_per_s"] > 0
        assert r["peak_mem_bytes"] is None
    assert rec["best"][MODEL]["batch_size"] in (2, 4)
    assert rec["model_load_s"][MODEL] >= 0
    for key in ("torch", "torch_cuda", "esm", "git_commit", "python", "residues"):
        assert key in rec


def test_pilot_output_feeds_the_plan_script(tmp_path):
    from conftest import load_script

    truth_path, rows = _truth(tmp_path)
    argv = ["--work-dir", str(tmp_path), "--device", "cpu", "--models", MODEL]
    assert throughput_pilot.main(argv + ["--n", "20", "--batch-sizes", "4"]) == 0
    import seqsets

    unique = seqsets.unique_rows({r["seq_sha256"]: r["sequence"] for r in rows})
    out = tmp_path / "phaseb"
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    plan = load_script("06_plan_embedding")
    assert plan.main(["--work-dir", str(tmp_path), "--models", MODEL]) == 0
    job = json.loads((out / "job_plan.json").read_text())
    assert job["rate_source"] == "J0" and job["batch_size"][MODEL] == 4


def test_gpu_cpu_diff_schema_with_cpu_on_both_sides(tmp_path):
    _truth(tmp_path)
    out = tmp_path / "diff.json"
    argv = ["--work-dir", str(tmp_path), "--out", str(out), "--models", MODEL]
    argv += ["--device-a", "cpu", "--device-b", "cpu", "--n", "10", "--n-long", "1"]
    assert gpu_cpu_diff.main(argv) == 0
    rec = json.loads(out.read_text())
    assert rec["schema"] == gpu_cpu_diff.SCHEMA and rec["n_windows"] == 11
    m = rec["models"][MODEL]
    assert m["max_abs_diff"] == 0.0 and m["repeat_identical_on_a"] is True
    assert m["min_cosine"] == pytest.approx(1.0, abs=1e-6)


@pytest.mark.skipif(torch.cuda.is_available(), reason="needs a host without a GPU")
def test_gpu_cpu_diff_stops_without_a_gpu(tmp_path, capsys):
    _truth(tmp_path)
    assert gpu_cpu_diff.main(["--work-dir", str(tmp_path), "--models", MODEL]) == 2
    assert "cuda" in capsys.readouterr().err
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phaseb_pilot.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'gpu_cpu_diff'`.

- [ ] **Step 3: Write `jobs/throughput_pilot.py`**

```python
#!/usr/bin/env python
"""J0: measure ESM-2 embedding throughput on a fixed sample of truth sequences.

Samples --n unique sequences (by seq_sha256) from truth_sequences.tsv.gz with a fixed seed,
cuts the N-terminal window (first 1,022 residues), and times
surface_glyco.embeddings.get_esm_embeddings for every model and batch size. Writes one JSON
file (schema SCHEMA) with proteins/s, residues/s, peak GPU memory, the GPU model and, per
model, the fastest batch size (`best`). 06_plan_embedding.py reads `best` and `model_load_s`.

A batch size that runs out of GPU memory is recorded with status `oom`. STOP (exit 2): the
input is missing, --device cuda without a GPU, or no model has a successful run.
"""

import argparse
import json
import random
import sys
import time
from pathlib import Path

import chunk_plan
import manifest
import paths
import runinfo
import seqwindow
import truth_table

SCHEMA = "step1-phaseb-j0-throughput/1"
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


class PilotError(RuntimeError):
    """The pilot cannot produce a usable throughput record."""


def sample_sequences(truth_rows, n: int, seed: int) -> list[tuple[str, str]]:
    """(seq_sha256, N-terminal window) for n unique sequences; same input and seed, same set."""
    unique = {r["seq_sha256"]: r["sequence"] for r in truth_rows}
    keys = sorted(unique)
    if n > len(keys):
        raise PilotError(f"--n {n} is larger than the {len(keys)} unique truth sequences")
    chosen = sorted(random.Random(seed).sample(keys, n))
    return [(k, seqwindow.nterm_window(unique[k])) for k in chosen]


def _sync(torch, device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize()


def time_model(model: str, sample, batch_sizes, device) -> tuple[list[dict], float]:
    import torch

    from surface_glyco import embeddings

    start = time.time()
    embeddings.get_cached_model(model, device)
    _sync(torch, device)
    load_s = time.time() - start
    records = [{"id": k, "sequence": s} for k, s in sample]
    residues = sum(len(s) for _, s in sample)
    embeddings.get_esm_embeddings(records[:4], model, 4, device)  # warm-up, not timed
    runs = []
    for size in batch_sizes:
        if device.type == "cuda":
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats(device)
        run = {"model": model, "batch_size": size}
        try:
            _sync(torch, device)
            t0 = time.time()
            arr, _, kept = embeddings.get_esm_embeddings(
                records, model, size, device, return_indices=True
            )
            _sync(torch, device)
            secs = time.time() - t0
        except torch.cuda.OutOfMemoryError:
            runs.append({**run, "status": "oom"})
            continue
        if len(kept) != len(records):
            runs.append({**run, "status": f"skipped_{len(records) - len(kept)}"})
            continue
        run.update(
            {
                "status": "ok",
                "seconds": round(secs, 3),
                "proteins_per_s": round(len(records) / secs, 2),
                "residues_per_s": round(residues / secs, 1),
                "peak_mem_bytes": int(torch.cuda.max_memory_allocated(device))
                if device.type == "cuda"
                else None,
                "dim": int(arr.shape[1]),
            }
        )
        runs.append(run)
    return runs, round(load_s, 3)


def best_runs(runs) -> dict:
    best = {}
    for r in runs:
        if r["status"] != "ok":
            continue
        cur = best.get(r["model"])
        if cur is None or r["residues_per_s"] > cur["residues_per_s"]:
            best[r["model"]] = {
                k: r[k] for k in ("batch_size", "residues_per_s", "proteins_per_s", "seconds")
            }
    return best


def run(truth_path: Path, out_path: Path, models, batch_sizes, n, seed, device) -> dict:
    import esm
    import torch

    sample = sample_sequences(truth_table.read_tsv(truth_path), n, seed)
    runs, load_s = [], {}
    for model in models:
        got, load_s[model] = time_model(model, sample, batch_sizes, device)
        runs += got
    best = best_runs(runs)
    if not best:
        raise PilotError("no model has a successful run")
    record = {
        "schema": SCHEMA,
        "device": device.type,
        "gpu_name": torch.cuda.get_device_name(device) if device.type == "cuda" else "cpu",
        "torch": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "esm": getattr(esm, "__version__", "unknown"),
        "n_proteins": len(sample),
        "residues": sum(len(s) for _, s in sample),
        "seed": seed,
        "sample_sha256": chunk_plan.members_digest([k for k, _ in sample]),
        "truth_sequences_sha256": manifest.sha256_file(truth_path),
        "batch_sizes": list(batch_sizes),
        "model_load_s": load_s,
        "runs": runs,
        "best": best,
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writers = {out_path.name: lambda p: runinfo.write_json(p, record)}
    runinfo.atomic_write_all(out_path.parent, writers)
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--out", default=None, help="default: phaseb/j0/throughput.json")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    parser.add_argument("--batch-sizes", default="8,16,32,64")
    parser.add_argument("--n", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    out = Path(args.out) if args.out else work / "phaseb" / "j0" / "throughput.json"
    try:
        import torch

        if args.device == "cuda" and not torch.cuda.is_available():
            raise PilotError("--device cuda but torch.cuda.is_available() is False")
        sizes = [int(x) for x in args.batch_sizes.split(",")]
        record = run(
            work / "truth_sequences.tsv.gz",
            out,
            args.models,
            sizes,
            args.n,
            args.seed,
            torch.device(args.device),
        )
    except (PilotError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(record["best"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Write `jobs/gpu_cpu_diff.py`**

```python
#!/usr/bin/env python
"""Spec 7 harness: difference between ESM-2 embeddings computed on two devices.

Samples --n unique truth sequences (fixed seed, the J0 sampler), adds the C-terminal window of
the first --n-long sequences longer than 1,022 aa, and embeds all windows with
surface_glyco.embeddings.get_esm_embeddings on --device-a and --device-b (default cuda and
cpu). It also embeds the windows twice on --device-a to check run-to-run repeatability.
Writes one JSON file per call: per model the maximum and mean absolute difference, the largest
difference relative to the largest absolute value, the lowest cosine similarity, and whether
the repeat on device A is bit-identical. No threshold is applied: the numbers are reported.

STOP (exit 2): a cuda device is requested and none is available; a sequence was not embedded.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import chunk_plan
import numpy as np
import paths
import runinfo
import seqwindow
import throughput_pilot
import truth_table

SCHEMA = "step1-phaseb-gpu-cpu-diff/1"
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


def windows_for(truth_rows, n: int, n_long: int, seed: int) -> list[str]:
    sample = throughput_pilot.sample_sequences(truth_rows, n, seed)
    unique = {r["seq_sha256"]: r["sequence"] for r in truth_rows}
    long_keys = sorted(k for k, s in unique.items() if seqwindow.needs_cterm(s))[:n_long]
    return [s for _, s in sample] + [seqwindow.cterm_window(unique[k]) for k in long_keys]


def embed(seqs, model, device, batch_size) -> np.ndarray:
    from surface_glyco.embeddings import get_esm_embeddings

    records = [{"id": str(i), "sequence": s} for i, s in enumerate(seqs)]
    arr, _, kept = get_esm_embeddings(records, model, batch_size, device, return_indices=True)
    if len(kept) != len(seqs):
        raise RuntimeError(f"{model} on {device}: {len(seqs) - len(kept)} sequences skipped")
    return np.asarray(arr, dtype=np.float32)


def compare(a: np.ndarray, b: np.ndarray) -> dict:
    diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
    na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    cos = (a * b).sum(axis=1) / np.maximum(na * nb, 1e-12)
    return {
        "max_abs_diff": float(diff.max()),
        "mean_abs_diff": float(diff.mean()),
        "max_rel_diff": float(diff.max() / max(float(np.abs(b).max()), 1e-12)),
        "min_cosine": float(cos.min()),
    }


def run(truth_path, out_path, models, dev_a, dev_b, n, n_long, seed, batch_size) -> dict:
    import torch

    seqs = windows_for(truth_table.read_tsv(truth_path), n, n_long, seed)
    result = {
        "schema": SCHEMA,
        "device_a": str(dev_a),
        "device_b": str(dev_b),
        "gpu_name": torch.cuda.get_device_name(dev_a) if dev_a.type == "cuda" else "cpu",
        "torch": torch.__version__,
        "n_windows": len(seqs),
        "windows_sha256": chunk_plan.members_digest(
            [hashlib.sha256(s.encode()).hexdigest() for s in seqs]
        ),
        "batch_size": batch_size,
        "models": {},
        "git_commit": runinfo.git_commit(),
    }
    for model in models:
        a1 = embed(seqs, model, dev_a, batch_size)
        a2 = embed(seqs, model, dev_a, batch_size)
        b = embed(seqs, model, dev_b, batch_size)
        result["models"][model] = {
            **compare(a1, b),
            "repeat_identical_on_a": bool(np.array_equal(a1, a2)),
            "repeat_max_abs_diff_on_a": float(np.abs(a1 - a2).max()),
        }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writers = {out_path.name: lambda p: runinfo.write_json(p, result)}
    runinfo.atomic_write_all(out_path.parent, writers)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--out", default=None, help="default: phaseb/j0/gpu_cpu_diff.json")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    parser.add_argument("--device-a", default="cuda")
    parser.add_argument("--device-b", default="cpu")
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--n-long", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    out = Path(args.out) if args.out else work / "phaseb" / "j0" / "gpu_cpu_diff.json"
    try:
        import torch

        devices = [torch.device(args.device_a), torch.device(args.device_b)]
        if any(d.type == "cuda" for d in devices) and not torch.cuda.is_available():
            raise RuntimeError("a cuda device was requested but torch.cuda.is_available() is False")
        result = run(
            work / "truth_sequences.tsv.gz",
            out,
            args.models,
            devices[0],
            devices[1],
            args.n,
            args.n_long,
            args.seed,
            args.batch_size,
        )
    except (RuntimeError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result["models"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phaseb_pilot.py -q`
Expected: `5 passed` (about 45 s on CPU). This is the CPU acceptance run of the J0 schema on 20 sequences.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/jobs/throughput_pilot.py analysis/step1_compare/jobs/gpu_cpu_diff.py \
  tests/step1_compare/test_phaseb_pilot.py
git commit -m "step1_compare: J0 throughput pilot and GPU-CPU difference harness"
```

### Task 7: SignalP and PredGPI parsers, PredGPI wrapper

**Files:**
- Create: `analysis/step1_compare/feature_parsers.py`, `analysis/step1_compare/jobs/predgpi_scores.py`, `tests/step1_compare/fixtures/phaseb/signalp_prediction_results.txt`, `tests/step1_compare/fixtures/phaseb/signalp_output.gff3`, `tests/step1_compare/fixtures/phaseb/predgpi_scores.tsv`
- Test: `tests/step1_compare/test_feature_parsers.py`

**Fixture provenance.** The two SignalP fixtures are four real rows of `analysis/cocci_repeats/signalp/CimmitisRS_FungiDB/prediction_results.txt` and the two matching `output.gff3` lines (SignalP 6.0i GPU build, job 29280458), with the long FASTA headers replaced by 64-character ids. The PredGPI fixture holds four real `predgpi_scores.py` rows from the 1,000-protein S288C run (YAL069W, YAL068W-A, YAR050W, YAL037C-A), with ids replaced. Nothing in the fixtures is invented except the ids.

**Interfaces:**
- Consumes: `gaf.open_text` (Phase A).
- Produces: `feature_parsers.SP_COLUMNS`, `GPI_COLUMNS`, `GPI_CALLS`, `OutputFormatError(ValueError)`, dataclass `SignalPCall(prediction, other_prob, sp_prob, cs_end, cs_prob)`, `parse_signalp(path) -> dict[str, SignalPCall]`, `parse_signalp_gff(path) -> dict[str, int]`, `check_signalp_consistency(calls, gff_ends) -> None`, dataclass `GpiCall(call, prob, omega, fpr, svm)`, `parse_predgpi_scores(path) -> dict[str, GpiCall]`, `ser_thr_fraction(sequence) -> float`. `jobs/predgpi_scores.py --fasta F --out O` writes columns `id, length, gpi_call, gpi_prob, omega, fpr, svm`.

- [ ] **Step 1: Write the fixtures**

```bash
/usr/bin/python3.12 - <<'EOF'
from pathlib import Path

rows = [
    ('# SignalP-6.0', 'Organism: Eukarya', 'Timestamp: 20260930100126'),
    ('# ID', 'Prediction', 'OTHER', 'SP(Sec/SPI)', 'CS Position'),
    ('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa', 'OTHER', '1.000000', '0.000000', ''),
    ('bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb', 'OTHER', '1.000000', '0.000014', ''),
    ('cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc', 'SP', '0.000293', '0.999663', 'CS pos: 24-25. Pr: 0.5487'),
    ('dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd', 'SP', '0.057443', '0.942539', 'CS pos: 22-23. Pr: 0.7016'),
]
path = Path("tests/step1_compare/fixtures/phaseb/signalp_prediction_results.txt")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text("".join("\t".join(r) + "\n" for r in rows))
EOF
```

```bash
/usr/bin/python3.12 - <<'EOF'
from pathlib import Path

rows = [
    ('## gff-version 3',),
    ('cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc', 'SignalP-6.0', 'signal_peptide', '1', '24', '0.9996627', '.', '.', '.'),
    ('dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd', 'SignalP-6.0', 'signal_peptide', '1', '22', '0.9425392', '.', '.', '.'),
]
path = Path("tests/step1_compare/fixtures/phaseb/signalp_output.gff3")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text("".join("\t".join(r) + "\n" for r in rows))
EOF
```

```bash
/usr/bin/python3.12 - <<'EOF'
from pathlib import Path

rows = [
    ('id', 'length', 'gpi_call', 'gpi_prob', 'omega', 'fpr', 'svm'),
    ('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa', '104', 'weakly', '0.55', '80', '0.0095988', '-0.708166'),
    ('bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb', '84', 'none', '0', '', '0.878876', '-1.68145'),
    ('cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc', '1537', 'highly_probable', '1.0', '1512', '2.07693e-08', '1.21144'),
    ('dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd', '30', 'too_short', '0', '', '', ''),
]
path = Path("tests/step1_compare/fixtures/phaseb/predgpi_scores.tsv")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text("".join("\t".join(r) + "\n" for r in rows))
EOF
```

- [ ] **Step 2: Write the failing test**

```python
import gzip
import os
import shutil
import subprocess

import feature_parsers as fp
import pytest

A, B, C, D = "a" * 64, "b" * 64, "c" * 64, "d" * 64


@pytest.fixture
def phaseb_fixtures(fixtures_dir):
    return fixtures_dir / "phaseb"


def test_signalp_rows_parse_by_header_name(phaseb_fixtures):
    calls = fp.parse_signalp(phaseb_fixtures / "signalp_prediction_results.txt")
    assert set(calls) == {A, B, C, D}
    assert calls[A] == fp.SignalPCall("OTHER", 1.0, 0.0, None, None)
    assert calls[B].sp_prob == pytest.approx(0.000014)
    assert calls[C] == fp.SignalPCall("SP", 0.000293, 0.999663, 24, 0.5487)
    assert calls[D].cs_end == 22 and calls[D].cs_prob == pytest.approx(0.7016)


def test_signalp_gff_agrees_with_prediction_results(phaseb_fixtures):
    calls = fp.parse_signalp(phaseb_fixtures / "signalp_prediction_results.txt")
    ends = fp.parse_signalp_gff(phaseb_fixtures / "signalp_output.gff3")
    assert ends == {C: 24, D: 22}
    fp.check_signalp_consistency(calls, ends)
    with pytest.raises(fp.OutputFormatError, match="disagree"):
        fp.check_signalp_consistency(calls, {C: 25, D: 22})


def test_signalp_reads_gzip(phaseb_fixtures, tmp_path):
    src = phaseb_fixtures / "signalp_prediction_results.txt"
    gz = tmp_path / "prediction_results.txt.gz"
    gz.write_bytes(gzip.compress(src.read_bytes()))
    assert fp.parse_signalp(gz) == fp.parse_signalp(src)


@pytest.mark.parametrize(
    "line, message",
    [
        (f"{A}\tSP\t0.1\t0.9\t\n", "SP row with CS field"),
        (f"{A}\tOTHER\t0.9\t0.1\tCS pos: 3-4. Pr: 0.2\n", "OTHER row with a CS field"),
        (f"{A}\tLIPO\t0.1\t0.9\t\n", "unexpected prediction"),
        (f"{A} OTHER 1.0 0.0\n", "fields"),
    ],
)
def test_signalp_malformed_rows_stop(tmp_path, line, message):
    path = tmp_path / "p.txt"
    path.write_text("# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n" + line)
    with pytest.raises(fp.OutputFormatError, match=message):
        fp.parse_signalp(path)


def test_signalp_missing_header_or_duplicate_id_stops(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(f"{A}\tOTHER\t1.0\t0.0\t\n")
    with pytest.raises(fp.OutputFormatError, match="before the '# ID' header"):
        fp.parse_signalp(path)
    path.write_text(
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        f"{A}\tOTHER\t1.0\t0.0\t\n{A}\tOTHER\t1.0\t0.0\t\n"
    )
    with pytest.raises(fp.OutputFormatError, match="duplicate id"):
        fp.parse_signalp(path)


def test_predgpi_scores_parse(phaseb_fixtures):
    calls = fp.parse_predgpi_scores(phaseb_fixtures / "predgpi_scores.tsv")
    assert calls[A] == fp.GpiCall("weakly", 0.55, 80, 0.0095988, -0.708166)
    assert calls[B] == fp.GpiCall("none", 0.0, None, 0.878876, -1.68145)
    assert calls[C].omega == 1512 and calls[C].prob == 1.0
    assert calls[D] == fp.GpiCall("too_short", 0.0, None, None, None)


def test_predgpi_wrong_columns_stop(tmp_path):
    path = tmp_path / "g.tsv"
    path.write_text("id\tgpi\n" + f"{A}\tyes\n")
    with pytest.raises(fp.OutputFormatError, match="columns"):
        fp.parse_predgpi_scores(path)


def test_ser_thr_fraction():
    assert fp.ser_thr_fraction("STSTAAAA") == 0.5
    assert fp.ser_thr_fraction("") == 0.0


PREDGPI_HOME = os.environ.get("PREDGPI_HOME")


@pytest.mark.skipif(not PREDGPI_HOME, reason="needs `module load predgpi/202001`")
def test_predgpi_wrapper_matches_cli(tmp_path):
    from conftest import STEP1_DIR

    python = shutil.which("python")  # the module's Python 3.9 with numpy
    fasta = tmp_path / "in.fasta"
    text = open(os.path.join(PREDGPI_HOME, "testdata", "test.fasta")).read()
    fasta.write_text(text + ">short1\nMKVLAAGIVALLLAAG\n")
    cli_out, wrap_out = tmp_path / "cli.gff3", tmp_path / "scores.tsv"
    cli = os.path.join(PREDGPI_HOME, "predgpi.py")
    subprocess.run([python, cli, "-f", str(fasta), "-o", str(cli_out), "-m", "gff3"], check=True)
    wrapper = STEP1_DIR / "jobs" / "predgpi_scores.py"
    subprocess.run(
        [python, str(wrapper), "--fasta", str(fasta), "--out", str(wrap_out)], check=True
    )
    scores = fp.parse_predgpi_scores(wrap_out)
    n = 0
    for line in cli_out.read_text().splitlines():
        f = line.split("\t")
        call = scores[f[0]]
        if f[2] == "GPI-anchor":
            assert call.call in ("highly_probable", "probable", "weakly")
            assert call.omega == int(f[3]) and call.prob == float(f[5])
        else:
            assert call.call in ("none", "too_short")
        n += 1
    assert n == len(scores) >= 3
    assert scores["short1"].call == "too_short"
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_feature_parsers.py -q`
Expected: FAIL, collection error `ModuleNotFoundError: No module named 'feature_parsers'`.

- [ ] **Step 4: Write `feature_parsers.py`**

```python
"""Parse SignalP 6 and PredGPI outputs, and compute Ser+Thr fraction. Standard library only.

SignalP 6 (`--organism eukarya --mode fast`) writes `prediction_results.txt`: two `#` header
lines, the second `# ID<TAB>Prediction<TAB>OTHER<TAB>SP(Sec/SPI)<TAB>CS Position`, then one TAB
separated row per protein. The CS column is `CS pos: 24-25. Pr: 0.5487` for SP rows and empty
for OTHER rows. Columns are found by header name, never by whitespace splitting (a whitespace
split mis-columns the CS field; see analysis/cocci_repeats/01_signalp.sh). `output.gff3` has
one `signal_peptide` line per SP protein (start 1, end = last residue of the signal peptide).

PredGPI rows come from jobs/predgpi_scores.py (columns id, length, gpi_call, gpi_prob, omega,
fpr, svm). Both readers accept plain or gzip files (gzip found by its magic bytes).
"""

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from gaf import open_text

SP_COLUMNS = ("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "CS Position")
_CS = re.compile(r"^CS pos: (\d+)-(\d+)\. Pr: ([0-9.]+)$")
GPI_COLUMNS = ("id", "length", "gpi_call", "gpi_prob", "omega", "fpr", "svm")
GPI_CALLS = ("highly_probable", "probable", "weakly", "none", "too_short")


class OutputFormatError(ValueError):
    """A tool output does not have the expected format."""


@dataclass(frozen=True, slots=True)
class SignalPCall:
    prediction: str
    other_prob: float
    sp_prob: float
    cs_end: int | None
    cs_prob: float | None


def parse_signalp(path: str | Path) -> dict[str, SignalPCall]:
    calls: dict[str, SignalPCall] = {}
    header = None
    with open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if line.startswith("# ID"):
                header = line[2:].split("\t")
                missing = [c for c in SP_COLUMNS if c not in header]
                if missing:
                    raise OutputFormatError(f"{path}: header lacks {missing}")
                continue
            if line.startswith("#") or not line:
                continue
            if header is None:
                raise OutputFormatError(f"{path}:{lineno}: data row before the '# ID' header")
            fields = line.split("\t")
            if len(fields) != len(header):
                raise OutputFormatError(f"{path}:{lineno}: {len(fields)} fields, not {len(header)}")
            row = dict(zip(header, fields, strict=True))
            pred, cs = row["Prediction"], row["CS Position"].strip()
            if pred not in ("SP", "OTHER"):
                raise OutputFormatError(f"{path}:{lineno}: unexpected prediction {pred!r}")
            cs_end = cs_prob = None
            if pred == "SP":
                match = _CS.match(cs)
                if not match:
                    raise OutputFormatError(f"{path}:{lineno}: SP row with CS field {cs!r}")
                cs_end, cs_prob = int(match.group(1)), float(match.group(3))
            elif cs:
                raise OutputFormatError(f"{path}:{lineno}: OTHER row with a CS field")
            if row["ID"] in calls:
                raise OutputFormatError(f"{path}:{lineno}: duplicate id {row['ID']}")
            calls[row["ID"]] = SignalPCall(
                pred, float(row["OTHER"]), float(row["SP(Sec/SPI)"]), cs_end, cs_prob
            )
    if header is None:
        raise OutputFormatError(f"{path}: no '# ID' header line")
    return calls


def parse_signalp_gff(path: str | Path) -> dict[str, int]:
    """id -> end of the signal_peptide feature."""
    ends: dict[str, int] = {}
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "signal_peptide":
                raise OutputFormatError(f"{path}: unexpected GFF line {raw[:60]!r}")
            ends[f[0]] = int(f[4])
    return ends


def check_signalp_consistency(calls: dict[str, SignalPCall], gff_ends: dict[str, int]) -> None:
    sp = {k: c.cs_end for k, c in calls.items() if c.prediction == "SP"}
    if sp != gff_ends:
        diff = sorted(set(sp.items()) ^ set(gff_ends.items()))[:5]
        raise OutputFormatError(f"prediction_results.txt and output.gff3 disagree: {diff}")


@dataclass(frozen=True, slots=True)
class GpiCall:
    call: str
    prob: float
    omega: int | None
    fpr: float | None
    svm: float | None


def parse_predgpi_scores(path: str | Path) -> dict[str, GpiCall]:
    calls: dict[str, GpiCall] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GPI_COLUMNS:
            raise OutputFormatError(f"{path}: columns {reader.fieldnames} != {GPI_COLUMNS}")
        for r in reader:
            if r["gpi_call"] not in GPI_CALLS:
                raise OutputFormatError(f"{path}: unknown gpi_call {r['gpi_call']!r}")
            if r["id"] in calls:
                raise OutputFormatError(f"{path}: duplicate id {r['id']}")
            calls[r["id"]] = GpiCall(
                r["gpi_call"],
                float(r["gpi_prob"]),
                int(r["omega"]) if r["omega"] else None,
                float(r["fpr"]) if r["fpr"] else None,
                float(r["svm"]) if r["svm"] else None,
            )
    return calls


def ser_thr_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("S") + sequence.count("T")) / len(sequence)
```

- [ ] **Step 5: Write `jobs/predgpi_scores.py`**

```python
#!/usr/bin/env python
"""PredGPI (module predgpi/202001) with its continuous scores, one TSV row per protein.

The PredGPI CLI (`predgpi.py -m gff3`) writes only a class: a `GPI-anchor` line with score
1.0 (FPR <= 0.0015), 0.70 (<= 0.005) or 0.55 (<= 0.01), else a `Chain` line. This wrapper
calls the same functions of predgpi.py (predGpipe, the HMM and SVM files of the installation)
and also writes the estimated false positive rate and the SVM output, which the hybrid
candidate H needs as a GPI score. The class rules and residue substitutions copy
predgpi.py main(); test_predgpi_wrapper_matches_cli checks the classes against the CLI.

Runs under the module's Python 3.9 (`module load predgpi/202001`; PREDGPI_HOME must be set).
Keep this file Python 3.9 compatible. Columns: id, length, gpi_call, gpi_prob, omega, fpr, svm.
"""

import argparse
import gzip
import os
import sys

COLUMNS = ("id", "length", "gpi_call", "gpi_prob", "omega", "fpr", "svm")


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def classify(fpr):
    if fpr <= 0.0015:
        return "highly_probable", "1.0"
    if fpr <= 0.005:
        return "probable", "0.70"
    if fpr <= 0.01:
        return "weakly", "0.55"
    return "none", "0"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    home = os.environ.get("PREDGPI_HOME")
    if not home:
        print("STOP: PREDGPI_HOME is not set; run `module load predgpi/202001`", file=sys.stderr)
        return 2
    sys.path.insert(0, home)
    import predgpi
    from predgpilib.hmm import HMM_IO
    from predgpilib.svm import SVMLike

    hmm = HMM_IO.get_hmm(os.path.join(home, "GPIDAT", "PHMM.TOT.ss.mod"))
    svm = SVMLike.getSVMLight(os.path.join(home, "GPIDAT", "MOD"))
    tmp = args.out + ".tmp"
    seen = set()
    with open(tmp, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for name, seq in read_fasta(args.fasta):
            if name in seen:
                print("STOP: duplicate id " + name, file=sys.stderr)
                os.remove(tmp)
                return 2
            seen.add(name)
            if len(seq) <= 40:
                out.write("\t".join([name, str(len(seq)), "too_short", "0", "", "", ""]) + "\n")
                continue
            seq_t = seq.replace("U", "C").replace("Z", "A").replace("B", "A").replace("X", "A")
            _, cut, svmout, fpr = predgpi.predGpipe(seq_t, svm, hmm)
            call, prob = classify(fpr)
            omega = str(len(seq) - cut) if call != "none" else ""
            row = [name, str(len(seq)), call, prob, omega, "%.6g" % fpr, "%.6g" % svmout]
            out.write("\t".join(row) + "\n")
    os.replace(tmp, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare/test_feature_parsers.py -q`
Expected: `11 passed, 1 skipped` (the wrapper test needs the module).
Run: `module load predgpi/202001 && $PY -m pytest tests/step1_compare/test_feature_parsers.py -q`
Expected: `12 passed`. The wrapper test runs the PredGPI CLI and the wrapper on `$PREDGPI_HOME/testdata/test.fasta` plus a 16 aa sequence and compares class, omega site and score line by line.

- [ ] **Step 7: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/feature_parsers.py analysis/step1_compare/jobs/predgpi_scores.py \
  tests/step1_compare/fixtures/phaseb/signalp_prediction_results.txt \
  tests/step1_compare/fixtures/phaseb/signalp_output.gff3 \
  tests/step1_compare/fixtures/phaseb/predgpi_scores.tsv tests/step1_compare/test_feature_parsers.py
git commit -m "step1_compare: SignalP 6 and PredGPI parsers; PredGPI wrapper with FPR and SVM scores"
```

### Task 8: Feature table and coverage report (07)

**Files:**
- Create: `analysis/step1_compare/07_build_features.py`
- Modify: `analysis/step1_compare/runinfo.py` (add `"d8_run.json": "03_triage_pm.py"` to `PRODUCERS`; final state shown in Task 2 Step 6)
- Test: `tests/step1_compare/test_phaseb_features.py`

**Interfaces:**
- Consumes: `feature_parsers.*` (Task 7); `seqsets.UNIQUE_COLUMNS`, `seqsets.MEMBER_COLUMNS` (Task 2); Phase A `truth_set_triaged.tsv.gz` and `d8_run.json`; J1 outputs `phaseb/signalp/part_*/prediction_results.txt.gz`, `output.gff3.gz`, `phaseb/predgpi/part_*.tsv.gz`.
- Produces: `07_build_features.py` with `OUTPUT_NAMES`, `FEATURE_COLUMNS`, `UNIQUE_FEATURE_COLUMNS`, `TRUTH_COLUMNS`, `MEMBER_FEATURE_COLUMNS`, `COVERAGE_COLUMNS`, `FeatureError`, `merge_parts`, `feature_rows`, `member_rows`, `coverage_rows`, `read_tool_outputs`, `run(work, allow_missing, provenance=None)`, `main(argv=None)`. Outputs: `features_unique.tsv.gz`, `features.tsv.gz`, `feature_coverage.tsv`, `features_run.json`.

- [ ] **Step 1: Write the failing test**

```python
import gzip
import json

import pytest
import seqsets
import truth_table
from conftest import load_script

A, B, C, D = "a" * 64, "b" * 64, "c" * 64, "d" * 64
SEQS = {A: "MKSTST" * 20, B: "M" * 84, C: "MS" * 600, D: "MKT" * 10}


def _work(tmp_path, fixtures_dir, drop_gpi=None):
    work = tmp_path / "w"
    out = work / "phaseb"
    unique = []
    crow = 0
    for i, h in enumerate(sorted(SEQS)):
        seq = SEQS[h]
        long = len(seq) > 1022
        unique.append(
            {
                "row": str(i),
                "seq_sha256": h,
                "length": str(len(seq)),
                "cterm_row": str(crow) if long else "",
                "sequence": seq,
            }
        )
        crow += long
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    members = [
        {"set_id": "truth", "source_id": "Scer_SGD", "gene_id": "S1", "seq_sha256": A},
        {"set_id": "truth", "source_id": "Calb_CGD", "gene_id": "C1", "seq_sha256": A},
        {"set_id": "truth", "source_id": "Scer_SGD", "gene_id": "S2", "seq_sha256": C},
        {"set_id": "Cimm", "source_id": "Cimm", "gene_id": "CIMG_1", "seq_sha256": B},
        {"set_id": "Cimm", "source_id": "Cimm", "gene_id": "CIMG_2", "seq_sha256": D},
    ]
    for m in members:
        m["length"] = str(len(SEQS[m["seq_sha256"]]))
    truth_table.write_tsv(out / "sequence_members.tsv.gz", seqsets.MEMBER_COLUMNS, members)
    truth = [
        {"source_id": "Scer_SGD", "gene_id": "S1", "label": "P-ext", "subset": "wall",
         "stratum": "P-gpi", "d8_class": "P-gpi", "homology_only": "no", "role": "train"},
        {"source_id": "Calb_CGD", "gene_id": "C1", "label": "N-int", "subset": "",
         "stratum": "N-int", "d8_class": "", "homology_only": "no", "role": "train"},
        {"source_id": "Scer_SGD", "gene_id": "S2", "label": "N-sec", "subset": "",
         "stratum": "N-sec", "d8_class": "", "homology_only": "yes", "role": "train"},
    ]  # fmt: skip
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", list(truth[0]), truth)
    (work / "d8_run.json").write_text(json.dumps({"all_sources": True}))
    (out / "prepare_run.json").write_text(json.dumps({"all_sources": True}))
    fx = fixtures_dir / "phaseb"
    sp = out / "signalp" / "part_000"
    sp.mkdir(parents=True)
    for name, dest in (("signalp_prediction_results.txt", "prediction_results.txt.gz"),
                       ("signalp_output.gff3", "output.gff3.gz")):  # fmt: skip
        (sp / dest).write_bytes(gzip.compress((fx / name).read_bytes()))
    gpi_lines = (fx / "predgpi_scores.tsv").read_text().splitlines(keepends=True)
    if drop_gpi:
        gpi_lines = [x for x in gpi_lines if not x.startswith(drop_gpi)]
    (out / "predgpi").mkdir()
    (out / "predgpi" / "part_000.tsv.gz").write_bytes(gzip.compress("".join(gpi_lines).encode()))
    return work, out


def test_features_join_members_truth_and_embedding_rows(tmp_path, fixtures_dir):
    work, out = _work(tmp_path, fixtures_dir)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 0
    feats = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "features_unique.tsv.gz")}
    assert len(feats) == 4
    assert feats[C]["sp_prediction"] == "SP" and feats[C]["sp_cs_end"] == "24"
    assert feats[A]["gpi_call"] == "weakly" and feats[A]["gpi_omega"] == "80"
    assert feats[A]["ser_thr_frac"] == "0.666667"
    assert feats[C]["cterm_row"] == "0" and feats[A]["cterm_row"] == ""
    rows = truth_table.read_tsv(out / "features.tsv.gz")
    assert len(rows) == 5
    s1 = next(r for r in rows if r["gene_id"] == "S1")
    c1 = next(r for r in rows if r["gene_id"] == "C1")
    assert s1["label"] == "P-ext" and s1["d8_class"] == "P-gpi" and c1["label"] == "N-int"
    assert s1["emb_row"] == c1["emb_row"] == feats[A]["row"]  # one embedding per sequence
    cimm = next(r for r in rows if r["gene_id"] == "CIMG_1")
    assert cimm["label"] == "" and cimm["sp_prediction"] == "OTHER"
    cov = {r["set_id"]: r for r in truth_table.read_tsv(out / "feature_coverage.tsv")}
    assert cov["truth"] == {
        "set_id": "truth",
        "members": "3",
        "unique_sequences": "2",
        "no_signalp": "0",
        "no_predgpi": "0",
        "over_1022": "1",
        "gpi_too_short": "0",
    }
    assert cov["Cimm"]["gpi_too_short"] == "1"
    log = json.loads((out / "features_run.json").read_text())
    assert log["all_sources"] is True and log["missing_signalp"] == 0
    assert log["sp_predictions"] == {"OTHER": 2, "SP": 2}


def test_missing_predgpi_call_stops_unless_allowed(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir, drop_gpi=B)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "1 no PredGPI call" in capsys.readouterr().err
    assert not (out / "features.tsv.gz").exists()
    assert build.main(["--work-dir", str(work), "--allow-missing-calls"]) == 0
    cov = {r["set_id"]: r for r in truth_table.read_tsv(out / "feature_coverage.tsv")}
    assert cov["Cimm"]["no_predgpi"] == "1" and cov["truth"]["no_predgpi"] == "0"


def test_stale_tool_output_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique[:3])
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work), "--allow-missing-calls"]) == 2
    assert "not unique sequences" in capsys.readouterr().err


def test_id_in_two_parts_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    src = out / "predgpi" / "part_000.tsv.gz"
    (out / "predgpi" / "part_001.tsv.gz").write_bytes(src.read_bytes())
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "two parts" in capsys.readouterr().err


def test_truth_member_without_truth_row_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    rows = truth_table.read_tsv(work / "truth_set_triaged.tsv.gz")[:2]
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", list(rows[0]), rows)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "Scer_SGD:S2" in capsys.readouterr().err


def test_partial_triage_is_refused(tmp_path, fixtures_dir, capsys):
    work, _ = _work(tmp_path, fixtures_dir)
    (work / "d8_run.json").write_text(json.dumps({"all_sources": False}))
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "03_triage_pm.py" in capsys.readouterr().err


@pytest.mark.parametrize("name", ["signalp", "predgpi"])
def test_missing_tool_output_folder_stops(tmp_path, fixtures_dir, capsys, name):
    import shutil

    work, out = _work(tmp_path, fixtures_dir)
    shutil.rmtree(out / name)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "no J1 output" in capsys.readouterr().err
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_features.py -q`
Expected (checked on the prototype with the file removed): `8 failed`, each with `FileNotFoundError` for `07_build_features.py` (raised by `conftest.load_script`).

- [ ] **Step 3: Write `07_build_features.py`**

```python
#!/usr/bin/env python3
"""Phase B step 4: merge SignalP 6, PredGPI and Ser+Thr into one feature table (D2).

Reads from $STEP1_WORKDIR/phaseb/: unique_sequences.tsv.gz, sequence_members.tsv.gz,
signalp/part_*/prediction_results.txt.gz and output.gff3.gz, predgpi/part_*.tsv.gz (all from
J1). Reads $STEP1_WORKDIR/truth_set_triaged.tsv.gz (03) for the truth columns. Writes to
$STEP1_WORKDIR/phaseb/:

  features_unique.tsv.gz  one row per unique sequence (row, seq_sha256, features)
  features.tsv.gz         one row per member (set_id, source_id, gene_id), truth columns for
                          truth members, the features, and the embedding rows emb_row and
                          emb_cterm_row (rows of <model>.nterm.npy and <model>.cterm.npy)
  feature_coverage.tsv    per set_id: members, unique sequences, no SignalP call,
                          no PredGPI call, longer than 1,022 aa, PredGPI too_short
  features_run.json       input hashes, counts, all_sources, git commit, arguments

STOP (exit 2, no output): a tool output id that is not a unique seq_sha256 (stale output); an id
in two parts; SignalP prediction_results and output.gff3 that disagree; a truth member with no
row in truth_set_triaged.tsv.gz; d8_run.json without all_sources: true (unless
--allow-partial-truth-set); a unique sequence without a SignalP or PredGPI call (unless
--allow-missing-calls, which writes empty feature fields and counts them in the coverage).
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import feature_parsers as fp
import manifest
import paths
import runinfo
import truth_table

OUTPUT_NAMES = (
    "features_unique.tsv.gz",
    "features.tsv.gz",
    "feature_coverage.tsv",
    "features_run.json",
)
FEATURE_COLUMNS = (
    "ser_thr_frac",
    "sp_prediction",
    "sp_prob",
    "sp_other_prob",
    "sp_cs_end",
    "sp_cs_prob",
    "gpi_call",
    "gpi_prob",
    "gpi_omega",
    "gpi_fpr",
    "gpi_svm",
)
UNIQUE_FEATURE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", *FEATURE_COLUMNS)
TRUTH_COLUMNS = ("label", "subset", "stratum", "d8_class", "homology_only", "role")
MEMBER_FEATURE_COLUMNS = (
    "set_id",
    "source_id",
    "gene_id",
    "seq_sha256",
    "length",
    *TRUTH_COLUMNS,
    *FEATURE_COLUMNS,
    "emb_row",
    "emb_cterm_row",
)
COVERAGE_COLUMNS = (
    "set_id",
    "members",
    "unique_sequences",
    "no_signalp",
    "no_predgpi",
    "over_1022",
    "gpi_too_short",
)


class FeatureError(ValueError):
    """The tool outputs do not cover the unique sequences as expected."""


def _fmt(value) -> str:
    return "" if value is None else str(value)


def merge_parts(parts: list[dict], what: str) -> dict:
    merged: dict = {}
    for part in parts:
        dup = merged.keys() & part.keys()
        if dup:
            raise FeatureError(f"{what}: id {sorted(dup)[0]} occurs in two parts")
        merged.update(part)
    return merged


def feature_rows(unique, sp_calls, gpi_calls, allow_missing: bool):
    known = {r["seq_sha256"] for r in unique}
    for what, calls in (("SignalP", sp_calls), ("PredGPI", gpi_calls)):
        stale = calls.keys() - known
        if stale:
            raise FeatureError(
                f"{what} output has {len(stale)} ids that are not unique sequences "
                f"(for example {sorted(stale)[0]}); re-run J1 on the current unique_sequences"
            )
    missing_sp = [r["seq_sha256"] for r in unique if r["seq_sha256"] not in sp_calls]
    missing_gpi = [r["seq_sha256"] for r in unique if r["seq_sha256"] not in gpi_calls]
    if (missing_sp or missing_gpi) and not allow_missing:
        raise FeatureError(
            f"{len(missing_sp)} sequences have no SignalP call and {len(missing_gpi)} no "
            "PredGPI call; use --allow-missing-calls to write the table anyway"
        )
    rows = []
    for u in unique:
        sp = sp_calls.get(u["seq_sha256"])
        gpi = gpi_calls.get(u["seq_sha256"])
        rows.append(
            {
                "row": u["row"],
                "seq_sha256": u["seq_sha256"],
                "length": u["length"],
                "cterm_row": u["cterm_row"],
                "ser_thr_frac": f"{fp.ser_thr_fraction(u['sequence']):.6f}",
                "sp_prediction": sp.prediction if sp else "",
                "sp_prob": _fmt(sp.sp_prob if sp else None),
                "sp_other_prob": _fmt(sp.other_prob if sp else None),
                "sp_cs_end": _fmt(sp.cs_end if sp else None),
                "sp_cs_prob": _fmt(sp.cs_prob if sp else None),
                "gpi_call": gpi.call if gpi else "",
                "gpi_prob": _fmt(gpi.prob if gpi else None),
                "gpi_omega": _fmt(gpi.omega if gpi else None),
                "gpi_fpr": _fmt(gpi.fpr if gpi else None),
                "gpi_svm": _fmt(gpi.svm if gpi else None),
            }
        )
    return rows


def member_rows(members, features_by_hash, truth_by_key):
    out = []
    for m in members:
        f = features_by_hash[m["seq_sha256"]]
        truth = {c: "" for c in TRUTH_COLUMNS}
        if m["set_id"] == "truth":
            t = truth_by_key.get((m["source_id"], m["gene_id"]))
            if t is None:
                raise FeatureError(
                    f"truth member {m['source_id']}:{m['gene_id']} has no row in "
                    "truth_set_triaged.tsv.gz; re-run 03_triage_pm.py and 05"
                )
            truth = {c: t.get(c, "") for c in TRUTH_COLUMNS}
        out.append(
            {
                **{c: m[c] for c in ("set_id", "source_id", "gene_id", "seq_sha256", "length")},
                **truth,
                **{c: f[c] for c in FEATURE_COLUMNS},
                "emb_row": f["row"],
                "emb_cterm_row": f["cterm_row"],
            }
        )
    return out


def coverage_rows(members, features_by_hash):
    by_set: dict[str, list] = {}
    for m in members:
        by_set.setdefault(m["set_id"], []).append(m)
    rows = []
    for set_id, ms in by_set.items():
        hashes = {m["seq_sha256"] for m in ms}
        feats = [features_by_hash[h] for h in hashes]
        rows.append(
            {
                "set_id": set_id,
                "members": str(len(ms)),
                "unique_sequences": str(len(hashes)),
                "no_signalp": str(sum(f["sp_prediction"] == "" for f in feats)),
                "no_predgpi": str(sum(f["gpi_call"] == "" for f in feats)),
                "over_1022": str(sum(f["cterm_row"] != "" for f in feats)),
                "gpi_too_short": str(sum(f["gpi_call"] == "too_short" for f in feats)),
            }
        )
    return rows


def read_tool_outputs(out: Path) -> tuple[dict, dict, dict]:
    sp_dirs = sorted((out / "signalp").glob("part_*"))
    gpi_files = sorted((out / "predgpi").glob("part_*.tsv.gz"))
    if not sp_dirs or not gpi_files:
        raise FeatureError(f"no J1 output under {out}/signalp or {out}/predgpi")
    hashes, sp_parts = {}, []
    for d in sp_dirs:
        pred, gff = d / "prediction_results.txt.gz", d / "output.gff3.gz"
        calls = fp.parse_signalp(pred)
        fp.check_signalp_consistency(calls, fp.parse_signalp_gff(gff))
        sp_parts.append(calls)
        hashes[str(pred.relative_to(out))] = manifest.sha256_file(pred)
    gpi_parts = []
    for f in gpi_files:
        gpi_parts.append(fp.parse_predgpi_scores(f))
        hashes[str(f.relative_to(out))] = manifest.sha256_file(f)
    return merge_parts(sp_parts, "SignalP"), merge_parts(gpi_parts, "PredGPI"), hashes


def run(work: Path, allow_missing: bool, provenance: dict | None = None):
    out = Path(work) / "phaseb"
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    members = truth_table.read_tsv(out / "sequence_members.tsv.gz")
    truth_by_key = {
        (r["source_id"], r["gene_id"]): r
        for r in truth_table.read_tsv(Path(work) / "truth_set_triaged.tsv.gz")
    }
    sp_calls, gpi_calls, hashes = read_tool_outputs(out)
    feats = feature_rows(unique, sp_calls, gpi_calls, allow_missing)
    by_hash = {f["seq_sha256"]: f for f in feats}
    mrows = member_rows(members, by_hash, truth_by_key)
    cov = coverage_rows(members, by_hash)
    log = {
        "tool_outputs_sha256": hashes,
        "unique_sequences": len(feats),
        "members": len(mrows),
        "missing_signalp": sum(f["sp_prediction"] == "" for f in feats),
        "missing_predgpi": sum(f["gpi_call"] == "" for f in feats),
        "sp_predictions": dict(Counter(f["sp_prediction"] for f in feats)),
        "gpi_calls": dict(Counter(f["gpi_call"] for f in feats)),
        "unique_sequences_sha256": manifest.sha256_file(out / "unique_sequences.tsv.gz"),
        "all_sources": None,
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
        "arguments": [],
    }
    log.update(provenance or {})
    runinfo.atomic_write_all(
        out,
        {
            "features_unique.tsv.gz": lambda p: truth_table.write_tsv(
                p, UNIQUE_FEATURE_COLUMNS, feats
            ),
            "features.tsv.gz": lambda p: truth_table.write_tsv(p, MEMBER_FEATURE_COLUMNS, mrows),
            "feature_coverage.tsv": lambda p: truth_table.write_tsv(p, COVERAGE_COLUMNS, cov),
            "features_run.json": lambda p: runinfo.write_json(p, log),
        },
    )
    return feats, mrows, cov, log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--allow-missing-calls", action="store_true")
    parser.add_argument("--allow-partial-truth-set", action="store_true")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        runinfo.require_full(
            work / "d8_run.json", "truth_set_triaged.tsv.gz", args.allow_partial_truth_set
        )
        provenance = {
            "all_sources": runinfo.says_all_sources(work / "d8_run.json")
            and runinfo.says_all_sources(work / "phaseb" / "prepare_run.json"),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        _, _, cov, log = run(work, args.allow_missing_calls, provenance)
    except (FeatureError, fp.OutputFormatError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in cov:
        print("\t".join(f"{k}={c[k]}" for k in COVERAGE_COLUMNS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Add the `d8_run.json` producer**

In `analysis/step1_compare/runinfo.py`, add the line `"d8_run.json": "03_triage_pm.py",` to `PRODUCERS` after the `keyword_tier_run.json` line (Task 2 Step 6 shows the result).

- [ ] **Step 5: Run the tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_features.py tests/step1_compare/test_runinfo.py tests/step1_compare/test_paths.py -q`
Expected: all pass (`test_phaseb_features.py`: 8 passed).

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/07_build_features.py analysis/step1_compare/runinfo.py \
  tests/step1_compare/test_phaseb_features.py
git commit -m "step1_compare: 07 merges SignalP, PredGPI and Ser+Thr into features.tsv.gz"
```

### Task 9: SLURM jobs J0, J1, J2

**Files:**
- Create: `analysis/step1_compare/jobs/j0_pilot.sh`, `analysis/step1_compare/jobs/j1_features.sh`, `analysis/step1_compare/jobs/j2_embed.sh`, `tests/step1_compare/fixtures/phaseb/stub_signalp/bin/signalp6` (executable), `tests/step1_compare/fixtures/phaseb/stub_signalp/modules/signalp/6-gpu`
- Test: `tests/step1_compare/test_phaseb_jobs.py`

**Interfaces:**
- Consumes: `jobs/throughput_pilot.py`, `jobs/gpu_cpu_diff.py` (Task 6); `jobs/predgpi_scores.py` (Task 7); `jobs/embed_chunks.py` (Task 4); `05_prepare_sequences.py` and `07_build_features.py` (tests, Tasks 2 and 8).
- Produces: environment contract of the jobs: `PROJ_ROOT`, `STEP1_WORKDIR` (required), `STEP1_ENV_PY` (optional), `J1_PARTS` (default 8), `J1_SIGNALP_MODULE` (default `signalp/6-gpu`), `J1_PREDGPI_MODULE` (default `predgpi/202001`), `J2_JOB_COUNT` (default 1). J1 writes `phaseb/signalp/part_NNN/{prediction_results.txt.gz,output.gff3.gz,region_output.gff3.gz,input.sha256}` and `phaseb/predgpi/part_NNN.tsv.gz`, `part_NNN.input.sha256`.

Notes found while prototyping (keep them):
- `module load signalp/6-gpu` runs `conda activate`; it works in a `#!/bin/bash -l` job. A login shell also resets `SCRATCH` (to the job's `/scratch/$USER/$SLURM_JOB_ID`) and `MODULEPATH`, which is why the tests call `bash j1_features.sh` (the `module` function is exported) and pass the stub module by its file path.
- Each subshell sets `set -euo pipefail` itself. A first prototype ran the SignalP loop after a top-level `set +e`; a failing `signalp6` then produced empty `.gz` files that were copied as "done". `test_j1_signalp_failure_leaves_no_done_marker` pins the fix.
- `prediction_results.txt.gz` is copied last; it is the done marker of a part.

- [ ] **Step 1: Write the SignalP stub (test only)**

```bash
mkdir -p tests/step1_compare/fixtures/phaseb/stub_signalp/bin \
  tests/step1_compare/fixtures/phaseb/stub_signalp/modules/signalp
```

`tests/step1_compare/fixtures/phaseb/stub_signalp/bin/signalp6`:

```bash
#!/bin/bash
# Test stub of SignalP 6: writes the real output layout and calls every protein OTHER,
# except ids that start with "c", which get a 24-residue signal peptide.
while [ $# -gt 0 ]; do
  case "$1" in
    --fastafile) fasta=$2; shift 2 ;;
    --output_dir) out=$2; shift 2 ;;
    *) shift ;;
  esac
done
if [ "${STUB_SIGNALP_FAIL:-0}" = 1 ]; then echo "stub signalp6 failure" >&2; exit 1; fi
mkdir -p "$out"
{
  printf '# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 0\n'
  printf '# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n'
  awk '/^>/ { id = substr($1, 2)
              if (id ~ /^c/) printf "%s\tSP\t0.001000\t0.999000\tCS pos: 24-25. Pr: 0.9000\n", id
              else printf "%s\tOTHER\t1.000000\t0.000000\t\n", id }' "$fasta"
} > "$out/prediction_results.txt"
{
  echo '## gff-version 3'
  awk '/^>/ { id = substr($1, 2)
              if (id ~ /^c/) printf "%s\tSignalP-6.0\tsignal_peptide\t1\t24\t0.999\t.\t.\t.\n", id }' "$fasta"
} > "$out/output.gff3"
echo '## gff-version 3' > "$out/region_output.gff3"
echo "stub signalp6 done"
```

`tests/step1_compare/fixtures/phaseb/stub_signalp/modules/signalp/6-gpu`:

```tcl
#%Module1.0
## Test stub: puts a fake signalp6 on PATH (tests/step1_compare/test_phaseb_jobs.py).
set here [file dirname [file dirname [file dirname $ModulesCurrentModulefile]]]
prepend-path PATH $here/bin
```

Run: `chmod +x tests/step1_compare/fixtures/phaseb/stub_signalp/bin/signalp6`

- [ ] **Step 2: Write the failing test**

```python
"""Static checks of the Phase B SLURM scripts and a J1 smoke run with a stub SignalP."""

import gzip
import json
import os
import subprocess

import paths
import pytest
import truth_table
from conftest import load_script

JOBS = paths.STEP1_DIR / "jobs"
SCRIPTS = sorted(JOBS.glob("*.sh"))
STUB = paths.STEP1_DIR.parents[1] / "tests/step1_compare/fixtures/phaseb/stub_signalp"


def test_three_job_scripts_exist():
    assert [s.name for s in SCRIPTS] == ["j0_pilot.sh", "j1_features.sh", "j2_embed.sh"]


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_job_script_conventions(script):
    text = script.read_text()
    assert "BASH_SOURCE" not in text
    assert '"${SCRATCH:?' in text
    assert ': "${PROJ_ROOT:?' in text and ': "${STEP1_WORKDIR:?' in text
    assert "#SBATCH -p exfab" in text and "#SBATCH --gres=gpu:1" in text
    assert "#SBATCH --time=" in text
    assert text.startswith("#!/bin/bash -l\n")
    assert "set -euo pipefail" in text
    subprocess.run(["bash", "-n", str(script)], check=True)


def test_j1_stub_signalp_is_executable():
    assert os.access(STUB / "bin" / "signalp6", os.X_OK)


needs_modules = pytest.mark.skipif(
    "BASH_FUNC_module%%" not in os.environ
    or not os.path.exists("/opt/linux/rocky/8.x/x86_64/modules/predgpi/202001"),
    reason="needs the HPCC module system and predgpi/202001",
)


def _tiny_work(tmp_path):
    work = tmp_path / "work"
    (work / "downloads").mkdir(parents=True)
    seqs = {
        "G1": "MKLLSVLALLLAAGSAQAS" + "ST" * 40 + "GAAAGLLSLLAALLF",
        "G2": "MSEEKKQ" * 12,
        "G3": "MKV" * 15,
    }
    (work / "downloads" / "p.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    sets = tmp_path / "sets.tsv"
    sets.write_text("set_id\tkind\tlocation\tnote\nP\tdownload\tp.fasta\t\n")
    prepare = load_script("05_prepare_sequences")
    argv = ["--sets", str(sets), "--work-dir", str(work), "--input-dir", str(work / "downloads")]
    assert prepare.main(argv) == 0
    return work


def _run_j1(work, tmp_path, **extra):
    env = {
        **os.environ,
        "SCRATCH": str(tmp_path / "scratch"),
        "PROJ_ROOT": str(paths.STEP1_DIR.parents[1]),
        "STEP1_WORKDIR": str(work),
        "J1_PARTS": "2",
        "J1_SIGNALP_MODULE": str(STUB / "modules" / "signalp" / "6-gpu"),
        **extra,
    }
    (tmp_path / "scratch").mkdir(exist_ok=True)
    script = JOBS / "j1_features.sh"
    return subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)


@needs_modules
def test_j1_runs_resumes_and_feeds_07(tmp_path):
    work = _tiny_work(tmp_path)
    first = _run_j1(work, tmp_path)
    assert first.returncode == 0, first.stdout[-2000:] + first.stderr[-2000:]
    out = work / "phaseb"
    for part in ("part_000", "part_001"):
        assert (out / "signalp" / part / "prediction_results.txt.gz").exists()
        assert (out / "predgpi" / f"{part}.tsv.gz").exists()
    again = _run_j1(work, tmp_path)
    assert again.returncode == 0 and again.stdout.count("skip") == 4
    header = "source_id\tgene_id\tlabel\tsubset\tstratum\td8_class\thomology_only\trole\n"
    (work / "truth_set_triaged.tsv.gz").write_bytes(gzip.compress(header.encode()))
    (work / "d8_run.json").write_text(json.dumps({"all_sources": True}))
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 0
    cov = truth_table.read_tsv(out / "feature_coverage.tsv")
    assert cov == [
        {
            "set_id": "P",
            "members": "3",
            "unique_sequences": "3",
            "no_signalp": "0",
            "no_predgpi": "0",
            "over_1022": "0",
            "gpi_too_short": "0",
        }
    ]


@needs_modules
def test_j1_signalp_failure_leaves_no_done_marker(tmp_path):
    work = _tiny_work(tmp_path)
    run = _run_j1(work, tmp_path, STUB_SIGNALP_FAIL="1")
    assert run.returncode != 0
    assert not list((work / "phaseb" / "signalp").glob("part_*/prediction_results.txt.gz"))
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_jobs.py -q`
Expected (checked on the prototype with the scripts removed): `2 failed, 2 passed, 1 skipped`. `test_three_job_scripts_exist` fails (`[] == ['j0_pilot.sh', ...]`); `test_j1_runs_resumes_and_feeds_07` fails (bash exit 127, no script); `test_job_script_conventions` is skipped (empty parameter set); the stub test and the failure test pass trivially before the scripts exist.

- [ ] **Step 4: Write `jobs/j0_pilot.sh`**

```bash
#!/bin/bash -l
# J0 (spec 8): ESM-2 8M and 35M throughput on 2,000 truth sequences, then the spec 7
# GPU-CPU difference harness (which also repeats the GPU run to check repeatability).
# Submit (see the plan's run section):
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=... -o <log> -e <log> j0_pilot.sh
# PROJ_ROOT comes from the environment, never from the script location (that is wrong in
# SLURM spool directories). Work happens in node-local $SCRATCH; results are copied to
# $STEP1_WORKDIR/phaseb/j0/ before the job ends.
# exfab: one node (gpu12, 2x ada6000), usually less contended than short_gpu. 1 h is far above
# the expected run time (assumption: under 15 min, spec 8); J0 measures it.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH -J step1_j0
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_j0"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/jobs"
OUT="$STEP1_WORKDIR/phaseb/j0"
mkdir -p "$TMP" "$OUT"

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv > "$TMP/nvidia_smi.csv"
"$ENV_PY" "$S1/jobs/throughput_pilot.py" --work-dir "$STEP1_WORKDIR" \
  --out "$TMP/throughput.json" --device cuda --n 2000 --batch-sizes 8,16,32,64
"$ENV_PY" "$S1/jobs/gpu_cpu_diff.py" --work-dir "$STEP1_WORKDIR" \
  --out "$TMP/gpu_cpu_diff.json" --device-a cuda --device-b cpu --n 200 --n-long 20

cp "$TMP/throughput.json" "$TMP/gpu_cpu_diff.json" "$TMP/nvidia_smi.csv" "$OUT/"
echo "J0 done: $OUT"
```

- [ ] **Step 5: Write `jobs/j1_features.sh`**

```bash
#!/bin/bash -l
# J1 (spec 8, D2): SignalP 6 (GPU build, fast mode) and PredGPI on every unique sequence.
# Submit (see the plan's run section):
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=... -o <log> -e <log> j1_features.sh
# Input: $STEP1_WORKDIR/phaseb/unique_sequences.fasta.gz (05). The FASTA is split into
# J1_PARTS parts (default 8; record k of N goes to part floor(k*parts/N)). SignalP runs the
# parts one after another on the GPU; PredGPI (CPU only, single-threaded) runs the parts in
# parallel at the same time. Outputs go to $STEP1_WORKDIR/phaseb/signalp/part_NNN/ and
# $STEP1_WORKDIR/phaseb/predgpi/part_NNN.tsv.gz. A part whose output exists and whose
# recorded input SHA-256 equals the current part is skipped, so a rerun resumes.
# One job, not one job per part: measured rates (SignalP 6 fast on gpu12 about 220 proteins/s,
# job 29280458; PredGPI about 31 proteins/s per core on c01) give minutes per part, and the
# global job-size rule says not to split work into jobs of minutes.
# The CPU build of SignalP 6 has no model weights on this cluster (verified 2026-09-30); the
# GPU build refuses to run without a CUDA device (verified 2026-10-01). So the job needs a GPU.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH -J step1_j1
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_j1"
S1="$PROJ_ROOT/analysis/step1_compare"
OUT="$STEP1_WORKDIR/phaseb"
PARTS="${J1_PARTS:-8}"
SIGNALP_MODULE="${J1_SIGNALP_MODULE:-signalp/6-gpu}"  # tests pass a stub modulefile path
PREDGPI_MODULE="${J1_PREDGPI_MODULE:-predgpi/202001}"
mkdir -p "$TMP/in" "$TMP/sp" "$TMP/gpi" "$OUT/signalp" "$OUT/predgpi"

zcat "$OUT/unique_sequences.fasta.gz" > "$TMP/all.fasta"
N=$(grep -c '^>' "$TMP/all.fasta")
awk -v n="$N" -v p="$PARTS" -v dir="$TMP/in" '
  /^>/ { k++; part = int((k - 1) * p / n) }
  { printf "%s\n", $0 > sprintf("%s/part_%03d.fasta", dir, part) }' "$TMP/all.fasta"
echo "J1: $N sequences in $PARTS parts"

part_sha() { sha256sum "$TMP/in/$1.fasta" | cut -d' ' -f1; }

# copy_atomic SRC DEST: copy to a temp name next to DEST, then rename.
copy_atomic() { cp "$1" "$(dirname "$2")/.tmp.$(basename "$2")" && \
  mv "$(dirname "$2")/.tmp.$(basename "$2")" "$2"; }

run_predgpi() {  # $1 = part name; runs in a subshell with only the PredGPI module
  local part="$1" sha dest="$OUT/predgpi/$1.tsv.gz"
  sha=$(part_sha "$part")
  if [ -s "$dest" ] && [ "$(cat "$OUT/predgpi/$part.input.sha256" 2>/dev/null)" = "$sha" ]; then
    echo "  PredGPI skip $part (done)"; return 0
  fi
  python "$S1/jobs/predgpi_scores.py" --fasta "$TMP/in/$part.fasta" --out "$TMP/gpi/$part.tsv"
  gzip -n -c "$TMP/gpi/$part.tsv" > "$TMP/gpi/$part.tsv.gz"
  echo "$sha" > "$TMP/gpi/$part.input.sha256"
  copy_atomic "$TMP/gpi/$part.input.sha256" "$OUT/predgpi/$part.input.sha256"
  copy_atomic "$TMP/gpi/$part.tsv.gz" "$dest"
  echo "  PredGPI done $part"
}
export -f run_predgpi part_sha copy_atomic
export TMP OUT S1

set +e
(
  set -euo pipefail
  module load "$PREDGPI_MODULE"
  ls "$TMP/in" | sed 's/\.fasta$//' | xargs -P "$PARTS" -I{} bash -c 'set -euo pipefail; run_predgpi "$1"' _ {}
) > "$TMP/predgpi.log" 2>&1 &
GPI_PID=$!

(
  set -euo pipefail
  module load "$SIGNALP_MODULE"
  command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH" >&2; exit 1; }
  for f in "$TMP"/in/part_*.fasta; do
    part=$(basename "$f" .fasta)
    sha=$(part_sha "$part")
    dest="$OUT/signalp/$part"
    if [ -s "$dest/prediction_results.txt.gz" ] && \
       [ "$(cat "$dest/input.sha256" 2>/dev/null)" = "$sha" ]; then
      echo "  SignalP skip $part (done)"; continue
    fi
    rm -rf "$TMP/sp/$part"
    signalp6 --fastafile "$f" --organism eukarya --output_dir "$TMP/sp/$part" \
      --format none --mode fast --write_procs 8 --torch_num_threads 8 2>&1 | tail -2
    mkdir -p "$dest"
    for name in prediction_results.txt output.gff3 region_output.gff3; do
      gzip -n -c "$TMP/sp/$part/$name" > "$TMP/sp/$part/$name.gz"
    done
    echo "$sha" > "$TMP/sp/$part/input.sha256"
    copy_atomic "$TMP/sp/$part/output.gff3.gz" "$dest/output.gff3.gz"
    copy_atomic "$TMP/sp/$part/region_output.gff3.gz" "$dest/region_output.gff3.gz"
    copy_atomic "$TMP/sp/$part/input.sha256" "$dest/input.sha256"
    # prediction_results.txt.gz last: it is the done marker of the part
    copy_atomic "$TMP/sp/$part/prediction_results.txt.gz" "$dest/prediction_results.txt.gz"
    echo "  SignalP done $part"
  done
)
SP_STATUS=$?

wait "$GPI_PID"
GPI_STATUS=$?
set -e
cat "$TMP/predgpi.log"
echo "J1 done: SignalP status $SP_STATUS, PredGPI status $GPI_STATUS"
exit $(( SP_STATUS != 0 || GPI_STATUS != 0 ))
```

- [ ] **Step 6: Write `jobs/j2_embed.sh`**

```bash
#!/bin/bash -l
# J2 (spec 8, D3): ESM-2 8M and 35M, layer 6, residue-mean pooling, N- and C-terminal windows.
# Submit (see the plan's run section) with the job count and time from job_plan.json:
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=...,J2_JOB_COUNT=<n_jobs> \
#          --array=0-<n_jobs-1> --time=<time_minutes> -o <log> -e <log> j2_embed.sh
# Array task i embeds chunks k with k mod J2_JOB_COUNT == i (06_plan_embedding.py). Each
# finished chunk is written to node-local $SCRATCH and copied at once to
# $STEP1_WORKDIR/phaseb/emb/<model>/, so a killed task loses at most the chunk in progress
# and a resubmission skips the finished chunks after a hash check.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH -J step1_j2
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
JOB_INDEX="${SLURM_ARRAY_TASK_ID:-0}"
JOB_COUNT="${J2_JOB_COUNT:-1}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_j2_${JOB_INDEX}"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/jobs"
LOGDIR="$STEP1_WORKDIR/phaseb/emb/logs"
mkdir -p "$TMP" "$LOGDIR"

nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used --format=csv,noheader -l 30 \
  > "$TMP/gpu_util.csv" &
SMI=$!
set +e
"$ENV_PY" "$S1/jobs/embed_chunks.py" --work-dir "$STEP1_WORKDIR" \
  --job-index "$JOB_INDEX" --job-count "$JOB_COUNT" --device cuda --scratch-dir "$TMP/emb"
STATUS=$?
set -e
kill "$SMI" || true
gzip -n -c "$TMP/gpu_util.csv" > "$LOGDIR/gpu_util.${SLURM_JOB_ID:-local}_${JOB_INDEX}.csv.gz"
echo "J2 task $JOB_INDEX of $JOB_COUNT: exit $STATUS"
exit "$STATUS"
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_jobs.py tests/step1_compare/test_paths.py -q`
Expected: all pass (`test_phaseb_jobs.py`: 7 passed on the HPCC; the two J1 tests skip where the module system is missing). The J1 smoke test runs the real PredGPI on three sequences, the stub SignalP, a rerun (4 skips), and 07 on the result.

- [ ] **Step 8: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/jobs/j0_pilot.sh analysis/step1_compare/jobs/j1_features.sh \
  analysis/step1_compare/jobs/j2_embed.sh tests/step1_compare/fixtures/phaseb/stub_signalp \
  tests/step1_compare/test_phaseb_jobs.py
git commit -m "step1_compare: SLURM jobs J0 (pilot), J1 (SignalP 6, PredGPI), J2 (embeddings)"
```

### Task 10: Documentation of the Phase B outputs

**Files:**
- Modify: `analysis/step1_compare/COLUMNS.md` (append), `analysis/step1_compare/README.md`
- Test: `tests/step1_compare/test_phaseb_docs.py`

**Interfaces:**
- Consumes: `OUTPUT_NAMES` of 05, 06, 07 (Tasks 2, 3, 8).
- Produces: a `## <name>` heading in `COLUMNS.md` for every Phase B output.

- [ ] **Step 1: Write the failing test**

```python
"""COLUMNS.md documents every Phase B output."""

import paths
from conftest import load_script


def test_phaseb_outputs_have_columns_md_headings():
    names = []
    for script in ("05_prepare_sequences", "06_plan_embedding", "07_build_features"):
        names += load_script(script).OUTPUT_NAMES
    names += [
        "j0/throughput.json",
        "j0/gpu_cpu_diff.json",
        "signalp/part_NNN/",
        "predgpi/part_NNN.tsv.gz",
        "emb/<model>/<chunk_id>.npy",
        "emb/<model>.nterm.npy",
        "emb/chunk_manifest.tsv",
        "emb/embedding_run.json",
    ]
    headings = [
        line[3:].split()[0]
        for line in (paths.STEP1_DIR / "COLUMNS.md").read_text().splitlines()
        if line.startswith("## ")
    ]
    missing = [n for n in names if n not in headings]
    assert not missing, f"COLUMNS.md has no heading for {missing}"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_phaseb_docs.py -q`
Expected: FAIL with `COLUMNS.md has no heading for ['sequence_members.tsv.gz', ...]`.

- [ ] **Step 3: Append the Phase B section to `COLUMNS.md`**

Append this text at the end of `analysis/step1_compare/COLUMNS.md`:

```markdown

# Phase B outputs (features and embeddings)

Scripts 05, 06 and 07 and the jobs in `jobs/` write to `$STEP1_WORKDIR/phaseb/`. Later phases
join on `(source_id, gene_id)` and on `seq_sha256`. `emb_row` and `emb_cterm_row` give the rows
of the embedding matrices.

## sequence_members.tsv.gz (05_prepare_sequences.py, one row per protein of an input set)

| Column | Meaning |
|---|---|
| set_id | Input set in `sequence_sets.tsv` (`truth`, `uniprot_kw`, `<species>_proteome`, ...). |
| source_id | For `truth`: the `source_id` of `truth_sequences.tsv.gz`. Else equal to `set_id`. |
| gene_id | For `truth`: the GAF gene ID. For `uniprot_kw`: the UniProt accession. Else the first FASTA header token. |
| seq_sha256 | SHA-256 of the cleaned sequence (`seqhash.seq_sha256`). |
| length | Length of the cleaned sequence. |

## unique_sequences.tsv.gz (05_prepare_sequences.py, one row per unique sequence)

| Column | Meaning |
|---|---|
| row | Position in seq_sha256 order (0-based). Row of `emb/<model>.nterm.npy`. |
| seq_sha256 | SHA-256 of the cleaned sequence. Sort key. |
| length | Length of the cleaned sequence. |
| cterm_row | For sequences longer than 1,022 aa: row of `emb/<model>.cterm.npy` (0-based, in seq_sha256 order). Empty otherwise. |
| sequence | The cleaned sequence. Every character is in the ESM-2 alphabet. |

## unique_sequences.fasta.gz (05_prepare_sequences.py)

The unique sequences as FASTA, in `row` order. The header is the `seq_sha256`. J1 reads it.

## prepare_run.json (05_prepare_sequences.py)

Keys: `inputs` (per set: `file`, `sha256`, `members`, `empty_records`), `members`,
`unique_sequences`, `unique_residues`, `over_max_residues`, `all_sources` (`true` when
`sequence_run.json` and `keyword_tier_run.json` both say `all_sources: true`), `git_commit`,
`python`, `arguments`. No time stamp.

## chunk_plan.tsv (06_plan_embedding.py, one row per embedding chunk)

| Column | Meaning |
|---|---|
| chunk_id | `<window>_<NNNN>`, numbered in plan order. |
| window | `nterm` (first 1,022 residues) or `cterm` (last 1,022 residues; sequences > 1,022 aa only). |
| n_seqs | Number of chunk members. |
| residues | Sum of the window lengths. |
| members_sha256 | SHA-256 of the member seq_sha256 values joined by newlines, in chunk order. |

## chunk_members.tsv.gz (06_plan_embedding.py, one row per chunk member)

| Column | Meaning |
|---|---|
| chunk_id | As in `chunk_plan.tsv`. |
| position | Row in the chunk array (0-based). |
| row | `row` of the sequence in `unique_sequences.tsv.gz`. |
| seq_sha256 | SHA-256 of the sequence. |

## job_plan.json (06_plan_embedding.py)

Keys: `models`, `rates_residues_per_s`, `model_load_s`, `batch_size` (per model),
`rate_source` (`J0` or `assumed`), `throughput_sha256`, `chunk_residues`, `chunks`,
`residues_per_model`, `unique_sequences`, `unique_sequences_sha256`, `total_seconds`,
`n_jobs`, `seconds_per_job`, `time_minutes`, `git_commit`, `python`, `arguments`. The formula is
in the docstring of `chunk_plan.py`.

## j0/throughput.json (jobs/throughput_pilot.py, J0)

Keys: `schema` (`step1-phaseb-j0-throughput/1`), `device`, `gpu_name`, `torch`, `torch_cuda`,
`esm`, `n_proteins`, `residues`, `seed`, `sample_sha256`, `truth_sequences_sha256`,
`batch_sizes`, `model_load_s` (per model), `runs` (one object per model and batch size:
`model`, `batch_size`, `status` (`ok`, `oom`, or `skipped_<n>`), `seconds`, `proteins_per_s`,
`residues_per_s`, `peak_mem_bytes`, `dim`), `best` (per model: the `ok` run with the highest
`residues_per_s`), `git_commit`, `python`.

## j0/gpu_cpu_diff.json (jobs/gpu_cpu_diff.py, J0)

Keys: `schema` (`step1-phaseb-gpu-cpu-diff/1`), `device_a`, `device_b`, `gpu_name`, `torch`,
`n_windows`, `windows_sha256`, `batch_size`, `models` (per model: `max_abs_diff`,
`mean_abs_diff`, `max_rel_diff`, `min_cosine`, `repeat_identical_on_a`,
`repeat_max_abs_diff_on_a`), `git_commit`. No threshold is applied.

## signalp/part_NNN/ (jobs/j1_features.sh, J1)

SignalP 6 output for one part of `unique_sequences.fasta.gz`: `prediction_results.txt.gz`
(the done marker of the part), `output.gff3.gz`, `region_output.gff3.gz` (gzip `-n`), and
`input.sha256` (SHA-256 of the part FASTA). The format is described in `feature_parsers.py`.

## predgpi/part_NNN.tsv.gz (jobs/j1_features.sh, J1)

`jobs/predgpi_scores.py` output for one part: columns `id` (seq_sha256), `length`, `gpi_call`
(`highly_probable`, `probable`, `weakly`, `none`, `too_short` for 40 aa or less), `gpi_prob`
(the PredGPI CLI score: 1.0, 0.70, 0.55, or 0), `omega` (omega site, GPI calls only), `fpr`
(PredGPI estimated false positive rate; lower is more GPI-like), `svm` (SVM output).
`part_NNN.input.sha256` holds the SHA-256 of the part FASTA.

## features_unique.tsv.gz (07_build_features.py, one row per unique sequence)

| Column | Meaning |
|---|---|
| row, seq_sha256, length, cterm_row | As in `unique_sequences.tsv.gz`. |
| ser_thr_frac | (count of S + count of T) / length, 6 decimals. |
| sp_prediction | SignalP 6: `SP` or `OTHER`. Empty if SignalP gave no call. |
| sp_prob | SignalP 6 `SP(Sec/SPI)` probability. |
| sp_other_prob | SignalP 6 `OTHER` probability. |
| sp_cs_end | Last residue of the signal peptide (`CS pos: X-Y`, X). SP only. |
| sp_cs_prob | Probability of the cleavage site. SP only. |
| gpi_call, gpi_prob, gpi_omega, gpi_fpr, gpi_svm | PredGPI columns `gpi_call`, `gpi_prob`, `omega`, `fpr`, `svm`. |

## features.tsv.gz (07_build_features.py, one row per member)

| Column | Meaning |
|---|---|
| set_id, source_id, gene_id, seq_sha256, length | As in `sequence_members.tsv.gz`. |
| label, subset, stratum, d8_class, homology_only, role | From `truth_set_triaged.tsv.gz` for `truth` members. Empty for other sets. |
| (feature columns) | As in `features_unique.tsv.gz`. |
| emb_row | Row of `emb/<model>.nterm.npy`. |
| emb_cterm_row | Row of `emb/<model>.cterm.npy`. Empty for sequences of 1,022 aa or less: use `emb_row` (their C-terminal window is the whole sequence). |

## feature_coverage.tsv (07_build_features.py, one row per set_id)

| Column | Meaning |
|---|---|
| set_id | Input set. |
| members | Members of the set. |
| unique_sequences | Distinct seq_sha256 values in the set. |
| no_signalp | Unique sequences of the set without a SignalP call. |
| no_predgpi | Unique sequences of the set without a PredGPI call. |
| over_1022 | Unique sequences of the set longer than 1,022 aa. |
| gpi_too_short | Unique sequences of the set of 40 aa or less (PredGPI gives no score). |

## features_run.json (07_build_features.py)

Keys: `tool_outputs_sha256` (per part file), `unique_sequences`, `members`, `missing_signalp`,
`missing_predgpi`, `sp_predictions` and `gpi_calls` (counts), `unique_sequences_sha256`,
`all_sources` (`true` when `d8_run.json` and `prepare_run.json` say `all_sources: true`),
`git_commit`, `python`, `arguments`.

## emb/<model>/<chunk_id>.npy (jobs/embed_chunks.py, J2)

float32 array, one row per chunk member in `chunk_members.tsv.gz` order, ESM-2 layer 6,
mean over residue tokens (`surface_glyco.embeddings.get_esm_embeddings`). The sidecar
`<chunk_id>.json` is the done marker: `chunk_id`, `model`, `window`, `members_sha256`,
`repr_layer`, `batch_size`, `device`, `seconds`, `shape`, `dtype`, `array_sha256` (SHA-256 of
the raw array bytes).

## emb/<model>.nterm.npy (jobs/assemble_embeddings.py)

float32, shape (unique sequences, dim): row `row` is the N-terminal window embedding of that
sequence. `emb/<model>.cterm.npy` has shape (sequences > 1,022 aa, dim) and row `cterm_row`.
`assemble_embeddings.window_matrix` builds the full C-terminal matrix for M8-C and M35-C.
dim is 320 for `esm2_t6_8M_UR50D` and 480 for `esm2_t12_35M_UR50D`.

## emb/chunk_manifest.tsv (jobs/assemble_embeddings.py, one row per model and chunk)

| Column | Meaning |
|---|---|
| model | ESM-2 model name. |
| chunk_id, window, n_seqs, members_sha256 | As in `chunk_plan.tsv`. |
| array_sha256 | SHA-256 of the raw bytes of `emb/<model>/<chunk_id>.npy` (as in its JSON). |

## emb/embedding_run.json (jobs/assemble_embeddings.py)

Keys: `models` (per model and window: `shape`, `dtype`, `array_sha256`), `unique_sequences`,
`chunks`, `unique_sequences_sha256`, `chunk_plan_sha256`, `git_commit`, `python`.
```

- [ ] **Step 4: Update `README.md`**

```diff
--- a/analysis/step1_compare/README.md
+++ b/analysis/step1_compare/README.md
@@ -9,12 +9,14 @@
 
 - Use Python 3.12 or later: `PY=/usr/bin/python3.12`.
 - The code uses `dataclass(slots=True)`. The default `python3` on the HPCC (3.9) fails.
-- Every module uses the Python standard library only. A test checks this with the `ast` module
-  (`test_every_module_imports_only_stdlib_or_local` in `tests/step1_compare/test_paths.py`).
-- No script uses `BASH_SOURCE`. The folder has no shell script. A test checks this
-  (`test_no_shell_script_uses_bash_source`).
-- No SLURM job is needed. Run every script on a login or interactive node.
-  If you use a job later, write to `${SCRATCH:?}` and copy the results to /bigdata.
+- Every module in this folder uses the Python standard library only. A test checks this with
+  the `ast` module (`test_every_module_imports_only_stdlib_or_local` in
+  `tests/step1_compare/test_paths.py`). The subfolder `jobs/` is outside this rule: its scripts
+  need numpy, torch and fair-esm (conda env `adhesionPred`) or the PredGPI module Python.
+- No script uses `BASH_SOURCE`. A test checks this (`test_no_shell_script_uses_bash_source`).
+- Scripts 00 to 07 need no SLURM job. Run them on a login or interactive node. The Phase B
+  jobs J0, J1 and J2 (`jobs/*.sh`) need a GPU. They write to `${SCRATCH:?}` and copy the
+  results to /bigdata (see "Phase B" below).
 - Data never go into the repository.
 
 ## Work directory and variables
@@ -215,3 +217,31 @@
 | Cneo_JEC21_GOA | 32 (0) |
 | Cneo_CRYD1 | 32 (0) |
 | Umay_MYCMD | 62 (10) |
+
+## Phase B: features and embeddings (D2, D3)
+
+The plan is `docs/superpowers/plans/2026-10-01-step1-features-embeddings.md`. `COLUMNS.md`
+lists every Phase B output. All Phase B outputs go to `$STEP1_WORKDIR/phaseb/`.
+
+| Step | Where | Command | Output |
+|---|---|---|---|
+| 05 | login node | `$PY 05_prepare_sequences.py` | members, unique sequences, FASTA |
+| J0 | exfab GPU | `sbatch ... jobs/j0_pilot.sh` | `j0/throughput.json`, `j0/gpu_cpu_diff.json` |
+| 06 | login node | `$PY 06_plan_embedding.py` | `chunk_plan.tsv`, `job_plan.json` |
+| J1 | exfab GPU | `sbatch ... jobs/j1_features.sh` | `signalp/`, `predgpi/` |
+| J2 | exfab GPU | `sbatch --array ... jobs/j2_embed.sh` | `emb/<model>/<chunk_id>.npy` |
+| 07 | login node | `$PY 07_build_features.py` | `features.tsv.gz`, `feature_coverage.tsv` |
+| assemble | login node | `$ENV_PY jobs/assemble_embeddings.py` | `emb/<model>.nterm.npy`, `.cterm.npy` |
+
+- The `jobs/` Python scripts run with the conda env Python
+  (`/rhome/jstajich/.conda/envs/adhesionPred/bin/python`) and
+  `PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/jobs`.
+  `jobs/predgpi_scores.py` runs with the Python of `module load predgpi/202001` (3.9).
+- Every sequence is embedded and scored once: members with the same `seq_sha256` share one row.
+- The chunk size and the J2 job count come from the measured J0 throughput (formula in
+  `chunk_plan.py`). `06_plan_embedding.py --rate` makes a dry plan from an assumed rate and
+  records `rate_source: assumed`.
+- J1 and J2 resume: a finished part or chunk is skipped after its hash is checked.
+- Embedding tests need torch and fair-esm:
+  `PYTHONPATH=src /rhome/jstajich/.conda/envs/adhesionPred/bin/python -m pytest tests/step1_compare -q`.
+  With `/usr/bin/python3.12` those test files are skipped.
```

- [ ] **Step 5: Run the whole suite**

Run: `$PY -m pytest tests/step1_compare -q`
Expected (prototype, 2026-10-01, `PROJ_ROOT` and `STEP1_WORKDIR` unset): `274 passed, 6 skipped`.
Run: `module load predgpi/202001 && PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q`
Expected (prototype): `293 passed` (220 Phase A tests + 73 new).

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/COLUMNS.md analysis/step1_compare/README.md \
  tests/step1_compare/test_phaseb_docs.py
git commit -m "step1_compare: document the Phase B outputs and the run order"
```

---

## Run (controller, after review)

Run these steps in this order. J1 does not depend on J0 and can run at the same time. Every command below was not run by the plan author (no job was submitted).

```bash
export PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-feat   # checkout of the reviewed commit
export STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
PY=/usr/bin/python3.12
ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
S1=$PROJ_ROOT/analysis/step1_compare
LOG=$STEP1_WORKDIR/logs
cd "$S1"
```

**R0. Preflight (login node).**

```bash
git -C "$PROJ_ROOT" status --short          # expect no output
$PY -m pytest "$PROJ_ROOT/tests/step1_compare" -q
grep -h '"all_sources"' "$STEP1_WORKDIR"/{extract_log,sequence_run,d8_run,keyword_tier_run}.json
```
Accept: tests pass; four lines `"all_sources": true`.

**R1. 05 (login node, about 15 s).**

```bash
$PY 05_prepare_sequences.py
```
Accept: `members=123567 unique=69941 residues=34524441 over_1022=5541` (the prototype values on the same inputs). A different value means an input changed: read `phaseb/prepare_run.json` before you continue.

**R2. J0 (exfab GPU).**

```bash
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR" \
  -o "$LOG/j0.%j.log" -e "$LOG/j0.%j.log" "$S1/jobs/j0_pilot.sh"
```
Accept:
```bash
$ENV_PY - <<'EOF'
import json, os
d = os.environ["STEP1_WORKDIR"] + "/phaseb/j0/"
t = json.load(open(d + "throughput.json"))
assert t["schema"] == "step1-phaseb-j0-throughput/1" and t["device"] == "cuda"
assert set(t["best"]) == {"esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D"}, t["best"]
print(t["gpu_name"], json.dumps(t["best"]), t["model_load_s"])
g = json.load(open(d + "gpu_cpu_diff.json"))
print(json.dumps(g["models"], indent=1))
EOF
```
Record in the run notes: the GPU name, the best batch size and residues/s per model, the run statuses (`oom` rows), `max_abs_diff`, `min_cosine` and `repeat_identical_on_a` per model. No threshold exists for the GPU-CPU difference (spec 7 says "max difference"); report the numbers.

**R3. 06 (login node).**

```bash
$PY 06_plan_embedding.py
```
Accept: `rate_source=J0`; note `chunks`, `chunk_residues`, `n_jobs`, `time_minutes`.

**R4. J1 (exfab GPU; independent of J0).**

```bash
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR" \
  -o "$LOG/j1.%j.log" -e "$LOG/j1.%j.log" "$S1/jobs/j1_features.sh"
```
Accept: the log ends with `J1 done: SignalP status 0, PredGPI status 0`; `ls $STEP1_WORKDIR/phaseb/signalp/part_00{0..7}/prediction_results.txt.gz $STEP1_WORKDIR/phaseb/predgpi/part_00{0..7}.tsv.gz` lists 16 files. If the job is killed, resubmit the same command; finished parts are skipped.

**R5. J2 (exfab GPU, array).**

```bash
N=$($PY -c "import json;print(json.load(open('$STEP1_WORKDIR/phaseb/job_plan.json'))['n_jobs'])")
T=$($PY -c "import json;print(json.load(open('$STEP1_WORKDIR/phaseb/job_plan.json'))['time_minutes'])")
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR",J2_JOB_COUNT="$N" \
  --array=0-$((N - 1)) --time="$T" -o "$LOG/j2.%A_%a.log" -e "$LOG/j2.%A_%a.log" "$S1/jobs/j2_embed.sh"
```
Accept: every array task log ends with `exit 0`; the number of `.json` files in `phaseb/emb/<model>/` equals `chunks` in `job_plan.json` for each model. A killed task: resubmit the same command; its finished chunks print `skip ... (done, hash verified)`.

**R6. 07 (login node).**

```bash
$PY 07_build_features.py
```
Accept: exit 0; every coverage line has `no_signalp=0 no_predgpi=0`; `features_run.json` has `members` 123567 and `unique_sequences` 69941. Cross-check with the earlier tracked run (same tool, same FASTA, same mode): the number of `Cimm_RS_proteome` members with `sp_prediction == "SP"` should be 460 (`analysis/cocci_repeats/signalp_summary.tsv`, 460 of 9,910). Report any difference; it is not a STOP.

```bash
$PY - <<'EOF'
import gzip, csv, os
p = os.environ["STEP1_WORKDIR"] + "/phaseb/features.tsv.gz"
rows = [r for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t") if r["set_id"] == "Cimm_RS_proteome"]
print(len(rows), sum(r["sp_prediction"] == "SP" for r in rows))
EOF
```

**R7. Assembly (login node, numpy).**

```bash
PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/jobs" $ENV_PY "$S1/jobs/assemble_embeddings.py"
```
Accept: exit 0 and the printed shapes `esm2_t6_8M_UR50D {'nterm': [69941, 320], 'cterm': [5541, 320]}` and `esm2_t12_35M_UR50D {'nterm': [69941, 480], 'cterm': [5541, 480]}` (if 05 gave the R1 numbers). The script itself checks: every chunk present and hash-verified; each matrix row filled exactly once; float32; no NaN or inf; row counts equal the unique and >1,022 aa counts. `chunk_manifest.tsv` has 2 x `chunks` rows.

**R8. Resumability and determinism on the real run.**

- Resubmit array task 0 of R5 (`--array=0`). Accept: its log shows only `skip` lines and `"done": 0`.
- GPU run-to-run determinism: `repeat_identical_on_a` and `repeat_max_abs_diff_on_a` in `j0/gpu_cpu_diff.json` (R2). CPU determinism and kill-and-resume are covered by `test_corrupted_chunk_is_recomputed_with_the_same_hash` and `test_kill_and_resume_gives_the_same_arrays`.

## Verification done while writing this plan

A prototype of every file was built in a scratch copy of the worktree (outside the repository) and the plan's code blocks are those files.

| Check | Result |
|---|---|
| `$PY -m pytest tests/step1_compare -q` (stdlib, `PROJ_ROOT`/`STEP1_WORKDIR` unset) | 274 passed, 6 skipped |
| `module load predgpi/202001; PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q` | 293 passed (Python 3.14.2) |
| New tests per file | seqwindow 10, phaseb_prepare 9, chunk_plan 10, phaseb_embed 9, phaseb_assemble 2, phaseb_pilot 5, feature_parsers 12, phaseb_features 8, phaseb_jobs 7, phaseb_docs 1 (73) |
| Cumulative after each task (full run) | 230, 239, 249, 258, 260, 265, 277, 285, 292, 293 |
| `pre-commit run --all-files` on the prototype | ruff and ruff-format pass |
| 05 on the production Phase A outputs | 123,567 members, 69,941 unique, 5,541 > 1,022 aa |
| 06 dry plan (`--rate 29500`) on the real unique set | 5 chunks, 8,850,000 residues per chunk, 1 job, 75 min |
| J0 pilot on CPU, 20 real truth sequences | schema written; 8M 364.3 and 35M 133.4 residues/s |
| PredGPI wrapper versus CLI, 1,000 S288C proteins | 0 mismatches |
| J1 script with stub SignalP and real PredGPI, 20 S288C proteins, 3 parts | exit 0; rerun skipped 6 parts; 07 coverage `no_signalp=0 no_predgpi=0` |
| Not verified | any GPU run (J0, J1, J2 on gpu12); the conda env torch on the gpu12 driver (the SignalP env with the same CUDA 12.8 build ran there in job 29280458) |

## Self-review

1. **Spec coverage.** D2: SignalP 6 and PredGPI job (Task 9, J1), parsers (Task 7), feature table for all truth species and *C. immitis* RS (Tasks 2, 8). D3: embedding job for 8M and 35M at layer 6 with residue-mean pooling through the package (Task 4), C-terminal windows (Tasks 1, 3, 5), throughput JSON (Task 6). Spec 7: `test_cterm_window_takes_last_1022` (Task 1), `gpu_cpu_diff.py` (Task 6). Spec 8: J0 to J2 on exfab with `${SCRATCH:?}`, gzip, `PROJ_ROOT` from the environment (Task 9); J3 is out of scope. Q10: separate C-terminal matrices, M8-C/M35-C use `window_matrix` (Task 5). Not covered by design: training, evaluation, gates, rule fitting (later phase).
2. **Placeholder scan.** No TBD or "similar to" steps; every code step holds a whole file or an exact diff.
3. **Type consistency.** `Chunk.members_sha256`, `chunk_plan.read_plan`, `embed_store.load_chunk -> (arr, meta)`, `job_plan.json["batch_size"][model]`, `throughput.json["best"][model]["residues_per_s"]` are used with the same names in Tasks 3 to 6 and 9. `emb_row`/`emb_cterm_row` in 07 equal `row`/`cterm_row` of 05, which the assembly uses as matrix rows.
4. **Review Focus.** Each of the five lines has a test in its owning task (listed in the section).

## Spec issues found while planning

1. Spec 8 says jobs copy results "to `analysis/step1_compare/`". That folder is in the repository, and the Phase A README says data never go into the repository. This plan copies to `$STEP1_WORKDIR/phaseb/`.
2. Spec 8 says J1 is "sized to 1 to 1.5 h". With the measured rates (SignalP 6 about 219 proteins/s on gpu12; PredGPI 31 proteins/s per core), all 69,941 sequences take about 10 minutes. One job is the minimum; it cannot reach 1 h.
3. Spec 5 lists a "GPI score" as an input of the hybrid H. The PredGPI CLI writes only a class (score 1.0, 0.70, 0.55, or a non-GPI line). This plan adds `jobs/predgpi_scores.py`, which writes the PredGPI FPR and SVM output with the CLI's classes (0 mismatches on 1,000 proteins).
4. Spec 8 assumes 8M and 35M are "at least" as fast as ESM-2 150M at 60 proteins/s (review 9.1). That number came from bf16, HuggingFace transformers and token-budget batching; `get_esm_embeddings` uses fp32 fair-esm with a fixed number of sequences per batch. The assumption cannot be checked before J0.
5. Spec 8 gives *S. pombe* 5,130 proteins (UniProt count). The PomBase `peptide.fa.gz` that Phase A uses has 5,126 records. The S288C (6,722), H99 (7,427), JEC21 (6,740), Af293 (9,647) and *C. immitis* RS (9,910) counts match the files.
6. Spec 7 names `gpu_cpu_diff.py` but gives no acceptance threshold. This plan reports the numbers and applies none.

## Assumptions

- The input universe includes whole proteomes (design section); the owner can remove a set row.
- J0 samples only `truth_sequences.tsv.gz` (as the brief says); the J2 length mix differs a little.
- `CHUNK_SECONDS = 300`, `SAFETY = 1.25`, `TARGET_SECONDS = 4500` are choices, not measurements.
- gpu12 CPU cores run PredGPI at the c01 rate (31 proteins/s per core).
- The adhesionPred env torch 2.10.0+cu128 runs on the gpu12 driver (not verified; J0 stops if CUDA is not available).
- SignalP 6 gives the same call for a sequence in any batch, so the R6 cross-check against 460 is expected to match (not verified).
- Embeddings are stored float32 and not compressed (binary data; compression not measured).
