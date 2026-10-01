# class2b_structure — is "class 2b" one structural class?

`docs/TOOL-ARCHITECTURE.md` defines **class 2b** as "small receptor-binding invasins"
(CalA, Ag2/PRA, PRA3) and says they are "not achievable from sequence; needs structure".

This directory tests the assumption behind that row. A class is only useful as a modelling
target if its members share something. The claim under test is:

> The class 2b exemplars share a compact fold with one binding interface.

The test is structural, not sequence-based: AlphaFold DB models, confident-core extraction by
pLDDT, and Foldseek comparison — first among the exemplars, then against the experimental PDB.

**This directory does not train a classifier, and cannot.** The class has 2 to 6 members,
depending on how it is drawn. See
`docs/reports/2026-09-29-class2b-structure.md` §1 for why that number rules classification out.

## Files

| file | what it is |
|---|---|
| `candidates.tsv` | the candidate set: 5 class 2b accessions plus 5 controls, with the role each control plays |
| `01_fetch_structures.py` | downloads the AlphaFold DB v6 mmCIF and the UniProt FASTA for each accession |
| `02_confidence_profile.py` | per-residue pLDDT; writes the pLDDT>=70 "confident core" as a PDB, and a confidence summary |
| `03_foldseek_allvall.sh` | all-vs-all Foldseek in TM-align mode over the confident cores |
| `04_foldseek_pdb.sh` | Foldseek search of the cores against the experimental PDB, in 3Di+AA and TM-align modes |
| `05_summarize.py` | the result tables, including the composition artifact check |
| `06_targeted_pairwise.sh` | pairwise TM-align against experimental PDB references, prefilter disabled -- finds what the whole-PDB search misses for short queries |
| `run.sh` | runs 01-02 locally and submits 03-04 to SLURM; submit 06 the same way |

## How to run

```bash
cd <repo root>
bash analysis/class2b_structure/run.sh
# when the SLURM job finishes:
/usr/bin/python3.12 analysis/class2b_structure/05_summarize.py \
    --workdir _workdir/class2b_structure
```

Outputs go to `_workdir/class2b_structure/` (gitignored; structures and the Foldseek PDB
database are bulky and regenerable). Nothing here writes to `analysis/cocci_repeats/`.

## Environment notes

- **Foldseek is not installed on this HPCC.** `run.sh` prints the one-line install for the
  static binary into `_workdir/bin/`. No root needed.
- **The AVX2 build dies with "Illegal instruction" on the older nodes** — `c01` has no AVX2.
  This is the same trap `docs/HANDOFF-HPCC.md` records for MMseqs2. `run.sh` sends the
  Foldseek steps to `epyc` via `sbatch`.
- **No structure prediction was needed.** Every candidate already has an AlphaFold DB model,
  so no GPU job was run. If the set is extended to Fungi_5k proteins with no UniProt
  accession, that changes — `alphafold/3.0.2`, `boltz/2.2.1` and `rosettafold2/1.0` are all
  available as modules.

## Reading the numbers

- **pLDDT < 70 is not a low-quality fold, it is not a fold.** Foldseek alignments through
  such residues are not evidence of homology. Only the confident cores are compared.
- **TM-score >= 0.5** is the conventional same-fold threshold. Values are normalised by query
  length (`qtmscore`).
- **Read LDDT next to it.** In these results every genuine fold relationship has LDDT 0.70-0.77
  and every spurious one 0.33-0.48. TM-score alone at 0.5 is not sufficient when the query is
  short and the target is long.
- **A whole-database Foldseek search can miss a true neighbour of a ~60-residue query** (the
  3Di prefilter drops it). Run `06_targeted_pairwise.sh` as well before concluding there is no
  fold relationship.
- **Small cysteine-rich domains and low-complexity regions both generate structural matches
  that are not fold homology.** `05_summarize.py` prints the Cys, Pro, Ser+Thr and
  low-complexity fraction of every core so a hit can be checked against that.
