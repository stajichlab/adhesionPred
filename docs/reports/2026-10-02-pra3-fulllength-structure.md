# Full-length PRA3 (E9CRM7): class 2b structure survey re-run

Date: 2026-10-02. Branch: `pra3-fulllength`. Analysis only. No function is inferred.

## Question

The earlier survey used UniProt Q2TVJ9 (153 aa). Its AlphaFold DB model gave a 38 aa
confident core and no Foldseek fold assignment. Q2TVJ9 matches the annotated PRA3 protein
(CIMG_02492, 220 aa, signal peptide 1-19) from 0-based offset 67. The 67 residues it lacks hold
4 Cys. Hypothesis: the "no fold" result is an artifact of the truncated model. This report
re-runs the survey on the full-length AlphaFold DB model of UniProt E9CRM7 (*Coccidioides
posadasii* Silveira, 220 aa, v6).

## What was run

- Candidate list: `analysis/class2b_pra3_fulllength/candidates_fulllength.tsv`. It is a copy of
  `analysis/class2b_structure/candidates.tsv` plus one row, E9CRM7 (label `PRA3_full`). The
  original file is unchanged. Q2TVJ9 stays in the list for comparison.
- Scripts `01`, `02`, `03`, `04`, `05`, `06` in `analysis/class2b_structure/` were used
  unchanged. Candidate list and work directory were passed on the command line.
- Work directory: `_workdir/class2b_pra3_fulllength/` (not committed).
- Foldseek 427df8a6b5d0ef78bee0f98cd3e6faaca18f172d (release 9-427df8a, avx2 build).
- Database: Foldseek `PDB` (pdb100), PDB date 250101 (`pdb.version`), already present. It was
  not downloaded again. Reference structures for step 06 (4Y7S, 5FID, 1LL7, 6GCJ) were fetched
  from RCSB.
- SLURM, partition `highclock`, node hz01, 8 CPU, 16 G, constraint cleared:
  - job 29341772: steps 03 and 04 (all-vs-all; whole-PDB search in 3Di+AA and TM-align modes),
    run time 32 s.
  - job 29341773: step 06 (pairwise TM-align, prefilter disabled), run time 1 s by `sacct`.
- Step 05 was run on the login node to build the tables.
- The UniProt FASTA fetch for E9CRM7 returned no file. The entry E9CRM7 is inactive in UniProt
  (deleted: "Not part of a reference proteome"). Its sequence is kept in UniParc as UPI0001A7D0B4.
  That 220 aa sequence is identical to the residue list of the AlphaFold mmCIF (checked
  2026-10-02). The AlphaFold model is therefore a model of a real, current sequence.

## Confidence profile

Cutoff pLDDT >= 70, as in the earlier survey. Positions are 1-based on the 220 aa protein.

| | E9CRM7 (full length) | Q2TVJ9 |
|---|---|---|
| length (aa) | 220 | 153 |
| Cys in whole model | 12 | 8 |
| mean pLDDT, whole model | 62.2 | 59.5 |
| pLDDT >= 70 core, residues | 67 | 38 |
| core ranges | 49-64 and 68-118 | 9-46 |
| Cys in core | 11 | 7 |
| mean pLDDT of core | 90.4 | 85.0 |
| core Cys, % | 16.4 | 18.4 |
| core Pro, % | 10.4 | 2.6 |

Notes:

- Core is two segments. Residues 65-67 (Y, P, S) have pLDDT 66.5, 59.2 and 66.5. Both segments
  are written to one core PDB.
- Q2TVJ9 core 9-46 maps to E9CRM7 residues 76-113. This lies inside the new core.
- The new core adds 29 residues over the old core. It extends to 49-75 and 114-118.
- The one Cys outside the core is Cys43, pLDDT 66.6.
- Core Cys positions (E9CRM7 numbering): 52, 56, 63, 69, 78, 79, 85, 91, 96, 97, 106.
- Signal peptide (1-19): mean pLDDT 47.8. Not part of the mature protein.
- Mature region (residues 20-220, 201 aa): mean pLDDT 63.6, median 55.1. 68 residues have
  pLDDT >= 70, 127 have >= 50, 74 have < 50. The 68 residues at >= 70 are the 67-residue core
  plus residue 46, which is isolated (a segment of 1 residue is dropped by script `02`). The
  confident core of the mature region is the same core.
- 10-residue window means: residues 21-40 about 49-51; 41-50 65.8; 51-110 80-95; 111-120
  82.6; 121-220 43-57. Residues 121-220 are Pro, Thr, Glu rich (for example
  `PEPTETEVEPTPTEEPTTPTIIP`) with low pLDDT.
- Because the core lies entirely in the mature region, a separate run on the mature core would
  repeat the same 67 residues. It was not done.

## Foldseek results for E9CRM7 (core, 67 aa)

TM is normalised by the query length (qTM). LDDT is the Foldseek LDDT of the alignment.
In the earlier survey, genuine relationships sat at LDDT 0.76-0.77 and artifacts at 0.33-0.48.

### All-vs-all over the cores (step 03)

| target | qTM | LDDT |
|---|---|---|
| CTS1 (392 aa) | 0.387 | 0.383 |
| Rbt5 | 0.330 | 0.289 |
| CalA | 0.332 | 0.281 |
| Csa2 | 0.308 | 0.268 |
| Ag2_PRA | 0.298 | 0.295 |
| Ag2 | 0.296 | 0.292 |
| PRA2 | 0.290 | 0.271 |
| RodA | 0.295 | 0.313 |
| PRA3 (Q2TVJ9 core, 38 aa) | 0.552 | 0.788 |

The last row is the same sequence region as the old core (identity 1.00 over 38 residues).
It is not an independent fold match. The reverse direction (query PRA3, 38 aa) gives TM 0.949 and
LDDT 0.945, which is expected for identical coordinates.

### Whole-PDB search (step 04)

TM-align mode returned 1 hit. 3Di+AA mode returned 4 hits (E-value 10 or below).

| mode | target | qTM | LDDT | E-value | aligned length | target length |
|---|---|---|---|---|---|---|
| TM-align | 2hf3-assembly1_A (*D. melanogaster*) | 0.391 | 0.361 | 0.215 | 113 | 364 |
| 3Di+AA | 8jxa-assembly1_A (*R. norvegicus*) | 0.250 | 0.340 | 3.5 | 72 | 1077 |
| 3Di+AA | 8gye-assembly1_A (*H. sapiens*) | 0.220 | 0.318 | 7.1 | 75 | 137 |
| 3Di+AA | 2l8l-assembly1_A (*M. tuberculosis*) | 0.195 | 0.285 | 5.1 | 56 | 139 |
| 3Di+AA | 6byv-assembly1_A (*H. sapiens*) | 0.180 | 0.286 | 7.5 | 58 | 121 |

The positive control CTS1 returned its own structure 1LL7 (qTM 0.999), so the search worked.
Controls Csa2, Rbt5, RodA and CalA found their known experimental structures
(4Y7S, 4Y7S, 6GCJ, 5FID) in the 3Di+AA mode (E-values 7e-20 to 2e-9).

### Targeted pairwise TM-align, prefilter disabled (step 06)

Best hit of E9CRM7 to each reference (all hits to that reference are in the work directory).

| reference | best qTM | LDDT |
|---|---|---|
| 6GCJ (RodA hydrophobin, 10 models) | 0.375 (model 8) | 0.338 |
| 4Y7S (CFEM, Csa2; chains A, B, C) | 0.368 (chain C) | 0.324 |
| 1LL7 (CTS1 chitinase) | 0.345 | 0.299 |
| 5FID (MoHrip2) | 0.307 | 0.348 |

Highest LDDT of any E9CRM7 hit in this table: 0.395 (6GCJ model 4, qTM 0.344).

For comparison, the same run gives Q2TVJ9 core (38 aa): best qTM 0.477 to 6GCJ model 2
(LDDT 0.47), 0.388 to 5FID (0.39), 0.381 to 1LL7 (0.32), 0.348 to 4Y7S (0.28). Controls in this
run: Ag2_PRA to 4Y7S 0.807 / 0.77, PRA2 0.799 / 0.76, Ag2 0.797 / 0.76, CalA to 5FID
0.832 / 0.74. These reproduce the earlier values.

## Comparison with the earlier PRA3 result

| | Q2TVJ9 (earlier) | E9CRM7 (this run) |
|---|---|---|
| core length | 38 | 67 |
| Cys in core | 7 | 11 |
| best TM to an experimental reference | 0.477 (6GCJ) | 0.375 (6GCJ) |
| LDDT of that hit | 0.47 | 0.34 |
| best TM to the CFEM structure 4Y7S | 0.348 | 0.368 |
| Foldseek whole-PDB hit with TM >= 0.5 | none | none |
| any hit at LDDT >= 0.70 | none | none |

## What this shows

- The full-length model has a larger confident core (67 aa against 38 aa) with 11 Cys.
- No hit reaches TM 0.5 in any search. The best qTM is 0.391 (whole-PDB, TM-align) and 0.375
  (targeted). All LDDT values are 0.40 or below. This is in the range of the earlier
  artifacts (0.33-0.48), not in the range of the genuine relationships (0.76-0.77).
- The full-length model does not gain a fold assignment. The hypothesis that the earlier
  null came from the truncated model is not supported by these data for this model and these
  databases and references.
- The CFEM relationship that Ag2/PRA and PRA2 show (TM about 0.80, LDDT 0.76-0.77) is not seen
  for E9CRM7 (qTM 0.308-0.368, LDDT 0.27-0.32 against Csa2/Rbt5/4Y7S).

## What this does not show

- It does not show that PRA3 has no fold. The AFDB model has mean pLDDT 62.2. The core is
  67 aa, and 66 percent of the mature protein (133 of 201 residues) is at low
  confidence. A short, Cys-rich core can fail structural search for reasons other than the
  absence of a fold.
- The search covered the Foldseek PDB (pdb100, 250101) and four reference structures. It did
  not cover the AlphaFold DB or other databases. A null here is not evidence of no
  relationship.
- The AlphaFold model is a prediction. Disulfide pairing of the 11 core Cys was not examined.
- No function is inferred for PRA3.
- One further caution: the Q2TVJ9 core is a sub-region of the E9CRM7 core (identical
  sequence), so the two results are not independent.

## Files

- `analysis/class2b_pra3_fulllength/candidates_fulllength.tsv`
- Outputs (not committed): `_workdir/class2b_pra3_fulllength/` (`cif/`, `core_pdb/`,
  `plddt_per_residue.tsv`, `confidence_summary.tsv`, `foldseek/*.tsv`, `logs/`).
- Command line used for the tables:
  `05_summarize.py --workdir _workdir/class2b_pra3_fulllength --candidates analysis/class2b_pra3_fulllength/candidates_fulllength.tsv`

## Follow-up (2026-10-02, later): AlphaFold DB search, disulfide pairing, ARB_05178

Questions from the first report: does a search against AlphaFold DB find anything, and does the
core keep a consistent disulfide pattern? Scripts: `analysis/class2b_pra3_fulllength/05_afdb_search.sh`
and `candidates_d4ali1.tsv`.

### Disulfide pairing in the model

All 12 Cys of E9CRM7 pair up. Six SG-SG contacts are 1.99 to 2.04 A: 43-56, 52-69, 63-78, 79-91,
85-96 and 97-106 (pLDDT of the Cys 67 to 97). Every Cys has exactly one partner within 3.0 A.
This is a property of the prediction. It shows that a full set of six disulfides fits the model.
It is not experimental evidence of disulfide bonds.

### Search against AlphaFold DB Swiss-Prot (v6)

Foldseek 9-427df8a, 3Di+AA and TM-align modes, query = core (67 aa) and full model. Controls ran in
the same job.

| Query | Best hit | qTM | LDDT | E-value |
|---|---|---|---|---|
| PRA3 core (E9CRM7) | D4ALI1, *T. benhamiae* ARB_05178 | 0.74 (3Di), 0.82 (TM-align) | 0.74, 0.76 | 3e-4 |
| PRA3 core, next hit | P28799 human progranulin | 0.44 | 0.48 | 5.3 |
| PRA3 full model | none above E 6 (human and mouse proteins) | 0.09 to 0.11 | 0.57 to 0.59 | 6.6 to 9.6 |
| Ag2/PRA core (CFEM control) | D4B2Q8, Q4WMA6 | 0.97 | 0.95 | 4e-6 to 2e-5 |
| RodA core (hydrophobin control) | RodA itself and *Penicillium* hydrophobins | 0.98 to 0.99 | 0.96 to 0.98 | 1e-9 to 3e-11 |

- D4ALI1 is the founding member of Pfam family PF28404 (see
  `docs/reports/2026-10-02-pf28404-family.md`). It is an uncharacterised secreted protein
  (Swiss-Prot A5178_ARTBC). So the one good hit of the PRA3 core is a member of the same family.
  The structure result agrees with the sequence result. It does not give PRA3 a fold or a function.
- Both controls find their families at LDDT 0.95 or higher, so the search works.
- The full model finds no structural neighbour. About two thirds of it is low-confidence
  (residues 121 to 220 are Pro/Thr/Glu rich).

### The ARB_05178 core as a second query

AlphaFold model of D4ALI1: 210 aa, mean pLDDT 61.2. Confident core 70 aa (residues 22-26, 51-72,
74-116). It is as low-confidence as PRA3.

- Swiss-Prot AFDB: itself only (qTM 1.0). Next hits are progranulin entries: P28797 (guinea pig,
  qTM 0.56 in TM-align mode, LDDT 0.53, E 0.3) and P28799-3 (human, qTM 0.48, LDDT 0.50, E 2.2).
  The PRA3 core also hit human progranulin P28799 (qTM 0.44, LDDT 0.48, E 5.3).
- pdb100: no hit above query TM 0.26 in 3Di mode, no TM-align hit.
- Reading: both PF28404 cores give weak hits to progranulin, which contains Cys-rich granulin
  domains. These hits are below the LDDT 0.7 level of the genuine relationships in this survey.
  Cys-rich small domains often match each other in structure search without a common origin. It is
  not evidence of homology.

### What was not done

- The large AlphaFold DB search (UniProt50-minimal, a 114 GB download) was not run. The download
  was stopped after the Swiss-Prot result. The reasons: the one good hit is the family founder, the
  cores are small and low-confidence, and the sequence search already finds the family members. It
  can be run later if a hit outside the family is wanted. It is unknown whether the "minimal"
  database has the C-alpha data that the TM-align mode needs.
- No claim is made about the function of PF28404 proteins.
