# Class 2b is not one structural class

**2026-09-29.** Continues `docs/reports/2026-09-28-pra-cfem-fold-survey.md`. That report found
that Ag2/PRA's only folded part is a CFEM domain, and flagged that it had **no structural
confirmation** because Foldseek was not installed. This report installs Foldseek, runs the
structural comparison, and extends it to the other class 2b members.

Reproduce with `analysis/class2b_structure/run.sh`.

## Summary

1. The class 2b positive set is **2 to 6 proteins**, depending on how the class is drawn. That
   is too few to train or validate a classifier. No classifier was built.
2. Class 2b is **not one fold**. Its members carry at least three unrelated structures.
3. The CFEM assignment for Ag2/PRA and PRA2 is now **structurally confirmed** — they superpose
   on the *Candida albicans* hemophores Csa2 and Rbt5 at TM-score 0.79–0.81 with 25–31%
   sequence identity.
4. **PRA3 is not CFEM.** Its only folded part is 38 residues, and it does not superpose on
   Ag2/PRA (TM 0.33) or on the CFEM references (TM 0.28–0.35).
5. **CalA does not superpose on any other class 2b member.** Its own fold is confirmed
   experimentally: TM 0.83 / LDDT 0.74 to the *Magnaporthe* elicitor MoHrip2 (PDB `5fid`).
6. The class 2a control (SOWgp) has **zero** residues above pLDDT 70. Class 2a and class 2b
   are different in kind, not in degree.
7. Methodological: a **whole-PDB Foldseek search misses** the Ag2/PRA–Csa2 relationship that
   direct pairwise TM-align finds. The 3Di prefilter drops true neighbours of ~60-residue
   queries. For a class defined by being small, "no Foldseek hit" is not evidence of no fold
   relationship. §6.1.

The practical consequence: the `docs/TOOL-ARCHITECTURE.md` row "2b. Small receptor-binding
invasins" describes a residual bucket, not a mechanism class. It should be split.

## 1. How many class 2b positives exist

Counted from `data/curated/adhesins/adhesins.tsv` (110 rows labelled `adhesin`):

| filter | count |
|---|---|
| labelled adhesin | 110 |
| labelled adhesin and <= 260 aa | **12** |
| of those 12: class I/II hydrophobins (class 2c, solved by HMM) | 4 |
| of those 12: flocculin or Hyr1 fragments, PA14 (class 2a families) | 4 |
| of those 12: plausibly class 2b | **2** (CalA `Q4WXJ1`, Ag2/PRA `A0A0E1RVD3`) |

PRA2 (`Q6K1L8`) and PRA3 (`Q2TVJ9`) are not in the curated adhesin set at all; they are
antigens. AGA2 (`P32781`, 87 aa) and PGA1 (`Q5ACL7`, 132 aa) could be argued in. So the
positive set is **2 named adhesins, or 6 if drawn as generously as the evidence allows.**

Those 6 span five different Pfam families: PF04681, PF05730, PF28404, PF17366, PF17056.

**A classifier cannot be built on this and should not be attempted.** Six positives cannot
support homology-grouped cross-validation, cannot estimate a precision, and cannot be split
into train and test. This is stated as a blocker, not worked around.

That is why the step taken here is a fold survey, and why its output is a statement about
class definition rather than a model.

## 2. Method

No structure prediction was needed. Every candidate already has an AlphaFold DB model, so no
GPU job was run.

| step | what |
|---|---|
| models | AlphaFold DB v6 mmCIF for 10 accessions |
| confident core | residues with pLDDT >= 70, in segments of >= 5 residues, written as PDB |
| comparison | Foldseek 9-427df8a, `--alignment-type 1` (TM-align), over the cores only |
| threshold | TM-score >= 0.5 for same fold, read together with LDDT and alignment length |

Comparing only the confident cores is the point. pLDDT below 70 is not a poor fold; it is not
a fold. An alignment that runs through such residues is not evidence of homology.

## 3. Model confidence

| protein | role | len | mean pLDDT | frac pLDDT>=70 | confident core |
|---|---|---|---|---|---|
| CalA (`Q4WXJ1`) | 2b | 177 | 86.6 | 0.82 | 145 aa, 33–177 |
| Ag2/PRA (`Q6QJA6`) | 2b | 194 | 62.0 | 0.33 | 64 aa, 23–86 |
| Ag2 (`Q12295`) | 2b | 194 | 60.7 | 0.33 | 63 aa, 23–85 |
| PRA2 (`Q6K1L8`) | 2b | 124 | 67.1 | 0.54 | 67 aa, 23–89 |
| **PRA3 (`Q2TVJ9`)** | 2b | 153 | 59.5 | **0.27** | **38 aa, 9–46** |
| Csa2 (`Q5A0X8`) | CFEM ref (PDB 4Y7S) | 147 | 87.7 | 0.77 | 113 aa, 33–145 |
| Rbt5 (`Q59UT4`) | CFEM ref | 241 | 70.4 | 0.47 | 112 aa, 25–136 |
| **SOWgp (`Q8NK60`)** | 2a control | 328 | **30.4** | **0.00** | **none** |
| RodA (`P41746`) | 2c control | 159 | 65.8 | 0.52 | 80 aa, 5 segments |
| CTS1 (`Q1E3R8`) | enzyme control | 427 | 93.9 | 0.92 | 392 aa, 36–427 |

Two things to read here.

**SOWgp has no folded part at all.** Not one residue reaches pLDDT 70. The class 2a repeat
proteins and the class 2b small proteins are not two points on a size axis; one class has a
fold and the other does not. This is consistent with the repeat-coverage measurement in
`docs/TOOL-ARCHITECTURE.md` and gives it a structural reading.

**PRA3's folded part is 38 residues.** Residues 9–46 carry 7 cysteines (18.4% Cys) at
pLDDT 80–90. Residues 47–140 are Pro/Thr/Glu low-complexity at pLDDT 35–50. Residues 141–153
are the hydrophobic GPI signal. PRA3 is a small disulfide knot on a disordered stalk, not a
globular domain.

## 4. All-vs-all TM-score over the confident cores

Normalised by query length. Rows are queries.

| | CalA | Ag2/PRA | Ag2 | PRA2 | PRA3 | Csa2 | Rbt5 | RodA | CTS1 |
|---|---|---|---|---|---|---|---|---|---|
| **CalA** | 1.000 | 0.223 | 0.222 | 0.228 | 0.144 | 0.302 | 0.280 | 0.194 | 0.362 |
| **Ag2/PRA** | 0.408 | 1.000 | 0.983 | 0.984 | 0.235 | **0.807** | **0.803** | 0.295 | 0.518 |
| **Ag2** | 0.406 | 0.999 | 1.000 | 0.984 | 0.242 | **0.797** | **0.793** | 0.298 | 0.506 |
| **PRA2** | 0.393 | 0.941 | 0.926 | 1.000 | 0.231 | **0.799** | **0.797** | 0.294 | 0.525 |
| **PRA3** | 0.403 | 0.334 | 0.336 | 0.335 | 1.000 | 0.284 | 0.345 | 0.340 | 0.382 |
| **Csa2** | 0.362 | 0.486 | 0.475 | 0.502 | 0.150 | 1.000 | 0.987 | 0.212 | 0.445 |
| **Rbt5** | 0.362 | 0.491 | 0.480 | 0.506 | 0.165 | 0.996 | 1.000 | 0.240 | 0.443 |
| **RodA** | 0.298 | 0.249 | 0.254 | 0.260 | 0.224 | 0.299 | 0.297 | 1.000 | 0.314 |
| **CTS1** | 0.168 | 0.108 | 0.112 | 0.103 | 0.064 | 0.156 | 0.154 | 0.093 | 1.000 |

### 4.1 Which of those numbers are real

TM-score alone is misleading when a short query is aligned into a long target. Reading LDDT
alongside separates the two cases cleanly:

| pair | TM | LDDT | seq id | verdict |
|---|---|---|---|---|
| Ag2/PRA — Csa2 | 0.807 | **0.77** | 0.26 | **real superposition** |
| Ag2/PRA — Rbt5 | 0.803 | **0.77** | 0.28 | **real superposition** |
| PRA2 — Csa2 | 0.799 | **0.76** | 0.28 | **real superposition** |
| Ag2/PRA — PRA2 | 0.984 | 0.97 | 0.81 | same family, as expected |
| Ag2/PRA — CTS1 | 0.518 | **0.41** | 0.01 | artifact: 64 aa query into a 392 aa TIM barrel |
| Ag2/PRA — CalA | 0.408 | **0.34** | 0.01 | artifact, same cause |
| PRA3 — CalA | 0.403 | **0.37** | 0.01 | artifact, same cause |

The three genuine relationships have LDDT 0.76–0.77. Every spurious pair sits at 0.33–0.41.
**A TM-score of 0.5 is not sufficient evidence for a small query; report LDDT with it.** This
is the "small compact proteins produce misleading structure results" check, and it fired.

### 4.2 Composition check

Low complexity does not explain any of the hits. Every core scores 0.00–0.07 on a 12-mer
low-complexity fraction. Cysteine content is high for the CFEM group (7–13%) and highest for
PRA3 (18.4%), which is expected for disulfide-stabilised secreted domains and is not what is
driving the Csa2/Rbt5 superposition — CalA at 3.4% Cys and RodA at 10.0% Cys do not superpose
on anything.

## 5. What this says about class 2b

**Ag2/PRA and PRA2 carry the CFEM hemophore fold.** This confirms the Pfam-based claim in the
2026-09-28 report with actual structural superposition, which that report explicitly could not
do — against the AlphaFold models here, and against the experimental crystal structure in
§6.1. Csa2 and Rbt5 are characterised heme-binding hemophores of *C. albicans*. So the folded
domain of Ag2/PRA is, structurally, a hemophore domain. Whether it binds heme was not tested
here and is not claimed.

That matters for the class definition. Class 2b is named "receptor-binding invasins". The
only 2b fold this project can now assign with confidence is a **hemophore** fold, which is an
iron-acquisition architecture, not a host-receptor-binding one. The name asserts a mechanism
the structure does not support.

**PRA3 is a separate thing.** 38 folded residues, 7 cysteines, no superposition on CFEM or on
Ag2/PRA. Pfam calls it PF28404 (`ARB_05178`, "uncharacterised secreted protein"), which is
consistent. PRA3 is also the best species-specificity candidate in
`docs/reports/2026-09-27-coccidioides-antigen-findings.md`, so its structure being unknown is a
real gap, not a curiosity.

**CalA is a third thing.** 145 confidently folded residues, a Bys1 domain (PF04681) with a
thaumatin-like fold,
no structural relationship to the PRA proteins. Its fold is confirmed against experimental
coordinates in §6.

So "class 2b" as written groups a hemophore-fold domain, a cysteine knot and a thaumatin
domain. They have in common that they are short and that the repeat detector misses them.
That is a property of the detector, not of the proteins.

## 6. Foldseek against the experimental PDB

Foldseek PDB database, 2.2 GB, downloaded 2026-09-29. Two passes: 3Di+AA (fast, E-value) and
TM-align. Best hit per protein:

| query | top PDB hit | what it is | E (3Di) | TM | LDDT |
|---|---|---|---|---|---|
| CTS1 (control) | `1ll7` | endochitinase CTS1, *C. immitis* — itself | 8.3e-87 | 0.999 | 1.00 |
| Csa2 (control) | `4y7s` | CFEM protein Csa2 — itself | 7.3e-20 | 0.965 | 0.98 |
| Rbt5 (control) | `4y7s` | CFEM protein Csa2 | 3.1e-15 | 0.963 | 0.97 |
| RodA (control) | `6gcj` | hydrophobin RodA, *A. fumigatus* — itself | 8.0e-09 | 0.943 | 0.89 |
| **CalA** | `5fid` | **elicitor MoHrip2, *Magnaporthe oryzae*** | **2.3e-09** | **0.795** | **0.70** |
| **Ag2/PRA** | `7ard` | ribosomal fragment, *Polytomella* | 6.6 | 0.311 | 0.42 | 
| **PRA2** | `8yzs` | human protein | 7.3 | 0.415 | 0.39 |
| **PRA3** | `8vu8` | wheat germ agglutinin domain D | 7.9 | 0.422 | 0.48 |

All four controls recover their own experimental structure. The search works.

**CalA's fold is confirmed experimentally.** Its top hit is the *Magnaporthe* elicitor MoHrip2
(TM 0.83 / LDDT 0.74 against the experimental coordinates), followed by plant thaumatins
(`2d8p`, `3x3s`, TM 0.79 / LDDT 0.66). CalA is a thaumatin-family protein, as the curated
annotation says, and its nearest structural relative is a fungal elicitor rather than a plant
thaumatin.

**Ag2/PRA, PRA2 and PRA3 return nothing above noise.** Best E-values 3.6–8.6, LDDT 0.38–0.48.

### 6.1 The whole-PDB search misses a relationship the pairwise comparison finds

This is the important methodological result. The whole-PDB search returns **no hit to `4y7s`**
for Ag2/PRA or PRA2 — in either mode — while a direct pairwise TM-align against the same
experimental structure, with the prefilter disabled (`--exhaustive-search 1`, script `06`),
gives:

| query | reference | TM | LDDT | seq id | aln |
|---|---|---|---|---|---|
| **Ag2/PRA** | **`4y7s` chain B (Csa2, experimental)** | **0.807** | **0.77** | 0.26 | 66 |
| **Ag2** | `4y7s` chain B | 0.797 | 0.76 | 0.25 | 65 |
| **PRA2** | `4y7s` chain B | 0.799 | 0.76 | 0.28 | 68 |
| CalA | `5fid` chain A (MoHrip2, experimental) | 0.832 | 0.74 | 0.22 | 150 |
| PRA3 | `6gcj` (RodA) — its best of four | 0.477 | 0.47 | 0.04 | 107 |
| Csa2 | `4y7s` chain B | 0.981 | 0.99 | 1.00 | 111 |
| CTS1 | `1ll7` chain A | 0.999 | 1.00 | 1.00 | 392 |

So **Ag2/PRA's CFEM assignment is now confirmed against an experimental crystal structure**,
which is exactly what `2026-09-28-pra-cfem-fold-survey.md` said it could not do.

And the 3Di prefilter drops the true neighbour of a 64-residue query. For class 2b — which is
*defined* by being small — a whole-database Foldseek search is not a reliable way to find fold
neighbours. Targeted pairwise comparison against a curated reference panel is. Any future 2b
fold search in this project should do both, and should not read "no Foldseek hit" as "no fold
relationship".

PRA3 still has no fold assignment. Its best pairwise score is 0.477 to the RodA hydrophobin,
below the 0.5 threshold and with LDDT 0.47, so no call is made. The whole-PDB search's weak
`8vu8` hit (wheat germ agglutinin domain D, a 43-aa hevein-like 8-cysteine knot) is
compositionally plausible for a 38-aa 7-cysteine core but is not evidence at E = 7.9.

## 7. What was not done, and why

- **The Coccidioides CFEM paralogs were not included.** The 2026-09-28 report found 7 CFEM
  proteins per *Coccidioides* genome in Fungi_5k. Those are Fungi_5k IDs with no UniProt
  accession, so they have no AlphaFold DB model and would need prediction
  (`alphafold/3.0.2` or `boltz/2.2.1` are available as modules). That is the obvious extension
  and is a GPU job, not a login-node job.
- **CalA's accession is not independently verified.** `data/curated/adhesins/README.md`
  already flags that the CalA accession was assigned from a locus tag recalled from memory.
  `Q4WXJ1` is `AFUA_3G09690`, 177 aa, annotated "Extracellular thaumatin domain protein,
  putative" — consistent with the curated description, but I did not resolve the locus tag
  from the primary paper (PMID 27841851). UniProt has no protein with gene name `calA` in
  *A. fumigatus*. Treat the CalA row as unverified.
- **No heme binding, no interface prediction, no electrostatics.** A binding-interface
  predictor was considered and rejected: with 2–6 positives there is nothing to validate it
  against, so its output would be unfalsifiable.

## 8. Recommended change to the architecture document

Replace the single class 2b row with what the structures support:

| revised class | members here | fold | right tool |
|---|---|---|---|
| CFEM-domain surface proteins | Ag2/PRA, PRA2, + ~7 per Onygenales genome | CFEM hemophore, confirmed | **PF05730 HMM**, like class 2c |
| small Cys-knot secreted proteins | PRA3 | 38 aa knot, uncharacterised | needs its own structural characterisation |
| Bys1-domain invasins | CalA | Bys1 (PF04681), thaumatin-*like* fold | **PF04681 HMM** |

Two of the three are HMM problems, not ML problems — the same conclusion standing rule 6 in
`docs/HANDOFF-HPCC.md` records for hydrophobins. The remaining open structural question is
PRA3, and it is a question about one protein, not a classification task.

---

## Addendum, 2026-09-30

Two corrections and one resolution, made when this report was reviewed.

**1. CalA's accession is resolved.** §7 listed `Q4WXJ1` as unverified. It is now verified from
the primary paper. Liu H, Lee MJ, Solis NV, Phan QT, Swidergall M, Ralph B, Ibrahim AS,
Sheppard DC, Filler SG. 2016. *Aspergillus fumigatus* CalA binds to integrin α5β1 and mediates
host cell invasion. *Nat Microbiol* 1:16211. https://doi.org/10.1038/nmicrobiol.2016.211
(PMID 27841851, PMC5495140). Its Methods state:

> "A split marker strategy was used to disrupt the 534 bp *calA* (Afu3g09690) protein coding
> sequence"

Afu3g09690 = AFUA_3G09690 = `Q4WXJ1`. Two independent checks agree: 534 bp is 177 aa plus a
stop codon, which matches `Q4WXJ1`'s length exactly; and the fold assignment in §6 is to an
experimental structure. **Every CalA conclusion in this report stands.**

**2. PF04681 is Bys1, not thaumatin.** This report called PF04681 a "thaumatin domain". PF04681
is `Bys1`, "*Blastomyces* yeast-phase-specific protein" (InterPro IPR006771). Thaumatin proper
is **PF00314**. UniProt's protein *name* for `Q4WXJ1` is "Extracellular thaumatin domain
protein, putative", and the structural neighbours found in §6 are thaumatin-like, so the
accurate phrasing is **Bys1 family, thaumatin-like fold**. The structural result is unchanged;
only the family label was wrong. §5 and §8 are corrected above.

**3. A specificity control that was missed.** UniProt returns **three** PF04681 proteins in the
Af293 proteome, not one: `Q4WXJ1` (calA, 177 aa), `Q4WBB5` (AFUA_8G01710, 167 aa, "Antigenic
thaumatin domain protein") and `Q4WFZ4` (AFUA_3G00510, 290 aa). These are calB and calC. They
are **not** in the curated set and were **not** tested here. A PF04681 HMM as recommended in §8
would return all three, so the recommendation "PF04681 HMM" is only as good as its behaviour on
calB and calC — which is unmeasured. The same question applies to PF05730 and the ~7 CFEM
paralogs per genome noted in §7. **An HMM that recovers the whole domain family is not by
itself a predictor of the phenotype.**

Sources for this addendum:

- Liu H et al. 2016. *Nat Microbiol* 1:16211. https://doi.org/10.1038/nmicrobiol.2016.211
- UniProt entry `Q4WXJ1`. https://rest.uniprot.org/uniprotkb/Q4WXJ1.txt
- Pfam/InterPro PF04681 (Bys1, IPR006771). https://www.ebi.ac.uk/interpro/entry/pfam/PF04681/
