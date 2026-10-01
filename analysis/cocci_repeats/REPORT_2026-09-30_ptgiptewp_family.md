# The PTGIPTEWP repeat family: universal, secreted, *Coccidioides*-specific — and not spherule-induced

**2026-09-30.** Type protein `CIMG_04070` (*C. immitis* RS), UniProt `J3KCQ3`.

Found by `32_class2a_anchor_sweep.py` as a family that **neither periodicity detector calls**.
This report follows it up on expression, prevalence, homology, structure and antigen status.

> **Status: computational only. No experimental validation.**

---

## 1. Summary

1. **It is the mirror image of SOWgp.** SOWgp is massively spherule-induced but fails the
   antigen prevalence filter. This family passes prevalence at **99.4%** and is
   *Coccidioides*-specific, but is **not spherule-induced** (log2FC **-1.08**) and is expressed
   at only 4-17 TPM.
2. **It was already in the antigen Tier 2 list** (`analysis/cocci_antigens/TIER2_candidates.tsv`),
   found independently by orthology and prevalence. Two unrelated methods converged on it, which
   is the strongest thing about it.
3. **A composition label in the candidate table is wrong.** It is recorded as
   `ProCys_rich(SOWgp_BAD1_type)`. The protein contains **zero cysteine**. It is Pro/Thr-rich.
4. **No experimental structure exists** — AlphaFold model only.
5. **Unit count varies 2-7 across the pangenome**, so it has the copy-number variation that
   motivated the antigenic-variation question, without the spherule expression that would make
   it a phase-specific target.

## 2. The protein

```
MQFKLSIALVAALAALSEASPYRLHPRRKAVSPRRVDFVTMTPSSPYPTGTITDIPTDVPP
TGIPTEWPPTGIPTEWPPTGIPTEWPPTEWPTGWPTGWPTGIPTELPTWFPSKDVPMETLT
FTYTLGKKPSQSVVTKTITRPAVAEPSEG
```

| | |
|---|---|
| length | 151 aa |
| annotation | "hypothetical protein" (FungiDB) / "Uncharacterized protein" (UniProt) |
| UniProt | `J3KCQ3` |
| repeat period | 9 aa, anchor `IPTEWP` at positions 63, 72, 81 |
| signal peptide | **yes**, cleavage 19-20, p = 0.987 |
| composition | Pro **18.5%**, Ser+Thr **23.2%**, Gly 6.0%, **Cys 0.0%** |

**Correction carried into `class2a_candidates_general.tsv`:** the `comp_class` column calls this
family `ProCys_rich(SOWgp_BAD1_type)`. With 0% cysteine that is wrong, and it implies a
relationship to SOWgp and BAD1 that nothing here supports. The classifier's Pro/Cys rule is
evidently driven by proline alone when cysteine is absent. **This affects the composition
class of any other cysteine-free Pro-rich family in that table.**

Scope of the mislabel, now measured: **6 of the 21 candidates labelled `ProCys_rich` have
under 2% cysteine.**

| protein | len | period | %Cys | %Pro | %Ser+Thr |
|---|---|---|---|---|---|
| `CPOS1038_006679-T1` | 156 | 9 | **0.0** | 19.2 | 23.1 |
| `CIMG_07912-t26_1-p1` | 105 | 7 | 1.0 | 20.0 | 17.1 |
| `CIB10637_006868-T1` | 105 | 7 | 1.0 | 21.0 | 17.1 |
| `CIB10992_007229-T1` | 105 | 7 | 1.0 | 21.0 | 17.1 |
| `CPOS1038_004507-T1` | 105 | 7 | 1.0 | 20.0 | 18.1 |
| `CPOS3700_005938-T1` | 105 | 7 | 1.0 | 20.0 | 18.1 |

`CPOS1038_006679-T1` is this family's Cpos1038 orthologue (period 9). The other five are a
**second** cysteine-free family, period 7, 105 aa, one per proteome — another orthologue set
also mislabelled Pro/Cys-rich, and also worth following up.

So the label is wrong for 6 of 21, not systemically for all of them. The other 15 are genuinely
cysteine-rich. The consequence is that "Pro/Cys-rich (SOWgp/BAD1 type)" in this project's tables
does **not** reliably mean the protein resembles SOWgp, and the `pct_cys` column should be read
directly rather than trusting `comp_class`.

## 3. Expression: not spherule-induced

*C. immitis* RS spherule/mycelium RNA-seq (`RS1_kallisto.TPM.csv`), TPM:

| gene | mycelia | spherule 48 h | spherule 8 d | log2FC (48 h) |
|---|---|---|---|---|
| **`CIMG_04070` (this family)** | 11.9 | **5.6** | 4.1 | **-1.08** |
| `CIMG_04613` (SOWgp, control) | 13.1 | 15000.4 | 4839.2 | **+10.15** |

It goes **down** in spherules, and its absolute level is low throughout. The Tier-1 antigen
criterion is log2FC > 1, so it fails decisively. This is a clean negative and it is the main
thing that limits the family's interest as a phase-specific antigen.

## 4. Pangenome prevalence: 99.4%

`34_newfam_pangenome.sh`, anchored search over 493 pangenome proteomes.

| | |
|---|---|
| hits | 491 |
| **passing the alignment member filter** | **491 (100%)** |
| proteomes with a member | **490 / 493 = 99.4%** |
| identity to `CIMG_04070` | 94.6-100% (median 98.6%) |
| length | 94-179 aa (median 152) |

Every single hit passes the member filter. That is unusual — for SOWgp the same pipeline
rejected nothing, but for the class-2a sweep as a whole it rejected 28 of 86. This family is
tight and unambiguous.

Model quality is also good: 364 `full_length`, 125 `short_allele`, **only 2 `truncated`** and
zero split models. Contrast SOWgp, where 91 of 474 hits were model problems.

### Unit-count variation

| anchored units | strains |
|---|---|
| 2 | 80 |
| 3 | 127 |
| 4 | **243** |
| 5 | 38 |
| 6 | 2 |
| 7 | 1 |

A modal 4-unit allele with a spread from 2 to 7. **This is real copy-number variation in a
universally present, secreted repeat protein** — the same phenomenon as SOWgp's 3-6 units, in a
family with far better gene models, so it is less confounded by annotation error. If the
antigenic-variation hypothesis is to be tested on unit count, this family is a cleaner
substrate than SOWgp.

## 5. It is already a Tier 2 antigen candidate

From `analysis/cocci_antigens/cocci_antigen_ranking.tsv`:

| field | value |
|---|---|
| orthogroup | OG0003601 |
| score / antigenicity | 2.4898 |
| **prevalence** | **0.9959** (486 / 488 strains) |
| copy_number_mean / CV | 1.00 / 0.045 |
| **n_confounder_genera** | **0** |
| max_fungal_crossreact_pid | **0.0** |
| human_homolog_pid | 0.0 |
| in TIER1 | no |
| **in TIER2** | **yes** |

Two independent routes converged. The antigen pipeline reached it through orthology and
prevalence; the anchored repeat search reached it through its 9 aa repeat. Neither knew about
the other. **My anchored prevalence (99.4%) independently reproduces the pipeline's 0.9959.**

It is in Tier 2 and not Tier 1 for exactly the reason section 3 measures: Tier 1 additionally
requires spherule induction, and this family is not induced.

## 6. Homology and structure

- **UniProt `J3KCQ3`**, *Coccidioides immitis* RS. "Uncharacterized protein". No function
  annotation, no GO terms of substance.
- **No experimental structure.** No PDB cross-reference. An AlphaFold DB model exists
  (`AF-J3KCQ3`). It has **not** been examined here — given SOWgp's AlphaFold model has zero
  residues above pLDDT 70, a Pro/Thr-rich repeat protein of this kind is likely to be largely
  disordered, and that should be checked before any structural claim.
- **Confounder fungi:** the antigen pipeline reports `n_confounder_genera = 0` and
  `max_fungal_crossreact_pid = 0.0`, i.e. no ortholog in *Histoplasma*, *Blastomyces*,
  *Paracoccidioides* or *Aspergillus*. That is whole-protein orthology against four genomes.
  **A wider Onygenales scan (71 genomes) and a unit-level comparison are running separately**
  and are needed before calling this family *Coccidioides*-specific with confidence.

## 7. What this family is, and is not

**Is:** a universally present (99.4%), secreted, *Coccidioides*-specific, Pro/Thr-rich tandem
repeat protein with variable unit count (2-7), clean gene models, no known function, no
experimental structure, and an existing Tier 2 antigen placement.

**Is not:** spherule-induced, highly expressed, cysteine-rich, related to SOWgp, or structurally
characterised.

**The practical reading.** As a *phase-specific* antigen it is weak — it is not induced when the
host meets the pathogen. As a **constitutive surface marker** it is much better placed than
SOWgp: universal where SOWgp is 92%, clean models where SOWgp's are split and truncated, and
specific by the same test. And as a **substrate for the unit-count variation question** it is
the better of the two, because its copy-number variation is not confounded by gene-model
failure.

## 8. Limits

1. Expression is from one RNA-seq experiment on one strain (RS), two replicates per condition.
   A low, flat profile in RS does not exclude induction in another condition or strain.
2. Unit counts are anchored counts from assembled proteomes. Tandem arrays collapse in
   short-read assemblies, so the 2-7 spread is a floor and the low end may be artifactual — the
   same caveat that applies to SOWgp, though the near-absence of truncated models here is
   reassuring.
3. Specificity rests on the antigen pipeline's four confounder genomes. The wider scan is not
   in yet.
4. Nothing here is experimental. No protein has been detected, no serum tested.
5. The AlphaFold model was not inspected.

## 9. Files

| file | contents |
|---|---|
| `newfam_CIMG_04070.faa` | the type protein |
| `34_newfam_pangenome.sh` | 493-proteome anchored search |
| `newfam_anchored_pangenome.tsv` | 491 member hits across 490 proteomes |
