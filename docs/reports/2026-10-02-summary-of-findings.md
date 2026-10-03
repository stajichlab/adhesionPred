# Summary of findings, 2026-10-02

Scope: the step 1 measurement (surface glycoproteins), the PRA3 / PF28404 work, and the SOWgp
presence/absence work. Each section names the full report. This page adds no new result. Where a
number comes from a report that I did not re-check, the text says so.

## 1. Step 1: surface glycoprotein prediction (rule versus ML)

Full report: `_workdir/step1_compare/phasec/report.md` (not committed, regenerate with
`analysis/step1_compare/phasec/12_report.py`). Summary: `docs/model-review/STATUS.md`.

- The shipped CLI is a surface-glycoprotein detector (package `surface_glyco`). It is not an
  adhesin predictor. No trained model ships.
- Phase C compared rules (R0, R1, R2), two baselines (B0, B1) and ESM-2 models (M8, M35, the
  C-terminal-window variants, and a hybrid H) on GO-based truth. Intervals are 95% cluster
  bootstrap.
- Result with the C grid widened to 0.0001-10 (ruling C-17):

| Test set | Label | R0 recall / FPR | R2 recall / FPR | B1 recall / FPR | M8 recall / FPR |
|---|---|---|---|---|---|
| S1:all (*S. cerevisiae* + *C. albicans*) | estimate | 0.603 / 0.037 | 0.418 / 0.006 | 0.763 / 0.130 | 0.772 / 0.084 |
| Eurotiomycetes (leave one clade out) | estimate | 0.727 / 0.010 | 0.227 / 0.005 | 0.859 / 0.168 | 0.898 / 0.034 |
| Basidiomycota (16 positives) | smoke test | 0.938 / 0.083 | 0.125 / 0.017 | 0.625 / 0.283 | 0.750 / 0.133 |

- No ML candidate beats both B1 and R2 on the false-positive rate for non-secreted proteins
  (findings b and c). R2 is precise and misses more than half of the surface proteins.
- Widening the C grid changed recall and FPR by at most a few points, except the Basidiomycota
  smoke test (M8 recall 0.875 to 0.750).
- Open: the choice of rule, ML or hybrid, and the gate values, are the owner's decisions.
  Basidiomycota truth and curated GPI rows do not exist. SignalP under-calling in *C. immitis*
  is untested.

## 2. PRA3 and the PF28404 family

Full reports: `docs/reports/2026-10-02-pra3-fulllength-structure.md` and
`docs/reports/2026-10-02-pf28404-family.md`.

**PRA3 structure.**
- The full-length AlphaFold model of PRA3 (UniProt E9CRM7, 220 aa) has a confident core of 67 aa
  with 11 Cys, against 38 aa for the earlier 153 aa model. No Foldseek search of pdb100 reaches
  query TM 0.5. All LDDT values are 0.40 or lower. The hypothesis that the earlier "no fold"
  result came from the truncated model is not supported. This does not show that PRA3 has no
  fold. The model has mean pLDDT 62.2.
- UniProt E9CRM7 is an inactive entry. Its UniParc sequence matches the AlphaFold model.
- All 12 Cys of the model form six SG-SG pairs at 1.99 to 2.04 A. This is a property of the model.
- In AlphaFold DB Swiss-Prot, the PRA3 core matches one entry well: D4ALI1, *Trichophyton
  benhamiae* ARB_05178 (query TM 0.74 to 0.82, LDDT 0.74 to 0.76). D4ALI1 is the founding member
  of PF28404. So the match is within the family. The D4ALI1 core matches only itself.
- Both cores give weak hits to progranulin (LDDT 0.48 to 0.53). I read these as Cys-rich
  artifacts. They are not evidence of homology.
- The 114 GB AlphaFold DB (UniProt50-minimal) search was not run.

**PF28404 family.**
- A search of 831 proteomes gave 408 proteins in 283 proteomes. 68 of 78 Onygenales proteomes have
  at least one copy (167 proteins). The family also occurs in Eurotiales, Chaetothyriales,
  Hypocreales, Helotiales and other orders. It is not specific to *Coccidioides*. Function is
  unknown. Pfam has no clan, no GO term and no experimental structure.
- In a tree of 217 proteins (IQ-TREE, 148 trimmed columns), the four *Coccidioides* paralogs (PRA3,
  CIMG_07303, CIMG_05560, CIMG_07843) fall in four separate groups. The closest clade of each group
  holds other Onygenales genera. The split between the groups is older than *Coccidioides*.
  Support is 100 and 99 for two groups and 68 and 79 for the other two.
- *Coccidioides* has 4 copies per proteome (6 in one, 3 in one). *Uncinocarpus* has 2. *Trichophyton*
  and *Arthroderma* have 1 to 3. *Histoplasma*, *Blastomyces* and *Paracoccidioides* have 0 or 1.
  All five *Ascosphaera* proteomes have none.
- Limits: orthology is not inferred. The deep relationships between the four groups are not
  resolved. Copy numbers depend on the Pfam threshold and on the gene models. *Coccidioides*
  pangenome strains were not searched.

## 3. SOWgp presence and absence

Reports: `analysis/cocci_repeats/REPORT_2026-10-01_sowgp_depth_vs_repeats.md`,
`REPORT_2026-10-01_sowgp_longread_annotation.md`, `CIMG_04613_coverage_summary.md`. The points
below come from the commit message of PR #43. I did not re-check the numbers.

- Read depth over CIMG_04613 was measured in 559 CRAMs (mosdepth), over the whole gene and over the
  repeat array.
- The 39 strains with no SOWgp gene model have normal depth (ratio 0.92 to 1.43). So the 8%
  "absent" in the pangenome is a gene-model gap, not a deletion. Seven other strains have low
  depth.
- Depth over the repeat array does not scale with unit count as predicted (slope 0.11 for
  *C. immitis*, 0.016 for *C. posadasii*, predicted 0.25).
- In the five UArizona long-read assemblies, 3 of 6 SOWgp loci have split gene models.
- This is computational. There is no experimental validation.

## 4. What is open

- Owner decisions: step 1 choice and gates; curation of Basidiomycota truth and GPI rows.
- Step 1 work: the Phase C comparison with R0, R1 and R2 at each rule's operating point is coded
  but needs a re-run on HPCC (steps 11 and 12) to give numbers.
- Step 2 and 3 of the pipeline (mechanism classes, antigen and biofilm layers) are not built.
- Issues #9, #10, #12 to #16, #19 and #26 remain open.
