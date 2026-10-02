# cys_candidates: short, Cys-rich, secreted protein candidates (Coccidioides)

This analysis lists short proteins that have a signal peptide and many Cys in a short window.
The aim is to find proteins that might behave like PRA3 (class "2b-ii" in
`docs/TOOL-ARCHITECTURE.md`). It is a read-only filter. It has no model, no training and no
accuracy claim. Python 3.12, standard library only.

## Files

| File | Purpose |
| --- | --- |
| `00_make_manifests.sh` | Writes `proteomes.tsv` (3 columns, for 01) and `manifest.tsv` (4 columns, for `cys_candidates.py`). |
| `cys_candidates.py` | Computes features and a tier for each protein. Writes TSV.gz, `candidates.tsv.gz`, `summary.tsv`, `run.json`. |
| `01_known_family_hmm.sh` | `hmmfetch` six Pfam models and run `hmmsearch --cut_ga` on each proteome. Writes `pfam_provenance.json`. |
| `02_pra3_orthologs.py` | One-off check. mmseqs search of UniProt Q2TVJ9 against each proteome, joined to the tiers. |
| `03_calibrate.py` | Measures the controls and the SP-called distribution. Source of the threshold numbers below. |
| `04_no_sp_profile.py` | Profiles the `cys_rich_no_sp` tier (Pfam hits, first residue, SP probability). |
| `tests/cys_candidates/` | pytest. Run `/usr/bin/python3.12 -m pytest tests/cys_candidates`. |

## How to run

```bash
P=/bigdata/stajichlab/jstajich/projects/adhesionPred
OUT=$P/_workdir/cys_candidates
Q=$P/_workdir/class2b_structure/fasta/PRA3_Q2TVJ9.fasta
# 1. manifests: proteomes.tsv (name, fasta, signalp) and manifest.tsv (+ domtbl path)
OUT_DIR=$OUT bash analysis/cys_candidates/00_make_manifests.sh
# 2. Pfam search (login node, about 25 s). PFAM_HMM defaults to the central database link.
OUT_DIR=$OUT bash -l analysis/cys_candidates/01_known_family_hmm.sh
# 3. tiers
/usr/bin/python3.12 analysis/cys_candidates/cys_candidates.py --manifest $OUT/manifest.tsv \
    --out-dir $OUT --pfam-provenance $OUT/pfam_provenance.json --repo .
# 4. one-off PRA3 ortholog check (mmseqs, about 1 minute)
/usr/bin/python3.12 analysis/cys_candidates/02_pra3_orthologs.py --manifest $OUT/manifest.tsv \
    --table-dir $OUT --query $Q --out $OUT/pra3_orthologs.tsv
# 5. calibration tables
/usr/bin/python3.12 analysis/cys_candidates/03_calibrate.py --manifest $OUT/manifest.tsv \
    --control-fasta-dir $P/_workdir/class2b_structure/fasta --pra3-id CIMG_02492-t26_1-p1
# 6. what the no_sp tier holds (see the section on that tier)
/usr/bin/python3.12 analysis/cys_candidates/04_no_sp_profile.py --manifest $OUT/manifest.tsv \
    --table-dir $OUT --write-fasta $OUT/no_sp.faa
bash -l -c "module load hmmer/3.4; hmmsearch --cut_ga --noali -o /dev/null \
    --domtblout $OUT/no_sp.domtbl <Pfam-A.hmm> $OUT/no_sp.faa"      # about 3 minutes
/usr/bin/python3.12 analysis/cys_candidates/04_no_sp_profile.py --manifest $OUT/manifest.tsv \
    --table-dir $OUT --domtbl $OUT/no_sp.domtbl
```

Single-proteome mode: `--name N --fasta F --signalp S [--domtbl D --pfam-provenance J]`. FASTA,
SignalP and domtbl files may be plain or gzip. A domtbl needs `--pfam-provenance`. Any input
problem prints `STOP: <reason>` to stderr and exits 2, including argument errors, a missing
manifest and an unwritable `--out-dir`. No output file is kept in that case. The FASTA header
(the full line after `>`) must match the SignalP ID column exactly. The domtbl target is the first
token of the header. If `--out-dir` holds a `*.tsv.gz` from an earlier run that this run does not
write (for example another proteome name), the script stops and names the file. Remove it or use
a new directory.

Sequence rule: lowercase letters are uppercased. One trailing `*` is removed. An internal `*`,
any other non-letter, or an empty sequence is a STOP that names the protein id.

Pfam source: `01_known_family_hmm.sh` uses the central database link
`/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm`, resolves it, and writes the
release directory name, the model accessions with versions and the sha256 of the fetched
`known_families.hmm` to `pfam_provenance.json`. `cys_candidates.py` copies that file into
`run.json` (field `pfam`). The run in this README used release `2026-01-27-Pfam38.2`. Models:
PF05730.17 (CFEM), PF04681.18 (Bys1), PF01185.24, PF06766.17 and PF28987.1 (DewD) as
hydrophobin, PF28404.1 (ARB_05178) as `pra3_like_family`. If the database cannot be read, 01
prints STOP and exits 2 and leaves no output.

## Columns (one row per protein)

`proteome`, `protein_id` (first token), `header`, `length`, `sp_call` (SP or OTHER), `sp_prob`
(SP(Sec/SPI) probability), `cs_end` (last residue of the signal peptide, from `CS pos: N-M`),
`mature_length` (length after SP removal; full length if no SP), `cys_count`, `cys_frac`
(both on the mature sequence), `max_cys_window` and `window_start` (max Cys in any window of W
residues and its first 1-based start in the mature sequence), `cc_pairs` (positions i with C at i
and i+1; CCC counts 2), `cxc` (positions i with C at i and i+2; x is any residue),
`pest_frac` (fraction of P, S, T, E), family flags `cfem`, `bys1`, `hydrophobin`,
`pra3_like_family` with `<family>_evalue` (best domain i-Evalue), and `tier`. A flag is `NA` when
no search was run for that family (option `--missing-models`; not used in the run below).

## Tiers and thresholds

| Tier | Rule |
| --- | --- |
| `cys_rich_sp_unassigned` | SP called, mature length <= L, max Cys in a W window >= K, no cfem/bys1/hydrophobin hit. `pra3_like_family` hits stay here and are flagged. |
| `cys_rich_sp_known_family` | Same, with a cfem, bys1 or hydrophobin hit. Reported, not a candidate. |
| `cys_rich_no_sp` | No SP call. Full length <= L and max Cys in a W window >= K on the full sequence. Family flags are reported but not used. This is NOT a candidate list; see below. |
| `other` | All else. |

Defaults: `--max-mature-len 300` (L), `--min-cys-window 8` (K), `--window 60` (W). No tier depends
on a protein ID (a test checks this). A PF28404 hit does not move a protein out of
`cys_rich_sp_unassigned`. A hydrophobin hit (PF01185, PF06766 or PF28987) does.

### How the defaults were chosen (measured with `03_calibrate.py`)

Controls, features on the full UniProt sequence (no SignalP run on them), W = 60:

| Control | Length | Cys | Max Cys in 60 aa | C-C | C-x-C |
| --- | --- | --- | --- | --- | --- |
| PRA3 Q2TVJ9 (UniProt) | 153 | 8 | 8 | 2 | 0 |
| PRA3 annotated, CIMG_02492, after SP removal | 201 (full 220) | 12 | 11 | 2 | 0 |
| Ag2/PRA Q6QJA6 (CFEM) | 194 | 8 | 8 | 0 | 1 |
| PRA2 Q6K1L8 (CFEM) | 124 | 9 | 8 | 0 | 1 |
| RodA P41746 (hydrophobin) | 159 | 8 | 5 | 2 | 0 |
| CalA Q4WXJ1 (Bys1) | 177 | 5 | 5 | 0 | 0 |
| SOWgp Q8NK60 | 328 | 21 | 6 | 0 | 0 |
| CTS1 Q1E3R8 | 427 | 2 | 1 | 0 | 0 |

SP-called proteins per proteome: 436 to 686 (460 in C. immitis RS). Median mature length is 341
to 358. The 95th percentile of max Cys in a 60-aa window is 8 in all seven proteomes. The median
is 3.

- K = 8: this is the highest K at which PRA3 (UniProt and annotated), Ag2/PRA and PRA2 all pass
  (their max values are 8, 8, 8, 11). K = 9 would drop three of these four forms. This is a choice
  made from the controls. The positive set is small (see limits below).
- W = 60: with W = 40, Ag2/PRA and PRA2 reach only 6, so K = 8 would drop them. With W = 60 the
  three controls without a Cys cluster (RodA 5, CalA 5, SOWgp 6) stay below 8. With W = 100,
  SOWgp reaches 10 and RodA 8.
  Cost of W = 60: more proteins pass. At L = 300 and K = 8, summed over the seven proteomes,
  150 SP-called proteins pass at W = 60 and 53 at W = 40. W = 40 with K = 8 still keeps all four
  PRA3-family members of C. immitis RS (8, 8, 8, 8 Cys in 40 aa; measured), but it drops the CFEM
  controls (6). W = 60 was chosen so that the CFEM controls reach K = 8 and are reported as a
  known family.
- L = 300: all small-protein controls have mature length <= 220. SOWgp (328) and CTS1 (427)
  fall outside. L = 300 is a round value above the controls, not a fitted value.

Cost at each tried threshold. SP-called proteins with mature length <= L and max Cys >= K
(W = 60). Columns are the seven proteomes in the order RS, CiB10637, CiB10992, VFC140, Cpos1038,
Cpos3700, Silveira.

| L | K >= 7 | K >= 8 | K >= 9 | K >= 10 |
| --- | --- | --- | --- | --- |
| 200 | 17, 11, 12, 14, 10, 15, 14 | 16, 11, 12, 14, 10, 15, 13 | 6, 4, 5, 5, 4, 7, 5 | 3, 2, 3, 3, 3, 4, 3 |
| 250 | 26, 22, 22, 22, 21, 28, 22 | 23, 18, 19, 19, 16, 24, 19 | 9, 7, 8, 8, 7, 11, 8 | 6, 5, 6, 6, 6, 8, 6 |
| 300 | 30, 23, 23, 26, 23, 30, 27 | 26, 19, 20, 21, 17, 26, 21 | 12, 8, 9, 10, 8, 13, 10 | 8, 6, 7, 8, 7, 9, 7 |
| 400 | 34, 27, 27, 30, 27, 40, 31 | 28, 20, 21, 23, 19, 30, 23 | 13, 8, 9, 11, 8, 13, 10 | 8, 6, 7, 9, 7, 9, 7 |

At the defaults (L 300, K 8) 17 to 26 SP-called proteins pass per proteome, 3.7 to 5.7 percent
of the SP-called set. The full grid (K 6 to 12) is in `03_calibrate.py` output.

## Result (default thresholds, 7 proteomes)

Counts of proteins per tier (`summary.tsv`):

| Proteome | Proteins | SP | sp_unassigned | sp_known_family | no_sp | other |
| --- | --- | --- | --- | --- | --- | --- |
| CimmitisRS_FungiDB | 9910 | 460 | 21 | 5 | 31 | 9853 |
| Coccidioides_immitis_CiB10637 | 9256 | 447 | 15 | 4 | 18 | 9219 |
| Coccidioides_immitis_CiB10992 | 9414 | 441 | 16 | 4 | 19 | 9375 |
| Coccidioides_immitis_VFC140 | 10850 | 492 | 18 | 3 | 25 | 10804 |
| Coccidioides_posadasii_Cpos1038 | 9694 | 460 | 13 | 4 | 17 | 9660 |
| Coccidioides_posadasii_Cpos3700 | 14686 | 686 | 22 | 4 | 23 | 14637 |
| CposadasiiSilveira2022_FungiDB | 8516 | 436 | 17 | 4 | 21 | 8474 |

`candidates.tsv.gz` has 304 rows (all tiers except `other`). The DewD model (PF28987) moves one
protein in RS, CIMG_03068 (94 aa, 68 aa mature, 8 Cys, 0.1176 Cys fraction), from
`cys_rich_sp_unassigned` to `cys_rich_sp_known_family`. A first run with Pfam 38.0 (no PF28404,
no PF28987) gave 22 and 4 for RS. Other proteomes did not change.

PRA3 and its relatives (`02_pra3_orthologs.py`, mmseqs easy-search of Q2TVJ9, ortholog = fident
>= 0.90 and alignment >= 100 aa; weak hit = other hit with e-value <= 1e-3):

- The PF28404 model (Pfam 38.2) flags exactly 30 proteins: 4 per proteome, 6 in Cpos3700. These
  are the same 30 proteins as the mmseqs ortholog and weak-hit list below (set comparison done
  once after the run). The two methods are independent. All 30 are `cys_rich_sp_unassigned`.
- Every proteome has an ortholog (fident 0.973 to 1.000 over 153 aa; Cpos3700 has two copies,
  CPOS3700_005120 and CPOS3700_009276). All are length 220, SP called, mature length 201, 12 Cys,
  max Cys in window 11, tier `cys_rich_sp_unassigned`.
- Three further paralog-like hits per proteome (four in Cpos3700, which has two copies of the
  third one) (for RS: CIMG_07303, CIMG_05560, CIMG_07843;
  fident 0.714 over 49 aa, 0.462 over 157 aa, 0.373 over 134 aa) have length 214 to 231, 11 to 12
  Cys, and are also `cys_rich_sp_unassigned`.
- Of the CFEM controls, CIMG_09696 is Ag2/PRA (mmseqs fident 1.000 over 194 aa) and CIMG_09560 is
  PRA2 (0.983 over 124 aa). Both are `cys_rich_sp_known_family`.

## PRA3 truncation observation

UniProt Q2TVJ9 (TrEMBL, `PE=2`) is 153 aa with 8 Cys and starts `MCYNPGRG`. SignalP was not run
on it. The `cys_rich_no_sp` tier assumes no SP. The test that puts it in that tier uses an
invented SignalP row (`OTHER`), not a SignalP result. The annotated ortholog CIMG_02492 is 220
aa. SignalP 6 calls an SP (SP probability 0.9996, cleavage after residue 19; mature 201 aa, 12
Cys).

The UniProt sequence matches the annotated proteins from 0-based offset 67 (the UniProt Met
aligns to residue 68), ungapped, and ends at the same C-terminal residue. Identical residues out
of 153, measured for each ortholog: QVM06151.1 (Silveira) 153, CPOS3700_005120 153, VFC140_000315
153, CPOS1038_002252 152, CPOS3700_009276 149, CIMG_02492 149, CIB10637_000324 149. The 67
residues before the offset contain 4 Cys (none in the 19-residue SP) and account for the
difference between 8 and 12 Cys.

Hypothesis, not tested here: the UniProt entry starts at a downstream Met and lacks about 67
N-terminal residues. Because Q2TVJ9 is `PE=2`, this analysis does not show which start is
correct. The earlier structure survey ("38 aa, 7 Cys" confident core, no fold) used the 153-aa
model. This observation does not overturn that result. Only a new fold of the 201-aa mature form
can test it.

## What the cys_rich_no_sp tier holds

Method: the 154 `cys_rich_no_sp` proteins of the seven proteomes were searched with
`hmmsearch --cut_ga` against the full Pfam-A of release 2026-01-27-Pfam38.2 (30,134 families),
using `04_no_sp_profile.py`. The zinc-finger / Fe-S / RING-type count is a keyword rule on the
Pfam name (`ZN_FE_RULE` in the script), not a curated annotation.

- 90 of 154 have at least one Pfam hit; 64 have none.
- 78 of 154 have a domain whose name matches zinc finger, Fe-S (Fer4), RING-type, BBOX, MYND,
  IBR or similar (for example zf-CCHC 16 proteins, zf-ANAPC11 12, zf-rbx1 12, Fer4 8, zf-C3HC4 8,
  zf-CCCH 7). 12 more have a Pfam hit of another kind (CYSTM 3, PLAC8 3, DUF3716 3, DUF7871,
  RhoGAP, Rcd1).
- All 154 start with Met. The maximum SP probability is 0.2932.

Reading: these are mostly intracellular Cys-rich proteins. This tier is not a candidate list. It
shows no sign of N-terminal truncation by two measures (first residue M, low SP probability). It
does not prove that no truncated secreted protein is in it.

## Next step: ORF input

The script is meant to be reused on ORF-translated sequences. Before that:

- Run SignalP on every ORF first. The script stops when a protein has no SignalP row.
- Remove internal stop symbols (`*`). The script stops on them.
- Expect the `cys_rich_no_sp` tier to grow, because ORFs include many non-genes.
- Run time and memory, measured for one proteome (C. immitis RS, 9,910 proteins, with domtbl):
  4.1 s wall time and 67 MB peak memory (`/usr/bin/time -v`). A run of 500,000 ORFs would need
  about 50 times that (estimate by linear scaling, not measured). All proteins are held in memory.

## What this cannot show

- No structure and no function. Cys-rich does not mean PRA3-like.
- The positive set is one protein family (PRA3), so there is no hit rate, no precision and no
  recall. The thresholds were set from 7 control sequences.
- SignalP gave an SP call to 4.5 to 5.1 percent of proteins here (436 of 8516 to 686 of 14686).
  No truth data show whether it under-calls in Coccidioides.
- Pfam releases: the first run of this analysis used the shared `funannotate_db` copy, which was
  Pfam 38.0 when run. I did not check other releases then. It has no PF28404 or PF28987. The run
  reported here uses the central release 2026-01-27-Pfam38.2. The funannotate_db copy now says
  38.2 too, but this analysis does not rely on it.
- Controls other than PRA3 had no SignalP run. Their features are on the full sequence.
- Max Cys in a window ignores the spacing pattern of the Cys.
