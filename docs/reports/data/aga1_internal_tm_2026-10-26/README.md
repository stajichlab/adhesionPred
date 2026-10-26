# GPI calls and internal transmembrane helices (2026-10-26)

**Why.** PredGPI can read the last transmembrane helix of a membrane protein as a GPI signal (seen in the hand-check sample: a ZIP-type transporter, WSC-type sensors). A GPI wall protein has no mature
helix except its own C-terminal signal. **Rule.** A helix that starts after residue 35 (past the signal peptide) and ends more than 30 residues before the C-terminus is *internal*. The `tm` module (now version 2) writes the count as
`n_tm_internal`; the new call `gpi_no_internal_tm_protein` is "GPI called and no internal helix"; `gpi_long_disorder_protein` (and so the broad and narrow flags) uses the same condition. `gpi_anchored_protein` is unchanged (the raw predictor call).
Where the window comes from: on the 1,099 GPI-called proteins of the ten proteomes the end of each mature helix lies 0 to 9 residues from the C-terminus in 231 cases and rarely 10 to 39 (25 cases), so 30 sits in a gap; the cut is a parameter of the module (`terminal_helix_window`).

**Effect on the ten proteomes** (engine runs `<species>/out_tm_internal`, config sha256 `a9dc86ce...`; counts of proteins called):

| species | gpi_anchored | gpi_no_internal_tm | gpi_long_disorder | broad flag | narrow flag |
|---|---|---|---|---|---|
| Scer_S288C | 82 | 77 | 10 (was 11) | 10 (was 11) | 5 |
| Calb_SC5314 | 135 | 123 | 32 | 25 | 6 |
| Cgla_CBS138 | 117 | 112 | 27 | 25 | 4 |
| Kmar_FIM1 | 80 | 70 | 10 | 10 | 2 |
| Caur_UniProt | 108 | 99 | 16 | 14 | 5 |
| Cpar_CDC317 | 127 | 116 | 28 (was 29) | 27 (was 28) | 8 (was 9) |
| Afum_Af293 | 150 | 137 | 7 | 7 | 3 |
| Cimm_RS | 91 | 83 | 8 (was 9) | 7 (was 8) | 2 |
| Cneo_H99 | 89 | 81 | 2 (was 3) | 1 (was 2) | 0 |
| Rtor_NP11 | 120 | 111 | 6 (was 7) | 2 | 1 |

90 of the 1,099 GPI calls (8%) have an internal helix, but only six change `gpi_long_disorder_protein` (most of them were not in a long disordered run anyway) and four change the broad flag:
P41809 (S. cerevisiae, one helix near residue 1,500), G8BH67 (C. parapsilosis, 8 helices; in the hand-check sample), J3KGV0 (C. immitis, 5 helices), T2BQH2 (C. neoformans, 5 helices; in the hand-check sample).
Rtor M7WLZ4 leaves `gpi_long_disorder_protein` but was not in the broad flag for another reason.

**Status files.** Both call status files of the acceptance record were re-measured: sensitivity unchanged (S. cerevisiae 0.75 with 8 positives, C. albicans 0.90 with 10), no tier A protein is lost; S. cerevisiae specificity 0.981 to 0.985.

**How.** `python3 -I analysis/sorting_hat_run/aga1_domain_modules.py --tm-only --out-name out_tm_internal --runs-root <runs> --taxdump nodes.dmp --taxa analysis/sorting_hat_run/aga1_acceptance_taxa.tsv --results <results> ...` then `cellsurface_sorting_hat_calibrate truth --call-status` as in the acceptance record README (with the note text naming the extra condition).
Runs without a `tm` table now give `not_assessable` for `gpi_long_disorder_protein` and the flags when GPI is called; the `tm` tables of older runs need the new column (the module writes it).
