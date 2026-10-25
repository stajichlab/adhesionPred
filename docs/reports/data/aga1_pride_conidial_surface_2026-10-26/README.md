# Conidial-surface proteomics evidence for *A. fumigatus* (PRIDE PXD031199, PXD044190), 2026-10-26

Public mass-spectrometry results joined to the AGA1-class pipeline calls for the Af293 proteome. The raw files are not used and not
in the repository; only protein-level results are, as accession-keyed tables. Script: `analysis/aga1_class/pride_surface_evidence.py`.
Source paper: Pinzan et al. 2024, Nat Microbiol, doi 10.1038/s41564-024-01782-y (PMID 39191887). PRIDE REST API
(https://www.ebi.ac.uk/pride/ws/archive/v2/), files from ftp.pride.ebi.ac.uk, downloaded 2026-10-25/26 into a directory outside the repo.

## What the data are

* **PXD031199**: trypsin shaving of conidia of *A. fumigatus* Af293 and three relatives (*A. oerlinghausenensis*, *A. fischeri*, *A. lentulus*), Q Exactive HF,
  searched with Proteome Discoverer against FungiDB (Af293 ids `Afu6g07910-T-p1`). 24 protein lists (`*.prot.xml`): runs 1-3 Af293, 4-6 *A. oerlinghausenensis*,
  7-9 *A. fischeri*, 10-12 *A. lentulus* (from the reference database named in each file), each with a second file `Ns`. The paper describes resting and swollen
  conidia; **mapping the `s` suffix to "swollen" is an inference from the file names and is recorded here only as `run_suffix`**. Only the six Af293 runs map to our
  Af293 proteome (Afu6g07910 to AFUA_6G07910, UniProt ordered-locus names); the other species' ids are listed in the run table but not mapped.
* **PXD044190**: MaxQuant 1.5.3.8 search of *A. fumigatus* A1163 (AFUB_ ids) plus mouse, LTQ Orbitrap, samples `MO_*` (an infection/macrophage experiment, not shaving).
  Of 1,695 protein groups, 83 contain an AFUB_ id; they were mapped to Af293 by the best shared 10-mers (at least 80% of the A1163 protein's 10-mers, unique best).
  23 Af293 proteins resulted, all intracellular (ribosomal proteins, Hsp90, metabolic enzymes), so this dataset adds no surface evidence here.

## Files

| File | Content |
|---|---|
| `pxd031199_run_detections.tsv` | one row per (run, protein): species, `run_suffix`, replicate, distinct peptides, spectral count, coverage, gene product (10,362 rows, 24 runs) |
| `afum_pride_surface_evidence.tsv` | one row per protein of the Af293 UniProt proteome (9,647): our predictor calls and tier flags, UniProt name, the number of Af293 runs it was detected in (of 6; plain and `s` runs), peptides, spectra, `detected_any`, `detected_2plus_runs`, and the PXD044190 columns `mq_*` |
| `enrichment.tsv` | Fisher exact tests: detection rate inside and outside each protein set (signal peptide, SP + GPI, broad rule, narrow rule, models) |
| `pxd044190_afumigatus_groups.tsv` | the A1163 proteins of PXD044190 with their Af293 match, k-mer share and per-sample MaxQuant intensity |

## Reading it

* **Detection is not presence.** Shaving finds exposed tryptic peptides on conidia; a protein not seen may be low in abundance, in another stage, or protected
  (Ser/Thr-rich O-mannosylated stalks resist trypsin, which affects exactly the class this project looks for). Only 380 of 9,647 Af293 proteins were detected in any of the six
  runs (288 in two or more).
* The protein lists are the authors' search results (Proteome Discoverer), not a re-analysis; the scores in the files are not on a 0-1 scale and are not used.
* A few UniProt accessions share one locus tag (8 rows); they receive the same detection (`accessions_for_locus`).

## Inputs not in the repository (sha256)

```
0f0c4f98efbbbd88b27c2744aa047dbdc2740cef73fbd22c34fc79e7fe637aaa  PXD031199/1.prot.xml
1eb25180e781527fab4bf1d020280603b9dd58cee13f706a2cbe32f6dd20b130  PXD031199/1s.prot.xml
c36e7ee355587c0e29cf38a4a7f6c9122458fc36b7eef1e2fe73e94ad9e48bda  PXD031199/2.prot.xml
4ed0d978c9005d94ecbedb345c025c23c18f72579173fe133bbf1b9f6a4bb888  PXD031199/2s.prot.xml
24fc8dac1f483755540eac73c5061e1e6290b86519381295e35c0b673da8048b  PXD031199/3.prot.xml
ecb3ff4457f4a420ab042109fafc60cc349d26f74680a94ceb0b328c301ca890  PXD031199/3s.prot.xml
c460fd2e40f078c02f127841003c6b12da71be68b026243285fcd2ea838d6489  PXD031199/4.prot.xml
9a4c7aab9a713532c934bacb179fb64fe6557dfa676c990d7e2b1d1ecee937ed  PXD031199/4s.prot.xml
9252f9ce62886f0527d18872da53d67946f7b2bf02aeb81972744585f02ad8c1  PXD031199/5.prot.xml
f8083dfdf449b1ab001f873a18ef4a94d0cef93143d50f00af0b41dec3ec466b  PXD031199/5s.prot.xml
362f313ef0d238570306b49dff2da0cea7d4edc0dff3b3151107a5925c0242b0  PXD031199/6.prot.xml
02f0f7dae342719bb4744a8f31acf2268a9b55c8bf415f6cb0a129d52b4f222e  PXD031199/6s.prot.xml
8b331d857297d341b5287af7c81c434bd153c34648ccb9d604d0abab1d94149d  PXD031199/7.prot.xml
66ae3ffeed6c6f44ba293004a3540940926a934c1e8b9b430ff97cf521e9e696  PXD031199/7s.prot.xml
115573856d9023b344f385101b57f29026ddc614b4c4f3e064e5c2c79e09d70b  PXD031199/8.prot.xml
ffd95385a79dc5d5af142421b0a9032d3a6206766d53187b4bf8db686e1f1300  PXD031199/8s.prot.xml
9d7e2082702b1169bdd25585b5a39ad0d9abdfb0b173e1dbfba7cff4b586b584  PXD031199/9.prot.xml
52ea7c0487e81876a3dbb5da657b1855b558f3719626b05a29dcd316af5a83b4  PXD031199/9s.prot.xml
9a9fc2c58eb255e49047380a9d4be7cbee8de087f88cdb82fd212f5729243e20  PXD031199/10.prot.xml
14eb9718a5beefd5a678a77391719361cfe4768224d9a80e5edee851b8b9946f  PXD031199/10s.prot.xml
68dc75e4f2c20b150572904426c8eb630010f53bfb3fed339bca4befa053282e  PXD031199/11.prot.xml
69186110656d15f9bf2556bcf8fea9fcb9d0bafbd5397f68b40843b7bf47b489  PXD031199/11s.prot.xml
518dee302bf145bd2c55482330efbad3fd8212c302c4bbd2995b6e1831c574eb  PXD031199/12.prot.xml
8c5bc5dc0152c5ab6ee80f8ab13ef52cb566441e8f4c74c85445bcf53dcb09e3  PXD031199/12s.prot.xml
16d6ab5210be433f7b46a71f6b3566f47dff67f954da8ab9e7201566f5566917  PXD044190/txt.rar
335b96c54454f41eb044556965c0c0f09a6b4e5da4106dbc4e982c82df1ba0a6  PXD044190/txt/proteinGroups.txt
ed4a23dbb7609cf7d41e16eb08358e4cbc0242bd2ce3f33597ca4420556608ef  afum_af293_uniprot_names.tsv
a66c6512871ce0463df62e653991da67f04c4449f09e5d07ef446ce1491f059d  afum_a1163_uniprot_names.tsv
74024269e9c8cad592fdf8cfc09a868dfce3d2d588647fe0e8301f22ed17aa19  Afum_A1163.UP000001699.faa
```

The Af293 FASTA used for the k-mer match is `Afum_Af293.UP000002530.faa` (hash in `docs/reports/data/aga1_nonyeast_test_2026-10-26/README.md`).
UniProt names were fetched on 2026-10-26 (release 2026_03).

## sha256 of the committed files (a test checks every line)

```
648768367eb784d36b6c432e69b58848a0f250a02e236a1741978bc0cda92543  afum_pride_surface_evidence.tsv
5750b7edd0ac091e8296ff4d15340b2213237eb4241488bcb19a7d58c7618abe  enrichment.tsv
7adb46d25dcb4936bb014d825a31a28b9462e853e78d07d231c986d5948c8a84  pxd031199_run_detections.tsv
bd63105af3260fc9968ef85c108c4610cce2bc61cf086593ba13591b37164d01  pxd044190_afumigatus_groups.tsv
```
