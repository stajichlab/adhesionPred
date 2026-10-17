# Family-matched non-allergen negatives (all-kingdom allergens)

- Date: 2026-10-10; UniProt release 2026_03 (relnotes.txt, 02-Sep-2026); Swiss-Prot only (reviewed:true)
- Seed 20261010; cap 10 per family per group; MMseqs2 easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5

## Positives

| Group | Positives | With Pfam | Without Pfam | Pfam families |
|---|---|---|---|---|
| plant | 913 | 868 | 45 | 81 |
| animal | 766 | 718 | 48 | 130 |
| total | 1679 | 1586 | 93 | 178 |

Families used by both groups: 33 (sampled separately per group).
Positives with more than one AllergenOnline type: 26. Group rule: plant if any type is in ['Aero Plant', 'Food Plant', 'Gliadin']; otherwise animal.
Animal-group positives whose types are only Bacteria/Protozoan/Vaccine: 10 (sampled from Metazoa per the rule).

## Negatives

| Group | Negatives | Families with >=1 negative |
|---|---|---|
| plant | 654 | 71 |
| animal | 893 | 102 |
| total | 1547 | 154 |

Families with negatives in both groups: 19.

## Candidate pools

| Group | Family-group pairs | Pairs with empty pool | Pairs with 10 sampled | Pairs with <10 sampled |
|---|---|---|---|---|
| plant | 81 | 10 | 58 | 23 |
| animal | 130 | 28 | 81 | 49 |

Candidate hits removed by exclusion (summed over pairs, a protein can count in several pairs): in positives 109, in kingdom negatives 215.

## Clusters (MMseqs2, positives and new negatives together)

| Cluster type | Count |
|---|---|
| positives only | 131 |
| negatives only | 242 |
| mixed | 123 |
| total | 496 |

## Per-family counts

| Group | Family | Family name | Pool after exclusion | Sampled |
|---|---|---|---|---|
| animal | PF00012 | Hsp70 protein | 104 | 10 |
| animal | PF00014 | Kunitz/Bovine pancreatic trypsin inhibitor domain | 366 | 10 |
| animal | PF00026 | Eukaryotic aspartyl protease | 68 | 10 |
| animal | PF00031 | Cystatin domain | 94 | 10 |
| animal | PF00036 | EF hand domain | 82 | 10 |
| animal | PF00042 | Globin | 980 | 10 |
| animal | PF00043 | Glutathione S-transferase, C-terminal domain | 91 | 10 |
| animal | PF00044 | Glyceraldehyde 3-phosphate dehydrogenase, NAD binding domain | 42 | 10 |
| plant | PF00044 | Glyceraldehyde 3-phosphate dehydrogenase, NAD binding domain | 43 | 10 |
| animal | PF00050 | Kazal-type serine protease inhibitor domain | 228 | 10 |
| animal | PF00061 | Lipocalin / cytosolic fatty-acid binding protein family | 180 | 10 |
| plant | PF00061 | Lipocalin / cytosolic fatty-acid binding protein family | 0 | 0 |
| animal | PF00062 | C-type lysozyme/alpha-lactalbumin family | 136 | 10 |
| animal | PF00079 | Serpin (serine protease inhibitor) | 191 | 10 |
| plant | PF00079 | Serpin (serine protease inhibitor) | 27 | 10 |
| plant | PF00080 | Copper/zinc superoxide dismutase (SODC) | 40 | 10 |
| animal | PF00081 | Iron/manganese superoxide dismutases, alpha-hairpin domain | 27 | 10 |
| animal | PF00082 | Subtilase family | 46 | 10 |
| plant | PF00082 | Subtilase family | 61 | 10 |
| animal | PF00085 | Thioredoxin | 121 | 10 |
| plant | PF00085 | Thioredoxin | 136 | 10 |
| animal | PF00089 | Trypsin | 627 | 10 |
| animal | PF00091 | Tubulin/FtsZ family, GTPase domain | 109 | 10 |
| animal | PF00093 | von Willebrand factor type C domain | 34 | 10 |
| animal | PF00094 | von Willebrand factor type D domain | 4 | 4 |
| animal | PF00112 | Papain family cysteine protease | 52 | 10 |
| plant | PF00112 | Papain family cysteine protease | 69 | 10 |
| animal | PF00113 | Enolase, C-terminal TIM barrel domain | 32 | 10 |
| animal | PF00121 | Triosephosphate isomerase | 27 | 10 |
| plant | PF00121 | Triosephosphate isomerase | 13 | 10 |
| animal | PF00128 | Alpha amylase, catalytic domain | 58 | 10 |
| animal | PF00149 | Calcineurin-like phosphoesterase | 151 | 10 |
| plant | PF00150 |  | 0 | 0 |
| animal | PF00151 | Lipase | 52 | 10 |
| animal | PF00160 | Cyclophilin type peptidyl-prolyl cis-trans isomerase/CLD | 82 | 10 |
| plant | PF00160 | Cyclophilin type peptidyl-prolyl cis-trans isomerase/CLD | 36 | 10 |
| animal | PF00168 | C2 domain | 183 | 10 |
| animal | PF00171 | Aldehyde dehydrogenase family | 89 | 10 |
| animal | PF00182 | Chitinase class I | 0 | 0 |
| plant | PF00182 | Chitinase class I | 70 | 10 |
| animal | PF00187 | Chitin recognition protein | 0 | 0 |
| plant | PF00187 | Chitin recognition protein | 68 | 10 |
| animal | PF00188 | Cysteine-rich secretory protein family | 108 | 10 |
| plant | PF00188 | Cysteine-rich secretory protein family | 14 | 10 |
| animal | PF00190 | Cupin | 0 | 0 |
| plant | PF00190 | Cupin | 145 | 10 |
| plant | PF00197 | Trypsin and protease inhibitor | 68 | 10 |
| animal | PF00210 | Ferritin-like domain | 50 | 10 |
| animal | PF00217 | ATP:guanido phosphotransferase, C-terminal catalytic domain | 49 | 10 |
| animal | PF00234 | Protease inhibitor/seed storage/LTP family | 0 | 0 |
| plant | PF00234 | Protease inhibitor/seed storage/LTP family | 132 | 10 |
| animal | PF00235 | Profilin | 25 | 10 |
| plant | PF00235 | Profilin | 15 | 10 |
| animal | PF00241 | Cofilin/tropomyosin-type actin-binding protein | 50 | 10 |
| plant | PF00251 | Glycosyl hydrolases family 32 N-terminal domain | 47 | 10 |
| animal | PF00261 | Tropomyosin | 52 | 10 |
| animal | PF00273 | Serum albumin family | 30 | 10 |
| animal | PF00274 | Fructose-bisphosphate aldolase class-I | 28 | 10 |
| animal | PF00280 | Potato inhibitor I family | 3 | 3 |
| plant | PF00280 | Potato inhibitor I family | 21 | 10 |
| plant | PF00295 | Glycosyl hydrolases family 28 | 24 | 10 |
| plant | PF00304 | Gamma-thionin family | 87 | 10 |
| animal | PF00307 | Calponin homology (CH) domain | 77 | 10 |
| animal | PF00314 | Thaumatin family | 0 | 0 |
| plant | PF00314 | Thaumatin family | 31 | 10 |
| animal | PF00326 | Prolyl oligopeptidase family | 23 | 10 |
| animal | PF00328 | Histidine phosphatase superfamily (branch 2) | 28 | 10 |
| animal | PF00332 | Glycosyl hydrolases family 17 | 0 | 0 |
| plant | PF00332 | Glycosyl hydrolases family 17 | 55 | 10 |
| animal | PF00363 | Casein | 27 | 10 |
| animal | PF00372 | Hemocyanin, copper containing domain | 58 | 10 |
| animal | PF00405 | Transferrin | 35 | 10 |
| animal | PF00407 | Pathogenesis-related protein Bet v 1 family | 0 | 0 |
| plant | PF00407 | Pathogenesis-related protein Bet v 1 family | 42 | 10 |
| animal | PF00435 |  | 0 | 0 |
| plant | PF00450 | Serine carboxypeptidase | 68 | 10 |
| animal | PF00477 |  | 0 | 0 |
| animal | PF00544 | Pectate lyase | 0 | 0 |
| plant | PF00544 | Pectate lyase | 28 | 10 |
| animal | PF00578 | AhpC/TSA family | 38 | 10 |
| plant | PF00578 | AhpC/TSA family | 23 | 10 |
| animal | PF00626 | Gelsolin repeat | 36 | 10 |
| plant | PF00626 | Gelsolin repeat | 7 | 7 |
| animal | PF00630 | Filamin/ABP280 repeat | 15 | 10 |
| animal | PF00657 | GDSL-like Lipase/Acylhydrolase | 6 | 6 |
| animal | PF00704 | Glycosyl hydrolases family 18 | 49 | 10 |
| plant | PF00704 | Glycosyl hydrolases family 18 | 21 | 10 |
| animal | PF00719 | Inorganic pyrophosphatase | 7 | 7 |
| plant | PF00732 | GMC oxidoreductase | 16 | 10 |
| plant | PF00903 | Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily | 12 | 10 |
| animal | PF00930 | Dipeptidyl peptidase IV (DPP IV) N-terminal region | 5 | 5 |
| animal | PF00967 |  | 0 | 0 |
| animal | PF00997 | Kappa casein | 40 | 10 |
| animal | PF01023 | S-100/ICaBP type calcium binding domain | 60 | 10 |
| animal | PF01042 | Endoribonuclease L-PSP | 4 | 4 |
| plant | PF01095 | Pectinesterase | 82 | 10 |
| animal | PF01099 | Uteroglobin family | 15 | 10 |
| animal | PF01123 |  | 0 | 0 |
| plant | PF01190 | Pollen protein Ole e 1 like | 6 | 6 |
| animal | PF01221 | Dynein light chain type 1 | 11 | 10 |
| animal | PF01273 | LBP / BPI / CETP family, N-terminal domain | 38 | 10 |
| animal | PF01277 | Oleosin | 0 | 0 |
| plant | PF01277 | Oleosin | 33 | 10 |
| animal | PF01347 |  | 0 | 0 |
| plant | PF01357 | Expansin C-terminal domain | 91 | 10 |
| animal | PF01372 | Melittin | 5 | 5 |
| plant | PF01373 | Glycosyl hydrolase family 14 | 20 | 10 |
| animal | PF01391 | Collagen triple helix repeat (20 copies) | 149 | 10 |
| animal | PF01395 | PBP/GOBP family | 82 | 10 |
| animal | PF01410 | Fibrillar collagen C-terminal domain | 2 | 2 |
| plant | PF01439 | Metallothionein | 47 | 10 |
| animal | PF01459 | Eukaryotic porin | 31 | 10 |
| plant | PF01565 | FAD binding domain | 81 | 10 |
| animal | PF01576 | Myosin tail | 14 | 10 |
| animal | PF01607 | Chitin binding Peritrophin-A domain | 20 | 10 |
| plant | PF01620 |  | 0 | 0 |
| animal | PF01630 | Hyaluronidase | 49 | 10 |
| animal | PF01717 | Cobalamin-independent synthase, Catalytic domain | 0 | 0 |
| plant | PF01717 | Cobalamin-independent synthase, Catalytic domain | 9 | 9 |
| animal | PF01734 | Patatin-like phospholipase | 13 | 10 |
| plant | PF01734 | Patatin-like phospholipase | 46 | 10 |
| animal | PF01826 | Trypsin Inhibitor like cysteine rich domain | 21 | 10 |
| animal | PF02098 | Tick histamine binding protein | 5 | 5 |
| plant | PF02209 |  | 0 | 0 |
| animal | PF02221 | ML domain | 11 | 10 |
| animal | PF02271 | Ubiquinol-cytochrome C reductase complex 14kD subunit | 6 | 6 |
| plant | PF02298 | Plastocyanin-like domain | 42 | 10 |
| animal | PF02520 | SXP/RAL-2 family protein Ani s 5-like, metal-binding domain | 1 | 1 |
| animal | PF02704 | Gibberellin regulated protein | 0 | 0 |
| plant | PF02704 | Gibberellin regulated protein | 17 | 10 |
| animal | PF02777 | Iron/manganese superoxide dismutases, C-terminal domain | 16 | 10 |
| animal | PF02798 | Glutathione S-transferase, N-terminal domain | 121 | 10 |
| animal | PF02800 | Glyceraldehyde 3-phosphate dehydrogenase, C-terminal domain | 34 | 10 |
| plant | PF02800 | Glyceraldehyde 3-phosphate dehydrogenase, C-terminal domain | 34 | 10 |
| animal | PF02806 | Alpha amylase, C-terminal all-beta domain | 40 | 10 |
| animal | PF02807 | ATP:guanido phosphotransferase, N-terminal domain | 39 | 10 |
| animal | PF02872 | 5'-nucleotidase, C-terminal domain | 11 | 10 |
| animal | PF02876 |  | 0 | 0 |
| animal | PF02987 |  | 0 | 0 |
| animal | PF03022 | Major royal jelly protein | 29 | 10 |
| plant | PF03157 |  | 0 | 0 |
| plant | PF03227 | Gamma interferon inducible lysosomal thiol reductase (GILT) | 1 | 1 |
| animal | PF03260 | Lepidopteran low molecular weight (30 kD) lipoprotein | 6 | 6 |
| plant | PF03330 | Lytic transglycolase | 82 | 10 |
| animal | PF03392 | Insect pheromone-binding family, A10/OS-D | 3 | 3 |
| animal | PF03722 | Hemocyanin, all-alpha domain | 48 | 10 |
| animal | PF03723 | Hemocyanin, ig-like domain | 38 | 10 |
| plant | PF03913 |  | 0 | 0 |
| animal | PF03952 | Enolase, N-terminal domain | 18 | 10 |
| animal | PF03953 | Tubulin C-terminal domain | 99 | 10 |
| animal | PF03973 | Triabin | 8 | 8 |
| plant | PF04043 | Plant invertase/pectin methylesterase inhibitor | 63 | 10 |
| animal | PF04702 | Vicilin N terminal region | 0 | 0 |
| plant | PF04702 | Vicilin N terminal region | 5 | 5 |
| plant | PF05129 | Transcription elongation factor Elf1 like | 2 | 2 |
| plant | PF05199 | GMC oxidoreductase | 6 | 6 |
| plant | PF05368 | NmrA-like family | 40 | 10 |
| animal | PF05755 |  | 0 | 0 |
| animal | PF05826 | Phospholipase A2 | 21 | 10 |
| animal | PF05922 | Peptidase inhibitor I9 | 11 | 10 |
| plant | PF05922 | Peptidase inhibitor I9 | 51 | 10 |
| animal | PF06757 | Insect allergen related repeat, nitrile-specifier detoxification | 1 | 1 |
| plant | PF07250 | Glyoxal oxidase N-terminus | 2 | 2 |
| animal | PF07966 | A1 Propeptide | 52 | 10 |
| plant | PF07983 | X8 domain | 17 | 10 |
| plant | PF08031 | Berberine and berberine like | 41 | 10 |
| animal | PF08212 | Lipocalin-like domain | 3 | 3 |
| plant | PF08244 | Glycosyl hydrolases family 32 C terminal | 37 | 10 |
| animal | PF08246 | Cathepsin propeptide inhibitor domain (I29) | 28 | 10 |
| plant | PF08246 | Cathepsin propeptide inhibitor domain (I29) | 41 | 10 |
| animal | PF08267 |  | 0 | 0 |
| plant | PF08267 |  | 0 | 0 |
| animal | PF08726 |  | 0 | 0 |
| plant | PF09118 |  | 0 | 0 |
| animal | PF09172 |  | 0 | 0 |
| animal | PF09175 | Vitellinogen, beta-sheet shell | 1 | 1 |
| animal | PF09252 | Allergen Fel d I-B chain | 2 | 2 |
| plant | PF09253 |  | 0 | 0 |
| animal | PF10417 | C-terminal domain of 1-Cys peroxiredoxin | 26 | 10 |
| plant | PF10417 | C-terminal domain of 1-Cys peroxiredoxin | 9 | 9 |
| animal | PF11642 |  | 0 | 0 |
| plant | PF11837 | Beta-fructofuranosidase, N-terminal domain | 6 | 6 |
| plant | PF12070 |  | 0 | 0 |
| plant | PF13016 | Cys-rich Gliadin N-terminal | 31 | 10 |
| animal | PF13202 | EF hand | 85 | 10 |
| animal | PF13405 | EF-hand domain | 26 | 10 |
| plant | PF13410 | Glutathione S-transferase, C-terminal domain | 25 | 10 |
| plant | PF13417 | Glutathione S-transferase, N-terminal domain | 10 | 10 |
| animal | PF13499 | EF-hand domain pair | 430 | 10 |
| plant | PF13499 | EF-hand domain pair | 228 | 10 |
| plant | PF13561 | Enoyl-(Acyl carrier protein) reductase | 67 | 10 |
| plant | PF13920 | Zinc finger, C3HC4 type (RING finger) | 36 | 10 |
| plant | PF14368 | Probable lipid transfer | 42 | 10 |
| animal | PF14497 | Glutathione S-transferase, C-terminal domain | 59 | 10 |
| plant | PF14541 | Xylanase inhibitor C-terminal | 27 | 10 |
| plant | PF14543 | Xylanase inhibitor N-terminal | 19 | 10 |
| plant | PF14547 | Hydrophobic seed protein | 10 | 10 |
| animal | PF14625 | Lustrin, cysteine-rich repeated domain | 1 | 1 |
| animal | PF16469 |  | 0 | 0 |
| plant | PF16845 | Aspartic acid proteinase inhibitor | 25 | 10 |
| animal | PF16984 |  | 0 | 0 |
| plant | PF17181 | Epidermal patterning factor proteins | 8 | 8 |
| animal | PF17499 | Ant venom peptides | 25 | 10 |
| plant | PF17766 | Fibronectin type-III domain | 40 | 10 |
| animal | PF22750 |  | 0 | 0 |
| plant | PF24068 | Tapetum determinant 1, C-terminal | 5 | 5 |
| plant | PF24300 | Kiwellin-1 | 8 | 8 |
| animal | PF24805 | EIF3I | 29 | 10 |
| animal | PF25001 | Aegyptin, C-terminal domain | 7 | 7 |
| animal | PF28923 |  | 0 | 0 |
| animal | PF29934 | Vitellogenin N-terminal domain | 3 | 3 |
