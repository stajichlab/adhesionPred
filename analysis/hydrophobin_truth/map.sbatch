#!/bin/bash
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J hyd_map
#SBATCH -o /bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/hydrophobin_truth/map.%j.log
set -euo pipefail
module load ncbi-blast/2.14.0+
cd /bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/hydrophobin_truth
S=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat
F5=/bigdata/stajichlab/shared/projects/Fungi_5k/input
/usr/bin/python3.12 map_to_proteomes.py --truth truth_all.tsv --out . \
  --proteome Afum_Af293=$S/Afum_Af293_UniProt.faa --proteome Afum_A1163=$S/Afum_A1163.faa --proteome Afum_W72310=$S/Afum_W72310.faa \
  --proteome Scer_S288C=$S/Scer_S288C.faa --proteome Calb_SC5314=$S/Calb_SC5314.faa --proteome Cimm_RS=$S/Cimm_RS.faa \
  --proteome Bder_ER3=$S/Bder_ER3.faa --proteome Cneo_H99=$S/Cneo_H99.faa \
  --proteome Bbas_ARSEF2860=$F5/Beauveria_bassiana_ARSEF_2860.proteins.fa --proteome Fvel_6-3=$F5/Flammulina_velutipes_6-3.proteins.fa \
  --proteome Fful_Race5=$F5/Fulvia_fulva_Race5_Kim.proteins.fa --proteome Fgra_PH-1=$F5/Fusarium_graminearum_PH-1.proteins.fa \
  --proteome Pexp_MD-8=$F5/Penicillium_expansum_MD-8.proteins.fa --proteome Post_PC9=$F5/Pleurotus_ostreatus_PC9.proteins.fa \
  --proteome Tasp_FT101=$F5/Trichoderma_asperellum_FT101.proteins.fa --proteome Tvir_Gv29-8=$F5/Trichoderma_virens_Gv29-8.proteins.fa
