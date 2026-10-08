#!/bin/bash
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J csh_hyd_disc
set -euo pipefail
cd /rhome/jstajich/projects/adhesionPred/_workdir/sorting_hat/calibration/hydrophobin_discovery
module load hmmer/3.4
PFAM=/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm
printf '%s\n' Hydrophobin Hydrophobin_2 Eas DewD Hyd1F Hydrophobin_D Hydrophobin_like HsbA Cerato-platanin Flocculin Flocculin_t3 Hyr1 PIR1-like_C ALS_M Candida_ALS Candida_ALS_N Hyphal_reg_CWP PIR > keys.txt
hmmfetch -f $PFAM keys.txt > models.hmm
grep -c '^NAME' models.hmm
for p in Afum_Af293_UniProt Afum_A1163 Afum_W72310 Scer_S288C Calb_SC5314 Cimm_RS Bder_ER3 Cneo_H99; do
  hmmsearch --cut_ga --cpu 8 --noali -o /dev/null --domtblout $p.domtbl models.hmm /rhome/jstajich/projects/adhesionPred/_workdir/sorting_hat/$p.faa
done
