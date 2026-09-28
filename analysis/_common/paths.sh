# Centralized site paths for adhesionPred shell/sbatch scripts. Source this file.
#
# PROJ_ROOT is NOT resolved from ${BASH_SOURCE[0]} -- that breaks in SLURM work
# directories (see CLAUDE.md). It defaults to this repo's fixed /bigdata location
# and can be overridden with the PROJ_ROOT env var if the repo is cloned elsewhere.
PROJ_ROOT="${PROJ_ROOT:-/bigdata/stajichlab/jstajich/projects/adhesionPred}"
SITE_YAML="$PROJ_ROOT/config/site.yaml"

_site_get() {
  awk -F': ' -v key="$1" '$1==key {print $2; exit}' "$SITE_YAML"
}

FUNGI5K_INPUT=$(_site_get fungi5k_input)
FUNGI5K_DUCKDB=$(_site_get fungi5k_duckdb)
FUNGI5K_SAMPLES=$(_site_get fungi5k_samples)
COCCI_PANGENOME=$(_site_get cocci_pangenome)
COCCI_LONGREAD=$(_site_get cocci_longread)
SPHERULE_RNASEQ=$(_site_get spherule_rnaseq)
RHODOTORULA_BIOFILM=$(_site_get rhodotorula_biofilm)
WORKDIR=$(_site_get workdir)
