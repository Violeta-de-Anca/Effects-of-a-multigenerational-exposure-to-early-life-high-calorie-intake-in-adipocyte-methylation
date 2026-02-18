#!/bin/bash -l
##SBATCH -A uppmax2025-2-302
#SBATCH -A uppmax2025-2-536
#SBATCH -p pelle
#SBATCH --mem 256GB
#SBATCH -t 5-00:00:00
#SBATCH -J mefisto
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/mefisto_sex_run.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/mefisto_sex_run.out

module load R/4.4.2-gfbf-2024a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load R-bundle-CRAN/2024.11-foss-2024a

R --no-save --quiet < mefisto_C57_female_adypocites.R
