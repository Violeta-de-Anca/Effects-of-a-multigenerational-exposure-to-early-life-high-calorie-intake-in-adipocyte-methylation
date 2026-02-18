#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
#SBATCH -p pelle
#SBATCH --mem 256GB
#SBATCH -t 3-00:00:00
#SBATCH -J qsea
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/qsea_alignment_load.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/qsea_alignment_load.out

module load R/4.4.2-gfbf-2024a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load R-bundle-CRAN/2024.11-foss-2024a

R --no-save --quiet < qsea_C57_adipocytes_female.R
