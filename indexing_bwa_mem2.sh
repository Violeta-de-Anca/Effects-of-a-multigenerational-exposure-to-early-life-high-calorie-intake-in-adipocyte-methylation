#!/bin/bash -l
#SBATCH -A uppmax2025-2-151
#SBATCH -p pelle
#SBATCH --mem 64GB
#SBATCH -t 10:00:00
#SBATCH -J mm39
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/indexing_mm39_bwa_mem2.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/indexing_mm39_bwa_mem2.out

module load bwa-mem2/2.2.1-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0

#working directories
reference_genome=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa


bwa-mem2 index $reference_genome

