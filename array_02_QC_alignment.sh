#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 1:00
#SBATCH -J array01
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_02.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_02.out

input_path=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/aligned

for F in $input_path/*_unique_sorted.bam; do
        echo $F
        sbatch --export=ALL,a=$F 02_alignment_QC.sh $F
done
