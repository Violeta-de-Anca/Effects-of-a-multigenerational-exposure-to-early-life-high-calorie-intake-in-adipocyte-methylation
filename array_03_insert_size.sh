#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 1:00
#SBATCH -J array03
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_03.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_03.out

input_path=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/aligned

for F in $input_path/*_unique_sorted.bam; do
        echo $F
        sbatch --export=ALL,a=$F 03_insert_size.sh $F
done


