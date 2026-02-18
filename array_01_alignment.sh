#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 1:00
#SBATCH -J array01
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_01.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_01.out

input_path=/proj/gbs_medip/TEMPLETON/Analysis/MK3850/Processed

for F in $input_path/Sample_WK-3850-AdipocyteMeDIP-72_R1_trimmed.fq.gz $input_path/Sample_WK-3850-AdipocyteMeDIP-74_R1_trimmed.fq.gz $input_path/Sample_WK-3850-AdipocyteMeDIP-75_R1_trimmed.fq.gz $input_path/Sample_WK-3850-AdipocyteMeDIP-76_R1_trimmed.fq.gz; do
	echo $F
	sbatch --export=ALL,a=$F 01_alignment.sh $F
done
