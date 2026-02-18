#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 1:00
#SBATCH -J array04
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_04.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_04.out

input_path=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/aligned

gtf_folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump

for F in $input_path/*_unique_sorted.bam; do
	for i in $gtf_folder/promoter_TSS_region_mus_musculus_mm39.bed $gtf_folder/CDS_mus_musculus_mm39.sorted.merged.bed $gtf_folder/intronic_mus_musculus_mm39.bed $gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/start_codon_mus_musculus_mm39.sorted.merged.bed $gtf_folder/stop_codon_mus_musculus_mm39.sorted.merged.bed; do
		echo $F
		echo $i
		sbatch 04_annotation_ROI.sh $F $i
	done
done
