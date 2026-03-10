#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 1:00
#SBATCH -J array06
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_06.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/array_06.out

input_path=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto

gtf_folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump

#for F in $input_path/100_top_windows_factor_*.bed; do
#        for i in $gtf_folder/promoter_TSS_region_mus_musculus_mm39.bed $gtf_folder/CDS_mus_musculus_mm39.sorted.merged.bed $gtf_folder/intronic_mus_musculus_mm39.bed $gtf_folder/3_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/5_UTR_mus_musculus_mm39.sorted.merged.bed $gtf_folder/start_codon_mus_musculus_mm39.sorted.merged.bed $gtf_folder/stop_codon_mus_musculus_mm39.sorted.merged.bed; do
#                echo $F
#                echo $i
#                sbatch 06_annotation_top_wind_from_factors_mefisto.sh $F $i
#        done
#done

#for F in $input_path/100_top_windows_factor_*.bed; do
#	sbatch 06.1_remap_annotation_mefisto.sh $F
#done

for F in factor_10 factor_8 factor_9; do
	sbatch 06.2_gene_names.sh $F
done
