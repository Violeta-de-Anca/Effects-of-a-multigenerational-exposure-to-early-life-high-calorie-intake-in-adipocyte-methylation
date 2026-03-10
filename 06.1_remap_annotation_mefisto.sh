#!/bin/bash -l
#SBATCH -A uppmax2026-1-34
##SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 30GB
#SBATCH -t 3-00:00:00
#SBATCH -J annot
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/annotation_remap_top_wind_from_mefisto.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/annotation_remap_top_wind_from_mefisto.out

annotated_folder=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs
mkdir -p $annotated_folder

module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0

sample=$1
feature=$2

a=${1##*/}
echo $a
prefix="100_top_windows_"
sufix=".bed"
base=${a#$prefix}
final=${base%$sufix}
echo $final

bedtools intersect -bed -wb -wa -a $sample -b /proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/remap/remap2022_eWAT_mm39.bed > $annotated_folder/mefisto_${final}_remap.bed


