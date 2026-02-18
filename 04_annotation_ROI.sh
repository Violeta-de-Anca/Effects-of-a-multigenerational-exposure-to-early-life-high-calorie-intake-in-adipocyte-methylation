#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 100GB
#SBATCH -t 3-00:00:00
#SBATCH -J annot
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/annotation_ROI_C57_fem.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/annotation_ROI_C57_fem.out

annotated_folder=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs
mkdir -p $annotated_folder

module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0

sample=$1
feature=$2

a=${1##*/}
base=$(basename $a _unique_sorted.bam)
echo $base

b=${2##*/}
b_feat=$(basename $b _mus_musculus_mm39.bed)
echo $b_feat

bedtools intersect -bed -wb -abam $sample -b $feature > $annotated_folder/${base}_${b_feat}.bed
