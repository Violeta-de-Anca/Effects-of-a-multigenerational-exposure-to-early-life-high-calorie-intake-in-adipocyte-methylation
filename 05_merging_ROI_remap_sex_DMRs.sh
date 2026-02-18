#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
#SBATCH -p pelle
#SBATCH --mem 100GB
#SBATCH -t 3-00:00:00
#SBATCH -J merge
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/merging_ROIs_remap.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/merging_ROIs_remap.out

module load BEDTools/2.31.1-GCC-13.3.0

ROIs_folder=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs
remap_folder=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/remap

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

#first merge in a single file all the ROIs
for i in $ROIs_folder/*; do
	touch $ROIs_folder/all_rois.bed
	cat $i >> $ROIs_folder/all_rois.bed
done

cat $remap_folder/remap2022_eWAT_mm39.bed >> $ROIs_folder/all_rois.bed

# now do the merge of all coordinates
sort -k1,1 -k2,2n -k6,6 $ROIs_folder/all_rois.bed > $ROIs_folder/all_rois_sorted.bed
bedtools merge -i $ROIs_folder/all_rois_sorted.bed -s -d 1 > $ROIs_folder/ROIs_remap_intersect.bed

rm $ROIs_folder/all_rois_sorted.bed
