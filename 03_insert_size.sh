#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 10GB
#SBATCH -t 3:00:00
#SBATCH -J QC_I
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/QC_insert_adypocites.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/QC_insert_adypocites.out

module load picard/3.4.0-Java-17

output_folder=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qc

a=${1##*/} #without path
b=$(basename $a _unique_sorted.bam)

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics --Histogram_FILE $output_folder/${b}_histogram.pdf --INPUT $1 --OUTPUT $output_folder/${b}.txt
