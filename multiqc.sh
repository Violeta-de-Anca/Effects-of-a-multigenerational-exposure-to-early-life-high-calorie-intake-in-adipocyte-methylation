#!/bin/bash -l
#SBATCH -A uppmax2025-2-151
#SBATCH -p core -n 1
#SBATCH -t 10:00:00
#SBATCH -J multiQC
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/multiQC_trimmed_adipocytes.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/multiQC_trimmed_adipocytes.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load MultiQC/1.9
module load samtools/1.14
module load FastQC/0.11.9

fastqc_dir=/proj/gbs_medip/TEMPLETON/Analysis/MK3850/Processed/FastQC_Trimmed
multiqc_dir=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qc

mkdir -p $multiqc_dir

#do a multiQC of the trimmed adipocytes


multiqc --export -outdir $multiqc_dir  $fastqc_dir/*
