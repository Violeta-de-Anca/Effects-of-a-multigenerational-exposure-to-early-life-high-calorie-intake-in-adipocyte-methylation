#!/bin/bash -l
##SBATCH -A uppmax2025-2-302
#SBATCH -A uppmax2025-2-536
#SBATCH -p pelle
##SBATCH --mem 2GB
#SBATCH --mem 256GB
#SBATCH -t 3-00:00:00
#SBATCH -J mm39
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/alignment_mm39.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/alignment_mm39.out

#module load bwa-mem2/2.2.1-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0

#working directories, and file names
a=${1##*/} #without path
sample_dir=/proj/gbs_medip/TEMPLETON/Analysis/MK3850/Processed
output_dir=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/aligned
reference_genome=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa

#temp files for the sam files
temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

#aligning
base=$(basename $a _R1_trimmed.fq.gz)
sample=${base#Sample_WK-3850-AdipocyteMeDIP-}
input2=${base}_R2_trimmed.fq.gz
echo $sample
echo $base
echo $output_dir/${sample}_unique.bam
echo $sample_dir/$input2
/sw/arch/eb/software/bwa-mem2/2.3-GCC-13.3.0/bin/bwa-mem2 mem -t 10 -M -R "@RG\tID:GBSICR\tSM:$sample\tPL:ILLUMINA" $reference_genome $1 $sample_dir/$input2 > $temp_dir/$sample.sam
samtools view -S -b $temp_dir/$sample.sam > $output_dir/$sample.bam
samtools sort -m 768M $output_dir/$sample.bam > $output_dir/${sample}_sorted.bam
samtools index $output_dir/${sample}_sorted.bam $output_dir/${sample}_sorted.bam.bai
samtools view -h $output_dir/${sample}_sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > $output_dir/${sample}_unique.bam
samtools sort -m 768M $output_dir/${sample}_unique.bam > $output_dir/${sample}_unique_sorted.bam
samtools index $output_dir/${sample}_unique_sorted.bam $output_dir/${sample}_unique_sorted.bam.bai
