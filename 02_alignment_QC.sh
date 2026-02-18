#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 10GB
#SBATCH -t 3-00:00:00
#SBATCH -J QC_A
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/QC_Alignment_adypocites.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/QC_Alignment_adypocites.out

module load bwa-mem2/2.2.1-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0
module load picard/3.4.0-Java-17
module load FastQC/0.12.1-Java-11


#working directories, and file names
a=${1##*/} #without path
b=${a%_unique_sorted.bam} #without the .bam extension
sample_dir=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/aligned
output_dir=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qc
reference_genome=/proj/naiss2024-23-57/reference_genomes/mus_musculus/uscs_ref/mm39.fa

java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics -I $sample_dir/$a -O $output_dir/${b}_QC_picard --HISTOGRAM_FILE $output_dir/${b}_alignment_hist_read_length.pdf

samtools flagstat $sample_dir/${a} > $output_dir/${b}_alignment.flagstat
samtools idxstats $sample_dir/${a} > $output_dir/${b}_alignment.idxstats
fastqc -o $output_dir $sample_dir/$a
#for average coverage
#samtools depth $sample_dir/${a} | awk '{sum+=$3} END {print "Average coverage:", sum/NR}' > $output_dir/average_coverage_${b}.txt
#for breadth of coverage
#samtools depth -a $sample_dir/${a} | awk '{c++; if($3>0) total+=1} END {print "Coverage Breadth:", (total/c)*100 "%"}' > $output_dir/breadth_coverage_${b}.txt

#for after having a summary of all breadths
#for i in breadth_coverage_*; do cat $i >> Breadth_coverage_female_adipocites_multigen.txt; done
