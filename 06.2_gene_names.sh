#!/bin/bash -l
#SBATCH -A uppmax2026-1-34
##SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 30GB
#SBATCH -t 3-00:00:00
#SBATCH -J annot
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/gene_names_top_wind_from_mefisto.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_adipocytes/log_files/gene_names_top_wind_from_mefisto.out

annotated_folder=/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs
mkdir -p $annotated_folder

module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.22-GCC-13.3.0

a=$1 #factor_X
echo $a

#first we need to have all the genomic features together
cat $annotated_folder/mefisto_${a}_3_UTR.bed $annotated_folder/mefisto_${a}_5_UTR.bed $annotated_folder/mefisto_${a}_CDS.bed \
$annotated_folder/mefisto_${a}_intronic_mus_musculus_mm39.bed.bed $annotated_folder/mefisto_${a}_promoter_TSS_region_mus_musculus_mm39.bed.bed \
$annotated_folder/mefisto_${a}_start_codon.bed \
$annotated_folder/mefisto_${a}_stop_codon.bed | cut -f 8,9,10 | bedtools sort -i - | bedtools merge -i - -d 6 > $annotated_folder/${a}_UCSC_genes.bed

ucsc_input=/proj/naiss2024-23-57/reference_genomes/mus_musculus/gtf_files/UCSC_dump/mm39.ncbiRefSeq.gtf.gz

zcat $ucsc_input | bedtools intersect -a $annotated_folder/${a}_UCSC_genes.bed -b stdin -wb | cut -f 12 > $annotated_folder/${a}_UCSC_gene_names.bed

cut -f 11 $annotated_folder/mefisto_${a}_remap.bed > $annotated_folder/transcript_fact_${a}_remap_name.txt

cat $annotated_folder/${a}_UCSC_gene_names.bed $annotated_folder/transcript_fact_${a}_remap_name.txt > $annotated_folder/all_genetic_features_names_${a}.txt
