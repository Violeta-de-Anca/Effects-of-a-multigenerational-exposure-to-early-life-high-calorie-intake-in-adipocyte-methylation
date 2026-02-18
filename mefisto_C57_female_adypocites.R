library(BiocManager)
library(basilisk)
# BiocManager::install("MOFA2", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
Sys.setenv(BASILISK_HOME="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/basilisk")
library(MOFA2, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(tidyverse)
library(pheatmap)
library(data.table)

setwd("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto")
set.seed(123)

#load the data
#load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/rpm_counts.rda")
#colnames(rpm_counts_C57)=sub("_rpm$","",colnames(rpm_counts_C57))
relational_table_C57=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_with_file_path.txt")
relational_table_C57[,sex := tolower(sex)]
relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)

#filter the windows to at least have 10 reads in at least half of the individuals (31)
#matrix_c57=as.matrix(rpm_counts_C57[,5:ncol(rpm_counts_C57)])
#keep=rowSums(matrix_c57>0) >=31
#rpm_counts_C57=rpm_counts_C57[keep,]

#so we need to put the dataframe in long format, there has to be only 3 columns, sample, feature and value
#long_rmp_counts_C57=pivot_longer(rpm_counts_C57,cols = c(5:76),values_to = "value", names_to = "sample")
#long_rmp_counts_C57$feature=paste0(long_rmp_counts_C57$chr,":",long_rmp_counts_C57$start,"-",long_rmp_counts_C57$end)
#long_rmp_counts_C57=long_rmp_counts_C57[,c("sample","feature","value")]
# sample_to_group=relational_table_C57%>%distinct(sample_name, group)%>%mutate(sample=sample_name)%>% select(sample,group)
# long_rmp_counts_C571=long_rmp_counts_C57%>%left_join(sample_to_group,by="sample")
# c57_mefisto=create_mofa(data=long_rmp_counts_C571,groups = "group")
# save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/first_c57_mefisto_object.rda")
# #this is to specifically see in between the C and SL differences
# relational_table_C571=relational_table_C57%>%
#   mutate(
#     generation=case_when(
#       generation=="F0"~0,
#       generation=="F1"~1,
#       generation=="F2"~2
#     ))
# relational_table_C571=relational_table_C571%>%mutate(family=as.numeric(family))
# long_relational_table_C57=pivot_longer(relational_table_C571,cols = c(5:6),values_to = "value", names_to = "covariate")
# long_relational_table_C57$sample=long_relational_table_C57$sample_name
# long_relational_table_C57=long_relational_table_C57[,c("sample","covariate","value")]
# c57_mefisto=set_covariates(c57_mefisto,covariates = as.data.frame(long_relational_table_C57))
# save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_object_with_covariates.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_object_with_covariates.rda")
# #now we need to set the options
 options_data=get_default_data_options(c57_mefisto)
 opts_model=get_default_model_options(c57_mefisto)
opts_model$likelihoods[]="poisson"
# #number of components, or factors as they call them here to calculate
 opts_model$num_factors=10
 opts_training=get_default_training_options(c57_mefisto)
 opts_training$maxiter=1000
 opts_training$convergence_mode="slow"
 opts_training$seed=123
 opts_mefisto=get_default_mefisto_options(c57_mefisto)
 opts_mefisto$model_groups=T
# 
 c57_mefisto=prepare_mofa(c57_mefisto,model_options = opts_model,
                          mefisto_options = opts_mefisto,
                          training_options = opts_training,
                          data_options = options_data)
# 
 c57_mefisto=run_mofa(c57_mefisto,
                      outfile = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mofa_run_ouput.hdf5",
                      use_basilisk = T)
 save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_ran_mofa.rda")

##################################################################################
##now let's specify the group for MOFA as the sex, we have seen that there is much differences with the sex in the QSEA analysis
#sample_to_group=relational_table_C57%>%distinct(sample_name, sex)%>%mutate(sample=sample_name)%>% select(sample,sex)
#long_rmp_counts_C571=long_rmp_counts_C57%>%left_join(sample_to_group,by="sample")
#c57_mefisto=create_mofa(data=long_rmp_counts_C571,groups = "sex")
#save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/first_c57_sex_groups_mefisto_object.rda")
#relational_table_C571=relational_table_C57%>%
#  mutate(
#    generation=case_when(
#      generation=="F0"~0,
#      generation=="F1"~1,
#      generation=="F2"~2
#    ))
#relational_table_C571=relational_table_C571%>%mutate(family=as.numeric(family))
#long_relational_table_C57=pivot_longer(relational_table_C571,cols = c(5:6),values_to = "value", names_to = "covariate")
#long_relational_table_C57$sample=long_relational_table_C57$sample_name
#long_relational_table_C57=long_relational_table_C57[,c("sample","covariate","value")]
#c57_mefisto=set_covariates(c57_mefisto,covariates = as.data.frame(long_relational_table_C57))
#save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_with_covariates.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_with_covariates.rda")
#now we need to set the options
options_data=get_default_data_options(c57_mefisto)
opts_model=get_default_model_options(c57_mefisto)
opts_model$likelihoods[]="poisson"

#number of components, or factors as they call them here to calculate
opts_model$num_factors=10
opts_training=get_default_training_options(c57_mefisto)
opts_training$maxiter=1000
opts_training$convergence_mode="slow"
opts_training$seed=123
opts_mefisto=get_default_mefisto_options(c57_mefisto)
opts_mefisto$model_groups=T

c57_mefisto=prepare_mofa(c57_mefisto,model_options = opts_model,
                         mefisto_options = opts_mefisto,
                         training_options = opts_training,
                         data_options = options_data)

c57_mefisto=run_mofa(c57_mefisto,
                     outfile = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mofa_sex_groups_run_ouput.hdf5",
                     use_basilisk = T)
save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_ran_mofa.rda")











