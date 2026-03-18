library(BiocManager)
library(basilisk)
# BiocManager::install("MOFA2", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
Sys.setenv(BASILISK_HOME="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/basilisk")
library(MOFA2, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(tidyverse)
library(plyr)
library(pheatmap)
library(purrr)
library(tibble)
library(data.table)
# BiocManager::install("karyoploteR", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(GenomicRanges)
library(karyoploteR, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
require(devtools)
# devtools::install_github("MiguelCastresana/anubix", lib="/crex/proj/naiss2024-23-57/reference_genomes")
library(ANUBIX , lib.loc = "/crex/proj/naiss2024-23-57/reference_genomes")
# install.packages("neat", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(neat, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
# install.packages("gprofiler2", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(gprofiler2, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(KEGGgraph)
library(KEGGREST)
# devtools::install_github("noriakis/ggkegg", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(ggkegg, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
# install.packages("ggfx", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(ggfx, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(ggraph)
library(igraph)
library(clusterProfiler)
library(tidygraph)
library(msigdbr)
library(org.Mm.eg.db)
# install.packages("pathfindR", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(pathfindR, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(biomaRt)
library(UpSetR)
setwd("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto")
set.seed(123)

#load the data
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/rpm_counts.rda")
# colnames(rpm_counts_C57)=sub("_rpm$","",colnames(rpm_counts_C57))
relational_table_C57=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_with_file_path.txt")
relational_table_C57[,sex := tolower(sex)]
relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)
# 
# #filter the windows to at least have 10 reads in at least half of the individuals (31)
# matrix_c57=as.matrix(rpm_counts_C57[,5:ncol(rpm_counts_C57)])
# keep=rowSums(matrix_c57>0) >=31
# rpm_counts_C57=rpm_counts_C57[keep,]
# 
# #so we need to put the dataframe in long format, there has to be only 3 columns, sample, feature and value
# long_rmp_counts_C57=pivot_longer(rpm_counts_C57,cols = c(5:76),values_to = "value", names_to = "sample")
# long_rmp_counts_C57$feature=paste0(long_rmp_counts_C57$chr,":",long_rmp_counts_C57$start,"-",long_rmp_counts_C57$end)
# long_rmp_counts_C57=long_rmp_counts_C57[,c("sample","feature","value")]
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
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_object_with_covariates.rda")
# # #now we need to set the options
#  options_data=get_default_data_options(c57_mefisto)
# # gaussian if you are using normalized counts, poisson for raw counts
#  opts_model=get_default_model_options(c57_mefisto)
# # opts_model$likelihoods[]="poisson"
# # #number of components, or factors as they call them here to calculate
#  opts_model$num_factors=10
#  opts_training=get_default_training_options(c57_mefisto)
#  opts_training$maxiter=1000
#  opts_training$convergence_mode="slow"
#  opts_training$seed=123
#  opts_mefisto=get_default_mefisto_options(c57_mefisto)
#  opts_mefisto$model_groups=T
# # 
#  c57_mefisto=prepare_mofa(c57_mefisto,model_options = opts_model,
#                           mefisto_options = opts_mefisto,
#                           training_options = opts_training,
#                           data_options = options_data)
# # 
#  c57_mefisto=run_mofa(c57_mefisto,
#                       outfile = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mofa_run_ouput.hdf5",
#                       use_basilisk = T)
#  save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_ran_mofa.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_ran_mofa.rda")
# 
# ########################################################333
# #downstream analysis, in this one we should get the sex differences ####
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_visualization_by_treatment.tiff", height = 2000, width = 2500, res = 150)
# plot_variance_explained(c57_mefisto)
# dev.off()
# #get the variance that each factor explains ####
# var_factors_treatment_grouping=get_variance_explained(c57_mefisto)
# 
# plot_factor_cor(c57_mefisto)
# #check if the components variation also varies during the different generations
# get_scales(c57_mefisto)
# 
# #let's plot the factors with the most variance explained and plot females and males to see the direction #####
# md=samples_metadata(c57_mefisto)
# sex_df=relational_table_C57%>%transmute(sample=sample_name,sex=sex)
# md=md%>%left_join(sex_df,by="sample")
# md$sex=factor(md$sex,levels = c("female","male"))
# samples_metadata(c57_mefisto)=md
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_evolution_treatment_grouping_color_sex.tiff", height = 2000, width = 2500, res = 150)
# plot_factors_vs_cov(c57_mefisto,covariates = "generation",color_by = "sex",factors = 1)
# dev.off()
# 
# #distribution of the factors per sex ####
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_evolution_treatment_grouping_color_sex.tiff", height = 2000, width = 2500, res = 150)
# plot_factors(c57_mefisto,color_by = "sex",factors = c(1,2,3,8))
# dev.off()

##################################################################################
##now let's specify the group for MOFA as the sex, we have seen that there is much differences with the sex in the QSEA analysis
# sample_to_group=relational_table_C57%>%distinct(sample_name, sex)%>%mutate(sample=sample_name, group=sex)%>% select(sample,group)
# long_rmp_counts_C571=long_rmp_counts_C57%>%left_join(sample_to_group,by="sample")
# c57_mefisto=create_mofa(data=long_rmp_counts_C571,groups = "group")
# save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/first_c57_sex_groups_mefisto_object.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/first_c57_sex_groups_mefisto_object.rda")
relational_table_C571=relational_table_C57%>%
 mutate(
   generation=case_when(
     generation=="F0"~0,
     generation=="F1"~1,
     generation=="F2"~2
   ),
   group=case_when(
     group=="C"~0,
     group=="SL"~1
   ))
relational_table_C571=relational_table_C571%>%mutate(family=as.numeric(family))
long_relational_table_C57=pivot_longer(relational_table_C571,cols = c(3,5:6),values_to = "value", names_to = "covariate")
long_relational_table_C57$sample=long_relational_table_C57$sample_name
long_relational_table_C57=long_relational_table_C57[,c("sample","covariate","value")]
c57_mefisto=set_covariates(c57_mefisto,covariates = as.data.frame(long_relational_table_C57))
save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_with_covariates.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_with_covariates.rda")
#now we need to set the options
options_data=get_default_data_options(c57_mefisto)
# gaussian if you are using normalized counts, poisson for raw counts
opts_model=get_default_model_options(c57_mefisto)
# opts_model$likelihoods[]="poisson"

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

load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_ran_mofa.rda")

#downstream analysis, in here we should get the treatment differences ####
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_visualization_by_sex.tiff", height = 2000, width = 2500, res = 150)
plot_variance_explained(c57_mefisto)
dev.off()
#get the variance that each factor explains ####
var_factors_sex_grouping=get_variance_explained(c57_mefisto)

plot_factor_cor(c57_mefisto)
#check if the components variation also varies during the different generations
get_scales(c57_mefisto)

#let's plot the factors with the most variance explained and plot females and males to see the direction #####
md=samples_metadata(c57_mefisto)
sex_df=relational_table_C57%>%transmute(sample=sample_name,sex=sex)
md=md%>%left_join(sex_df,by="sample")
md$sex=factor(md$sex,levels = c("female","male"))
samples_metadata(c57_mefisto)=md
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_evolution_treatment_grouping_color_sex.tiff", height = 2000, width = 2500, res = 150)
plot_factors_vs_cov(c57_mefisto,covariates = "generation",color_by = "group_scaled",factors = 1)
dev.off()

#distribution of the factors per sex ####
female_57=subset_groups(c57_mefisto,groups = "female")
male_57=subset_groups(c57_mefisto,groups = "male")
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_evolution_sex_grouping_color_treatment_onlyfemales.tiff", height = 2000, width = 2500, res = 150)
plot_factors(female_57,color_by = "group_scaled",factors = 1:10)
dev.off()
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_evolution_sex_grouping_color_treatmentallfactors_onlymales.tiff", height = 2000, width = 2500, res = 150)
plot_factors(male_57,color_by = "group_scaled",factors = 1:10)
dev.off()

#get the scores of the factors so we can do our own density plots and test ####
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_ran_mofa.rda")
relational_table_C57=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_with_file_path.txt")
relational_table_C57[,sex := tolower(sex)]
relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)
factors_sex_grouping=get_factors(c57_mefisto)
female_factors_sex_grouping=as.data.frame(factors_sex_grouping["female"])
female_factors_sex_grouping=female_factors_sex_grouping%>%mutate(sample_name=row.names(female_factors_sex_grouping))
female_factors_sex_grouping=left_join(female_factors_sex_grouping,relational_table_C57,by = "sample_name")
male_factors_sex_grouping=as.data.frame(factors_sex_grouping["male"])
male_factors_sex_grouping=male_factors_sex_grouping%>%mutate(sample_name=row.names(male_factors_sex_grouping))
male_factors_sex_grouping=left_join(male_factors_sex_grouping,relational_table_C57,by = "sample_name")
#mu_4=ddply(male_factors_sex_grouping,"group",summarise,grp.mean=mean(male.Factor4))
# male_factors_sex_grouping=pivot_longer(male_factors_sex_grouping,cols = 1:10,values_to = "score",names_to = "Factors")

tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factor_10_male_density.tiff", height = 1500, width = 2500, res = 150)
male_factor_density=ggplot(male_factors_sex_grouping,aes(x=male.Factor10, color=group,fill = group))+
  geom_density(alpha=0.4,linewidth= 3)+
  theme_classic()+theme(
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 28),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28)
  )+
  scale_colour_manual(values = c("#85bc37","#e74269"),
                      breaks = c("C","SL"),
                      labels=c("Control","Small Litter"))+
  labs(x="Male factor 10",
       y="Density of the score",
       color="Treatment")+
  scale_fill_manual(values = c("#85bc37","#e74269"),
                    breaks = c("C","SL"),
                    labels=c("Control","Small Litter"),
                    guide="none")
male_factor_density
dev.off()

#now we are going to see which of the factor, for each of the sex, have different distributions, with the kolmogorov-smirnov test ####
ks_dist=function(df,group_col = "group", g1= "C", g2= "SL",
                 factor_pattern= "\\.Factor\\d+$"){
  factor_cols=grep(factor_pattern,names(df),value = T)
  map_dfr(factor_cols,function(col){
    x=df%>% filter(.data[[group_col]]==g1)%>%pull(.data[[col]])%>%na.omit()
    y=df%>% filter(.data[[group_col]]==g2)%>%pull(.data[[col]])%>%na.omit()
    kt=suppressWarnings(stats::ks.test(x,y,exact=F))
    tibble(
      factor=col,
      n_C=length(x),
      n_SL=length(y),
      D= unname(kt$statistic),
      p_value=kt$p.value
    )
  })
}
female_ks_dist=ks_dist(female_factors_sex_grouping)
male_ks_dist=ks_dist(male_factors_sex_grouping)

#now we are going to see which of the factor, for each of the sex, have the same distributions, but different medians, with the mann whitney test ####
mw_dist=function(df,group_col = "group", g1= "C", g2= "SL",
                 factor_pattern= "\\.Factor\\d+$"){
  factor_cols=grep(factor_pattern,names(df),value = T)
  map_dfr(factor_cols,function(col){
    x=df%>% filter(.data[[group_col]]==g1)%>%pull(.data[[col]])%>%na.omit()
    y=df%>% filter(.data[[group_col]]==g2)%>%pull(.data[[col]])%>%na.omit()
    mw=suppressWarnings(stats::wilcox.test(x,y))
    tibble(
      factor=col,
      n_C=length(x),
      n_SL=length(y),
      D= unname(mw$statistic),
      p_value=mw$p.value
    )
  })
}
female_mw_dist=ks_dist(female_factors_sex_grouping)
male_mw_dist=mw_dist(male_factors_sex_grouping)

#plot the only female factor that is borderline significant ####
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factor_9_female_density.tiff", height = 1500, width = 2500, res = 150)
female_factor_density=ggplot(female_factors_sex_grouping,aes(x=female.Factor9, color=group,fill = group))+
  geom_density(alpha=0.4,linewidth= 3)+
  theme_classic()+theme(
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 28),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28)
  )+
  scale_colour_manual(values = c("#85bc37","#e74269"),
                      breaks = c("C","SL"),
                      labels=c("Control","Small Litter"))+
  labs(x="Female factor 9",
       y="Density of the score",
       color="Treatment")+
  scale_fill_manual(values = c("#85bc37","#e74269"),
                    breaks = c("C","SL"),
                    labels=c("Control","Small Litter"),
                    guide="none")
female_factor_density
dev.off()

# now we need to see if the significant factors go in the same direction or not ####
# we know they are independent of time ####
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factor_correlation_female.tiff", height = 1500, width = 1500, res = 150)
plot_factor_cor(female_57,"spearman")
dev.off()

tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factor_correlation_male.tiff", height = 1500, width = 1500, res = 150)
plot_factor_cor(male_57,"spearman")
dev.off()

#get the top 100 windows which contributes the most to the relevant factors ####
male_features=get_weights(male_57,as.data.frame = T)
male_features=male_features%>%filter(factor%in% c("Factor8", "Factor9", "Factor10"))
male_features=male_features%>%pivot_wider(names_from = "factor",values_from = "value")
factors=c("Factor8", "Factor9", "Factor10")
top_features_males=factors%>%set_names()%>%
  map(~male_features%>%
        transmute(
          feature,
          value =.data[[.x]]
        )%>%
        arrange(desc(abs(value)))%>%
        slice_head(n=100)
      )
save(top_features_males,file="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_features_males_relevant_factors.rda")

tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_weight_sign_factor10_male.tiff",
     height = 3000, width = 1500, res = 150)
plot_top_weights(male_57,nfeatures = 100,factors = 10)
dev.off()

female_features=get_weights(female_57,as.data.frame = T)
female_features=female_features%>%filter(factor%in% c("Factor9"))
female_features=female_features%>%pivot_wider(names_from = "factor",values_from = "value")
factors=c("Factor9")
top_features_females=factors%>%set_names()%>%
  map(~female_features%>%
        transmute(
          feature,
          value =.data[[.x]]
        )%>%
        arrange(desc(abs(value)))%>%
        slice_head(n=100)
  )
save(top_features_females,file="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_features_females_relevant_factors.rda")

tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_weight_sign_factor9_female.tiff",
     height = 3000, width = 1500, res = 150)
plot_top_weights(female_57,nfeatures = 100,factors = 9)
dev.off()

load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_features_males_relevant_factors.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_features_females_relevant_factors.rda")

#let's plot the top weights to see if we have more incidence in a chromosome #####
factor8=as.data.frame(top_features_males$Factor8)
factor8=factor8%>%mutate(feature=as.character(feature))%>%
  separate(feature,into = c("chr","pos"),sep = ":",remove = F)%>%
  separate(pos, into = c("start","end"),sep = "-",convert = T)
factor8$factor="factor_8"
gr8=GRanges(seqnames = factor8$chr,ranges = IRanges(start = factor8$start,end = factor8$end),
            factor=factor8$factor,value=factor8$value)
#guardalo como tabla para anotar
bed_8=data.frame(
  chrom=factor8$chr,
  chromStart=factor8$start,
  chromEnd=factor8$end,
  score=0,
  strand=".",
  value=factor8$value,
  factor=factor8$factor
)
write.table(bed_8,file = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/100_top_windows_factor_8.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)


factor9=as.data.frame(top_features_males$Factor9)
factor9=factor9%>%mutate(feature=as.character(feature))%>%
  separate(feature,into = c("chr","pos"),sep = ":",remove = F)%>%
  separate(pos, into = c("start","end"),sep = "-",convert = T)
factor9$factor="factor_9"
gr9=GRanges(seqnames = factor9$chr,ranges = IRanges(start = factor9$start,end = factor9$end),
            factor=factor9$factor,value=factor9$value)
#guardalo como tabla para anotar
bed_9=data.frame(
  chrom=factor9$chr,
  chromStart=factor9$start,
  chromEnd=factor9$end,
  score=0,
  strand=".",
  value=factor9$value,
  factor=factor9$factor
)
write.table(bed_9,file = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/100_top_windows_factor_9.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

factor10=as.data.frame(top_features_males$Factor10)
factor10=factor10%>%mutate(feature=as.character(feature))%>%
  separate(feature,into = c("chr","pos"),sep = ":",remove = F)%>%
  separate(pos, into = c("start","end"),sep = "-",convert = T)
factor10$factor="factor_10"
gr10=GRanges(seqnames = factor10$chr,ranges = IRanges(start = factor10$start,end = factor10$end),
            factor=factor10$factor,value=factor10$value)
#guardalo como tabla para anotar
bed_10=data.frame(
  chrom=factor10$chr,
  chromStart=factor10$start,
  chromEnd=factor10$end,
  score=0,
  strand=".",
  value=factor10$value,
  factor=factor10$factor
)
write.table(bed_10,file = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/100_top_windows_factor_10.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_weight_genome_loc_male.tiff",
     height = 3000, width = 1500, res = 150)
cyto=getCytobands(genome = "mm39")
cyto$gieStain[cyto$name=="p"]="pArm"
cyto$gieStain[cyto$name=="q"]="qArm"
pp=getDefaultPlotParams(plot.type = 2)
pp$data2height=30
bw=c(
  gneg= "white",
  qArm= "grey80",
  gpos50 = "grey60",
  gpos75 = "grey40",
  gpos100 = "black",
  gvar = "grey70",
  stalk = "grey50",
  acen = "black"
)
kp=plotKaryotype(genome = "mm39",plot.type = 1,plot.params = pp,
                 ideogram.plotter = NULL,cytobands = mm10,cex=2)
kpAddCytobands(kp,color.table = bw,color.schema = "biovizbase")
kpAddBaseNumbers(kp)
kpPlotRegions(kp,data = gr8,col = "#0072B2",border = "#0072B2",avoid.overlapping = T)
kpPlotRegions(kp,data = gr9,col = "#D55E00",border = "#D55E00",avoid.overlapping = T)
kpPlotRegions(kp,data = gr10,col = "#009E73",border = "#009E73",avoid.overlapping = T)
legend("right",
       legend = c("Factor 8","Factor 9","Factor 10"),
       fill = c("#0072B2","#D55E00","#009E73"),
       cex = 2,
       bty = "n")

dev.off()

#get also the heatmaps of the counts to show#####

plot_data_heatmap(male_57,view = "single_view",factor = 8,features = 100, show_rownames =F,
                  show_colnames=T,cluster_rows=T,cluster_cols=F,max.value = 0.10)
tiff(filename = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/top_weight_factor9_male_counts_heatmap.tiff",
     height = 1500, width = 1700, res = 150)
plot_data_heatmap(male_57,view = "single_view",factor = 9,features = 100, show_rownames =F,
                  show_colnames=T,cluster_rows=T,cluster_cols=T,annotation_samples = "group_scaled")
dev.off()
plot_data_heatmap(male_57,view = "single_view",factor = 10,features = 100, show_rownames =F,
                  show_colnames=T,cluster_rows=T,cluster_cols=F)
plot_data_heatmap(female_57,view = "single_view",factor = 9,features = 100, show_rownames =F,
                  show_colnames=T,cluster_rows=T,cluster_cols=F)

#we did with bedtools the gene annotation and the transcript factor annotation #####
#import the gene list of all factors
gene_names_factor_8=fread("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs/final_unique_all_gene_names_TF_UCSC_mm39_factor_8.txt",col.names = "SYMBOL")
gene_names_factor_9=fread("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs/final_unique_all_gene_names_TF_UCSC_mm39_factor_9.txt",col.names = "SYMBOL")
gene_names_factor_10=fread("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/annotated_ROIs/final_unique_all_gene_names_TF_UCSC_mm39_factor_10.txt",col.names = "SYMBOL")

#we then need to switch to uniprot IDs, this are the reviewed #####
listEnsembl()
ensembl=useEnsembl(biomart = "ensembl")
listDatasets(mart = ensembl)
grep("unipr",listAttributes(mart = ensembl_id))
ensembl_id=useEnsembl("genes",dataset = "mmusculus_gene_ensembl")
uniprot_F8=getBM(
  attributes = c("mgi_symbol","uniprotswissprot"),
  filters = "mgi_symbol",
  values = gene_names_factor_8$SYMBOL,
  mart = ensembl_id
)%>%filter(uniprotswissprot!="")%>%dplyr::rename(
  SYMBOL=mgi_symbol,
  UNIPROT=uniprotswissprot
)

uniprot_F9=getBM(
  attributes = c("mgi_symbol","uniprotswissprot"),
  filters = "mgi_symbol",
  values = gene_names_factor_9$SYMBOL,
  mart = ensembl_id
)%>%filter(uniprotswissprot!="")%>%dplyr::rename(
  SYMBOL=mgi_symbol,
  UNIPROT=uniprotswissprot
)

uniprot_F10=getBM(
  attributes = c("mgi_symbol","uniprotswissprot"),
  filters = "mgi_symbol",
  values = gene_names_factor_10$SYMBOL,
  mart = ensembl_id
)%>%filter(uniprotswissprot!="")%>%dplyr::rename(
  SYMBOL=mgi_symbol,
  UNIPROT=uniprotswissprot
)


# Then for the gene names we don't get in the review uniprot database ######
missing_f8=gene_names_factor_8%>%distinct(SYMBOL)%>%anti_join(uniprot_F8%>%distinct(SYMBOL),by = c("SYMBOL"))
not_primary_F8=AnnotationDbi::select(org.Mm.eg.db,keys = as.character(missing_f8$SYMBOL),
                  columns = c("UNIPROT","SYMBOL"),keytype = "SYMBOL")%>%distinct(UNIPROT,SYMBOL,.keep_all = T)%>%
  filter(!is.na(UNIPROT))
to_add_f8=not_primary_F8%>%group_by(SYMBOL)%>%filter(
  case_when(
    SYMBOL == "Raph1" ~T,
    SYMBOL %in% c("Il1rapl1","Lhfpl3","Sntg1") ~ row_number() == 2,
    SYMBOL == "Lrrc70" ~UNIPROT == "Q80TE7",
    SYMBOL == "Scml2" ~ UNIPROT == "B1AVB3",
    T ~ row_number() == 1
  )
)%>%ungroup()
uniprot_F8=rbind(uniprot_F8,to_add_f8)
uniprot_F8$factor="factor_8"


missing_f9=gene_names_factor_9%>%distinct(SYMBOL)%>%anti_join(uniprot_F9%>%distinct(SYMBOL),by = c("SYMBOL"))
not_primary_F9=AnnotationDbi::select(org.Mm.eg.db,keys = as.character(missing_f9$SYMBOL),
                      columns = c("UNIPROT","SYMBOL"),keytype = "SYMBOL")%>%distinct(UNIPROT,SYMBOL,.keep_all = T)%>%
  filter(!is.na(UNIPROT))
to_add_f9=not_primary_F9%>%group_by(SYMBOL)%>%filter(
  case_when(
    SYMBOL %in% c("Il1rapl1") ~ row_number() == 2,
    T ~ row_number() == 1
  )
)%>%ungroup()
uniprot_F9=rbind(uniprot_F9,to_add_f9)
uniprot_F9$factor="factor_9"

missing_f10=gene_names_factor_10%>%distinct(SYMBOL)%>%anti_join(uniprot_F10%>%distinct(SYMBOL),by = c("SYMBOL"))
not_primary_F10=AnnotationDbi::select(org.Mm.eg.db,keys = as.character(missing_f10$SYMBOL),
                      columns = c("UNIPROT","SYMBOL"),keytype = "SYMBOL")%>%distinct(UNIPROT,SYMBOL,.keep_all = T)%>%
  filter(!is.na(UNIPROT))
to_add_f10=not_primary_F10%>%group_by(SYMBOL)%>%filter(
  case_when(
    SYMBOL == "Abi3bp" ~ UNIPROT == "A0A338P6S8",
    SYMBOL == "Opcml" ~ UNIPROT == "G5E8G3",
    SYMBOL == "Nlrp4g" ~ UNIPROT == "",
    T ~ row_number() == 1
  )
)%>%ungroup()
uniprot_F10=rbind(uniprot_F10,to_add_f10)
uniprot_F10$factor="factor_10"

# funcoup=fread("/proj/naiss2024-23-57/reference_genomes/neural_network/FC6.0_M.musculus_full")
funcoup_compact=fread("/proj/naiss2024-23-57/reference_genomes/neural_network/FC6.0_M.musculus_compact")
#filtar solo las de calidad, min 0.95, high confidence
# funcoup_compact=funcoup_compact%>%filter(`5:PPV`>0.9)

#and now we can start with the enrichment analysis #######
factor_list=list(
  F8=uniprot_F8[,c("UNIPROT","factor")],
  F9=uniprot_F9[,c("UNIPROT","factor")],
  F10=uniprot_F10[,c("UNIPROT","factor")]
)
db_list=c("CP:KEGG","CP:REACTOME","CP:WIKIPATHWAYS")
funcoup=as.data.frame(funcoup_compact[,c(1:2,6)])

prepare_pathways=function(subcat){
  pathways=msigdbr(species = "Mus musculus",category ="C2",subcategory = subcat)%>%
    distinct(gs_name,gene_symbol,.keep_all = T)
  mapped_symbols_to_unipro=getBM(
    attributes = c("mgi_symbol","uniprotswissprot"),
    filters = "mgi_symbol",
    values = pathways$gene_symbol,
    mart = ensembl_id
  )%>%filter(uniprotswissprot!="")%>%dplyr::rename(
    gene_symbol=mgi_symbol,
    UNIPROT=uniprotswissprot
  )
  pathways=pathways%>%left_join(mapped_symbols_to_unipro,by="gene_symbol",relationship = "many-to-many")
  uniprot_names_pathways=pathways[,c("UNIPROT","gs_name")]%>%
    distinct(UNIPROT,gs_name,.keep_all = T)%>%filter(!is.na(UNIPROT))
  return(uniprot_names_pathways)
}

results= list()

for(db in db_list){
  cat("Running:",db,"\n")
  uniprot_names_pathways=prepare_pathways(db)
  anubix_links_matrix=ANUBIX::anubix_links(
    network = funcoup,
    pathways = as.data.frame(uniprot_names_pathways)
  )
  cat("Anubix finished","\n")
  results[[db]]=list()
  for(f_name in names(factor_list)){
    cat(" Running", f_name,"\n")
    geneset_df=factor_list[[f_name]]
    trans_res=anubix_transitivity(network= funcoup, 
                                  pathways= as.data.frame(uniprot_names_pathways),
                                  links_matrix = anubix_links_matrix, 
                                  genesets = geneset_df)%>%
      filter(`q-value`<0.05)
    cat(" Running normal anubix","\n")
    anubix_res=anubix(network= funcoup, 
                      pathways= as.data.frame(uniprot_names_pathways),
                      links_matrix = anubix_links_matrix, 
                      genesets = geneset_df)%>%
      filter(`q-value`<0.05)
    results[[db]][[f_name]]=list(transitivity=trans_res,
                                 anubix_normal=anubix_res)
  }
  cat("Finished","\n")
}

save(results,file = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/anubix_kegg_reactome_wikipathways_all_factors_annotated.rda")
load("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/anubix_kegg_reactome_wikipathways_all_factors_annotated.rda")

enrichment_f8=fread("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/FunCoup_enrichment_factor8.tsv")
enrichment_f9=fread("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/FunCoup_enrichment_factor9.tsv")
enrichment_f10=fread("/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/FunCoup_enrichment_factor10.tsv")

enrichment_f10$factor="factor_10"
enrichment_f9$factor="factor_9"
enrichment_f8$factor="factor_8"

enrichment_all_factors=full_join(enrichment_f8,enrichment_f9,by = "Pathway_ID")
enrichment_all_factors=full_join(enrichment_all_factors,enrichment_f10,by = "Pathway_ID")

#plot the enrichment score by pathway ID ####
pathway_common_all_f=enrichment_all_factors%>%filter(if_all(everything(), ~ !is.na(.)))


enrichment_all_factors_1=enrichment_all_factors%>%transmute(
  Pathway_ID,
  Pathway_Name_1=Pathway_Name,
  Pathway_Name_2= Pathway_Name.x,
  Pathway_Name_3=Pathway_Name.y,
  factor_1=factor,
  factor_2=factor.x,
  factor_3=factor.y,
  ANUBIX_FDR_1=ANUBIX_FDR,
  ANUBIX_FDR_2=ANUBIX_FDR.x,
  ANUBIX_FDR_3=ANUBIX_FDR.y
)%>%
  mutate(
  Pathway_Name = coalesce(Pathway_Name_1,Pathway_Name_2,Pathway_Name_3))%>%
  dplyr::select(-Pathway_Name_1,-Pathway_Name_2,-Pathway_Name_3)%>%pivot_longer(
    cols = c(factor_1,factor_2,factor_3,ANUBIX_FDR_1,ANUBIX_FDR_2,ANUBIX_FDR_3),
    names_to = c(".value","set"),
    names_pattern = "(factor|ANUBIX_FDR)_(\\d)"
  )%>%
  dplyr::select(Pathway_ID,Pathway_Name,factor,ANUBIX_FDR)%>%filter(!is.na(factor))

all_paths=ggplot(enrichment_all_factors_1,aes(x=ANUBIX_FDR,y=Pathway_Name))+
  geom_point(aes(color=factor),size=5)+
  labs(x="FDR",y="",title = "",color="Significant factors\nfrom MEFISTO")+
  theme_minimal()+
  scale_x_reverse()+
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.title =  element_text(size = 18)
  )+
  scale_colour_manual(values = c("factor_10" = "#0072B2","factor_9" = "#D55E00", "factor_8" = "#009E73"))

tiff(filename = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/enrichment_fdr_all_factors.tiff",
     height = 4000, width = 2500, res = 150)
all_paths
dev.off()


overlap_all=list(
  factor_8=enrichment_f8$Pathway_Name,
  factor_9=enrichment_f9$Pathway_Name,
  factor_10=enrichment_f10$Pathway_Name
)

tiff(filename = "/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/enrichment_overlapp_all_factors.tiff",
     height = 1500, width = 1700, res = 150)
upset(fromList(overlap_all),order.by = "freq",text.scale = 3)
dev.off()

#overlapp factor 8 with factor 9
overlap_f8andf9=inner_join(enrichment_f8,enrichment_f9,by = "Pathway_ID")
overlap_f8andf10=inner_join(enrichment_f8,enrichment_f10,by = "Pathway_ID")
overlap_f10andf9=inner_join(enrichment_f10,enrichment_f9,by = "Pathway_ID")

#Now I want to plot the insulin and metabolism related pathways: ####
mm39_pathways=get_gene_sets_list(org_code =  "mmu")
# 4911 (Insulin secretion) this is in overlapp 8 and 10 #### 
insulin_4911 = names(mm39_pathways$descriptions)[mm39_pathways$descriptions == "Insulin secretion"]
insulin_4911_id=mm39_pathways$gene_sets[[insulin_4911]]
insulin_4911_id=sub("mmu:","",insulin_4911_id)
#get the genes in uniprot ID
insulin_4911_genes=bitr(insulin_4911_id,
                        fromType = "ENTREZID",
                        toType = "UNIPROT",
                        OrgDb = org.Mm.eg.db)

funcoup_uniprotA=data.frame(UNIPROT = funcoup_compact$`0:ProteinA`,
                            other_UNIPROT = funcoup_compact$`1:ProteinB`)
insulin_4911_unipro_funcoupA=left_join(insulin_4911_genes,funcoup_uniprotA,by = "UNIPROT")
#which genes are affected in factor 8
links_insulin_4911_f8=uniprot_F8%>%inner_join(insulin_4911_unipro_funcoupA%>%
                                            pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                         names_to = "link_source",
                                                         values_to = "uniprot_genes")%>%
                                            distinct(uniprot_genes,.keep_all = T),
                                          by = c("UNIPROT" = "uniprot_genes"))
#which genes are affected in factor 10
links_insulin_4911_f10=uniprot_F10%>%inner_join(insulin_4911_unipro_funcoupA%>%
                                                pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                             names_to = "link_source",
                                                             values_to = "uniprot_genes")%>%
                                                distinct(uniprot_genes,.keep_all = T),
                                              by = c("UNIPROT" = "uniprot_genes"))

# 4931 (Insulin resistance) this is in overlapp 8 and 9 #####
insulin_4931 = names(mm39_pathways$descriptions)[mm39_pathways$descriptions == "Insulin resistance"]
insulin_4931_id=mm39_pathways$gene_sets[[insulin_4931]]
insulin_4931_id=sub("mmu:","",insulin_4931_id)
insulin_4931_genes=bitr(insulin_4931_id,
                        fromType = "ENTREZID",
                        toType = "UNIPROT",
                        OrgDb = org.Mm.eg.db)
insulin_4931_unipro_funcoupA=left_join(insulin_4931_genes,funcoup_uniprotA,by = "UNIPROT")
#which genes are affected in factor 8
links_insulin_4931_f8=uniprot_F8%>%inner_join(insulin_4931_unipro_funcoupA%>%
                                                pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                             names_to = "link_source",
                                                             values_to = "uniprot_genes")%>%
                                                distinct(uniprot_genes,.keep_all = T),
                                              by = c("UNIPROT" = "uniprot_genes"))
#which genes are affected in factor 9
links_insulin_4931_f9=uniprot_F9%>%inner_join(insulin_4931_unipro_funcoupA%>%
                                                pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                             names_to = "link_source",
                                                             values_to = "uniprot_genes")%>%
                                                distinct(uniprot_genes,.keep_all = T),
                                              by = c("UNIPROT" = "uniprot_genes"))

# 1522 (Endocrine resistance) in overlapp 8 and 9 #####
endocrine_1522 = names(mm39_pathways$descriptions)[mm39_pathways$descriptions == "Endocrine resistance"]
endocrine_1522_id=mm39_pathways$gene_sets[[endocrine_1522]]
endocrine_1522_id=sub("mmu:","",endocrine_1522_id)
endocrine_1522_genes=bitr(endocrine_1522_id,
                        fromType = "ENTREZID",
                        toType = "UNIPROT",
                        OrgDb = org.Mm.eg.db)
endocrine_1522_unipro_funcoupA=left_join(endocrine_1522_genes,funcoup_uniprotA,by = "UNIPROT")
#which genes are affected in factor 8
links_endocrine_1522_f8=uniprot_F8%>%inner_join(endocrine_1522_unipro_funcoupA%>%
                                                pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                             names_to = "link_source",
                                                             values_to = "uniprot_genes")%>%
                                                distinct(uniprot_genes,.keep_all = T),
                                              by = c("UNIPROT" = "uniprot_genes"))
#which genes are affected in factor 9
links_endocrine_1522_f9=uniprot_F9%>%inner_join(endocrine_1522_unipro_funcoupA%>%
                                                pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                             names_to = "link_source",
                                                             values_to = "uniprot_genes")%>%
                                                distinct(uniprot_genes,.keep_all = T),
                                              by = c("UNIPROT" = "uniprot_genes"))


# 4919 (Thyroid hormone signaling pathway) only in factor 10 #####
insulin_4919 = names(mm39_pathways$descriptions)[mm39_pathways$descriptions == "Thyroid hormone signaling pathway"]
insulin_4919_id=mm39_pathways$gene_sets[[insulin_4919]]
insulin_4919_id=sub("mmu:","",insulin_4919_id)
insulin_4919_genes=bitr(insulin_4919_id,
                        fromType = "ENTREZID",
                        toType = "UNIPROT",
                        OrgDb = org.Mm.eg.db)


# 4930 (Type II diabetes mellitus) supposedly in fator 8 ####
diabetes= names(mm39_pathways$descriptions)[mm39_pathways$descriptions == "Type II diabetes mellitus"]
diabetes_id=mm39_pathways$gene_sets[[diabetes]]
diabetes_id=sub("mmu:","",diabetes_id)
diabetes_genes=bitr(diabetes_id,
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Mm.eg.db)

# diabetes_path=pathway(diabetes)
# to print the actual pathway
# pathway(diabetes) |>
#   activate(nodes) |>
#   mutate(convert_mmu=convert_id("mmu"),
#          convert_map=convert_id("pathway")) |>
#   ggraph(x=x, y=y)+
#   geom_edge_parallel(arrow = arrow(length = unit(1,"mm")),
#                      aes(linetype=subtype_name),
#                      end_cap=circle(7.5,"mm"))+
#   geom_node_rect(aes(filter = type=="gene",
#                      fill = I(bgcolor)),
#                  color="black")+
#   geom_node_text(aes(label = convert_mmu),
#                  size=5,family="serif")+
#   theme_void()

diabetes_genes=bitr(diabetes_id,
                    fromType = "ENTREZID",
                    toType = "UNIPROT",
                    OrgDb = org.Mm.eg.db)

funcoup_uniprotA=data.frame(UNIPROT = funcoup_compact$`0:ProteinA`,
                            other_UNIPROT = funcoup_compact$`1:ProteinB`)

diabetes_unipro_funcoupA=left_join(diabetes_genes,funcoup_uniprotA,by = "UNIPROT")
links_diabetes_f8=uniprot_F8%>%inner_join(diabetes_unipro_funcoupA%>%
                                            pivot_longer(cols = c(UNIPROT,other_UNIPROT),
                                                                                  names_to = "link_source",
                                                                                  values_to = "uniprot_genes")%>%
                                            distinct(uniprot_genes,.keep_all = T),
                                          by = c("UNIPROT" = "uniprot_genes"))


#############################################
#this is for doing it individually ######
# get the metabolic pathways ###
# collections=msigdbr_collections()
# print(collections,n = 40)
# metabolic_pathways=msigdbr(species = "Mus musculus",category ="C2",subcategory = "CP:KEGG")%>%
#   distinct(gs_name,gene_symbol,.keep_all = T)
# #the gene symbol is different from ncbi to uniprot that is found in funcoup, so we need to liftover the ncbi gene names to uniprot
# mm_gene_symbol=metabolic_pathways$entrez_gene
# mapped_symbols_to_unipro=select(org.Mm.eg.db,keys = as.character(mm_gene_symbol),columns = c("ENTREZID","UNIPROT","SYMBOL"),keytype = "ENTREZID")
# mapped_symbols_to_unipro$gene_symbol=mapped_symbols_to_unipro$SYMBOL
# metabolic_pathways=metabolic_pathways%>%left_join(mapped_symbols_to_unipro,by="gene_symbol")
# uniprot_names_pathways=metabolic_pathways[,c("UNIPROT","gs_name")]%>%
#   distinct(UNIPROT,gs_name,.keep_all = T)
# 
# #network es funcoup columnas 1 y 2, pathways es kegg los q te interesen
# anubix_links=anubix_links(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                           pathways= as.data.frame(uniprot_names_pathways))
# 
# 
# clustering_factor_8=anubix_clustering(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                                       pathways= as.data.frame(uniprot_names_pathways),
#                                       links_matrix = anubix_links, 
#                                       genesets = uniprot_F8[,c("UNIPROT","factor")])
#   
# clustering_factor_9=anubix_clustering(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                                       pathways= as.data.frame(uniprot_names_pathways),
#                                       links_matrix = anubix_links, 
#                                       genesets = uniprot_F9[,c("UNIPROT","factor")])
# 
# clustering_factor_10=anubix_clustering(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                                       pathways= as.data.frame(uniprot_names_pathways),
#                                       links_matrix = anubix_links, 
#                                       genesets = uniprot_F10[,c("UNIPROT","factor")])
# 
#   
# transivity_factor_8=anubix_transitivity(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                                         pathways= as.data.frame(uniprot_names_pathways),
#                                         links_matrix = anubix_links, 
#                                         genesets = uniprot_F8[,c("UNIPROT","factor")])
# transivity_factor_8=transivity_factor_8%>%filter(`q-value`<0.05)
# anubix_F8=anubix(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                  pathways= as.data.frame(uniprot_names_pathways),
#                  links_matrix = anubix_links, 
#                  genesets = uniprot_F8[,c("UNIPROT","factor")])
# anubix_F8=anubix_F8%>%filter(`q-value`<0.05)
# 
# transivity_factor_9=anubix_transitivity(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                                         pathways= as.data.frame(uniprot_names_pathways),
#                                         links_matrix = anubix_links, 
#                                         genesets = uniprot_F9[,c("UNIPROT","factor")])
# transivity_factor_9=transivity_factor_9%>%filter(`q-value`<0.05)
# 
# anubix_F9=anubix(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                  pathways= as.data.frame(uniprot_names_pathways),
#                  links_matrix = anubix_links, 
#                  genesets = uniprot_F9[,c("UNIPROT","factor")])
# anubix_F9=anubix_F9%>%filter(`q-value`<0.05)
# 
# transivity_factor_10=anubix_transitivity(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                                         pathways= as.data.frame(uniprot_names_pathways),
#                                         links_matrix = anubix_links, 
#                                         genesets = uniprot_F10[,c("UNIPROT","factor")])
# transivity_factor_10=transivity_factor_10%>%filter(`q-value`<0.05)
# 
# anubix_F10=anubix(network= as.data.frame(funcoup_compact[,c(1:2,6)]), 
#                  pathways= as.data.frame(uniprot_names_pathways),
#                  links_matrix = anubix_links, 
#                  genesets = uniprot_F10[,c("UNIPROT","factor")])
# anubix_F10=anubix_F10%>%filter(`q-value`<0.05)



  
  
  
  
  
  
  
  
  
########################################################
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/first_c57_no_groups_mefisto_object.rda")
# relational_table_C571=relational_table_C57%>%
#  mutate(
#    generation=case_when(
#      generation=="F0"~0,
#      generation=="F1"~1,
#      generation=="F2"~2
#    ),
#    sex=case_when(
#      sex=="male"~0,
#      sex=="female"~1
#    ),
#    group=case_when(
#      group=="C"~0,
#      group=="SL"~1
#    ))
# relational_table_C571=relational_table_C571%>%mutate(family=as.numeric(family))
# long_relational_table_C57=pivot_longer(relational_table_C571,cols = c(3:6),values_to = "value", names_to = "covariate")
# long_relational_table_C57$sample=long_relational_table_C57$sample_name
# long_relational_table_C57=long_relational_table_C57[,c("sample","covariate","value")]
# c57_mefisto=set_covariates(c57_mefisto,covariates = as.data.frame(long_relational_table_C57))
# save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_no_groups_with_covariates.rda")
#load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_sex_groups_with_covariates.rda")
#now we need to set the options
# options_data=get_default_data_options(c57_mefisto)
# # gaussian if you are using normalized counts, poisson for raw counts
# opts_model=get_default_model_options(c57_mefisto)
# # opts_model$likelihoods[]="poisson"
# 
# #number of components, or factors as they call them here to calculate
# opts_model$num_factors=10
# opts_training=get_default_training_options(c57_mefisto)
# opts_training$maxiter=1000
# opts_training$convergence_mode="slow"
# opts_training$seed=123
# opts_mefisto=get_default_mefisto_options(c57_mefisto)
# opts_mefisto$model_groups=T
# 
# c57_mefisto=prepare_mofa(c57_mefisto,model_options = opts_model,
#                          mefisto_options = opts_mefisto,
#                          training_options = opts_training,
#                          data_options = options_data)
# 
# c57_mefisto=run_mofa(c57_mefisto,
#                      outfile = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mofa_no_groups_run_ouput.hdf5",
#                      use_basilisk = T)
# save(c57_mefisto,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_no_groups_ran_mofa.rda")

# 
# #downstream analysis, this is the analysis without any group definition ####
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/c57_mefisto_no_groups_ran_mofa.rda")
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_visualization_no_grouping.tiff", height = 2000, width = 2500, res = 150)
# plot_variance_explained(c57_mefisto)
# dev.off()
# #get the actual values
# var_factors_no_grouping=get_variance_explained(c57_mefisto)
# 
# plot_factor_cor(c57_mefisto)
# #check if the components variation also varies during the different generations
# get_scales(c57_mefisto)
# #see how the factors evolve during the generations
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_evolution_no_grouping.tiff", height = 2000, width = 2500, res = 150)
# plot_factors_vs_cov(c57_mefisto,covariates = "generation")
# dev.off()
# 
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_no_grouping_sexseparation_color_sex.tiff", height = 2000, width = 2500, res = 150)
# plot_factors(c57_mefisto,factors = c(2,3,4),color_by = "sex",color_name = "Sex")
# dev.off()
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/mefisto/factors_no_grouping_treatmentseparation_color_sex.tiff", height = 2000, width = 2500, res = 150)
# plot_factors(c57_mefisto,factors = 1:10,color_by = "group_scaled",color_name = "Treatment")
# dev.off()
# 
# #extract the coordinates of the top coordinates that contributes to the factors
# plot_top_weights(c57_mefisto,factors = 1,view = 1,nfeatures = 50)








