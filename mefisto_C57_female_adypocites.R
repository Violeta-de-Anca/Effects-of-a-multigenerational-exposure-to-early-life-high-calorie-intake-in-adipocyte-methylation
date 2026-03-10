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








