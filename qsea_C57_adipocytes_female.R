#BiocManager::install("qsea", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(qsea, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(BiocManager)
library(BSgenome) #we have the BSgenome.Mmusculus.UCSC.mm39
#install("BSgenome.Mmusculus.UCSC.mm39", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(BSgenome.Mmusculus.UCSC.mm39, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(GenomicRanges)
library(data.table)
library(dplyr)
setwd("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results")

#load the data
relational_table_C57=fread("qsea_with_file_path.txt")
relational_table_C57[,sex := tolower(sex)]
relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)

#create a qsea object with all metadata
# C57_adipo=createQseaSet(sampleTable = relational_table_C57,BSgenome = "BSgenome.Mmusculus.UCSC.mm39",window_size = 300)
# # #add the coverage for each of the 300bps windows we divided the genome
# C57_adipo=addCoverage(C57_adipo,uniquePos = T,paired = T)
# save(C57_adipo,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_bams_loaded.rda")
# # load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_bams_loaded.rda")
# #control for CNV. as we don't have WGS data it going to be infered from the MeDIP data
# C57_adipo=addCNV(C57_adipo,file_name = "file_name",paired = T,MeDIP = T)
# save(C57_adipo,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_cnvs_loaded.rda")
# #normalizing the libraries with the TMM approach
# C57_adipo=addLibraryFactors(C57_adipo)
# #estimation of CpG density
# C57_adipo=addPatternDensity(C57_adipo,"CG",name = "CpG")
# C57_adipo=addOffset(C57_adipo,enrichmentPattern = "CpG")
# save(C57_adipo,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized.rda")
# #enrichment efficiency, this can be supplied with validated data or to try and estimate from the data
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized.rda")
# ee=which(getRegions(C57_adipo)$CpG_density>1 & getRegions(C57_adipo)$CpG_density<15)
# signal=(15-getRegions(C57_adipo)$CpG_density[ee])*0.55/15+0.25
# set_blind=addEnrichmentParameters(C57_adipo,enrichmentPattern = "CpG",windowIdx = ee, signal = signal)
# save(set_blind,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_added.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_added.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized.rda")
# #get the normalized count matrix
# norm_counts_c57=qsea:::getNormalizedValues(C57_adipo,methods = normMethod("rpm"))
# regs_c57=getRegions(C57_adipo)
# coords=data.frame(chr=as.character(seqnames(regs_c57)),
#                   start=start(regs_c57),
#                   end=end(regs_c57),
#                   strand=as.character(strand(regs_c57)))
# rpm_counts_C57=cbind(coords,as.data.frame(norm_counts_c57))
# save(rpm_counts_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/rpm_counts.rda")

#fraction_noise=summary(getOffset(C57_adipo,scale="fraction"))
#write.table(fraction_noise,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/fraction_of_noise_MeDIP_qsea.txt")
#plotEPmatrix(set_blind) #don't run yourself
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/CNV_plot_C57_adipocites.tiff", height = 2000, width = 2500, res = 150)
# plotCNV(set_blind)
# dev.off()
pca_subset=getPCA(set_blind,norm_method = "rpm")
save(pca_subset,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_all_ind_object.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_object.rda")
# variance_PC=pca_subset@svd
# prop_var=variance_PC$d^2/sum(variance_PC$d^2)
# percen_var=100*prop_var
# groups_C57=getSampleTable(C57_adipo)
# groups_C57$color=ifelse(groups_C57$group == "C", "#0072B2", ifelse(groups_C57$group == "SL","#D55E00", NA))
# #do also colors for the sex, IT WAS SEPARATED BY SEX
# groups_C57$color_sex=ifelse(groups_C57$sex == "male", "#009E73", ifelse(groups_C57$sex == "female","#CC79A7", NA))
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_1_2_BY_SEX_plot.tiff",
#      height = 1200, width = 2500, res = 150)
# plotPCA(pca_subset,plotComponents=c(1,2),bg=groups_C57$color_sex,main="C57 adipocites")
# dev.off()

design_C57=model.matrix(~group*sex*generation*family,getSampleTable(C57_adipo))
glm_c57=fitNBglm(set_blind,design_C57,norm_method = "rpm") #you need to run this with the object from addEnrichmentParameters()
save(glm_c57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_C57_all_ind.rda")
glm_c57=addContrast(set_blind,glm_c57,coef = 2, name = "treatment")
glm_c57=addContrast(set_blind,glm_c57,coef = 3, name = "sex")
glm_c57=addContrast(set_blind,glm_c57,coef = 4, name = "interaction_treatment_and_sex")

save(glm_c57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_C57_with_contrast.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_C57_with_contrast.rda")
significant_sex=isSignificant(glm_c57,contrast = "sex",fdr_th = .05)
significant_treatment=isSignificant(glm_c57,contrast = "treatment",fdr_th = .05)
significant_interaction=isSignificant(glm_c57,contrast = "interaction_treatment_and_sex",fdr_th = .05)
significant_peaks_C57=makeTable(set_blind,glm = glm_c57,groupMeans = getSampleGroups(set_blind),keep = significant_sex, norm_methods = "beta")
save(significant_peaks_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/significant_peaks_c57.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/significant_peaks_c57.rda")
# write.table(significant_peaks_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/significant_peaks_c57.txt",sep = "\t",quote = F,row.names = F)

#let's exclude from the analysis the regions that are associated with sex differences ####
#we need to know the indices of the sex DMRs, and we need to add the generation and family variable
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_added.rda")
# 
# design_C57=model.matrix(~group*sex,getSampleTable(C57_adipo))
# glm_c57=fitNBglm(set_blind,design_C57,norm_method = "rpm")





#####################################################################################################
#now as we have clear very different signals from females and males, let's fit again the model but using only the females or the males #####
# first subset females
females_C57_relational_table=relational_table_C57%>%filter(sex=="female")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_added.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized.rda")
#we need to subset CNV, sampleTable, libraries and enrichment from the set_blind object
C57_female=C57_adipo
keep=intersect(females_C57_relational_table$sample_name,getSampleNames(C57_adipo))
C57_female@count_matrix=getCounts(C57_adipo,samples=keep)
st=as.data.frame(getSampleTable(C57_adipo))
st=st[match(keep,st$sample_name), , drop=F]
rownames(st)=st$sample_name
C57_female@sampleTable=st
if(!is.null(C57_female@libraries)&&length(C57_female@libraries)>0){
  if(!is.null(rownames(C57_female@libraries))&&all(keep%in%rownames(C57_female@libraries))){
    C57_female@libraries=C57_female@libraries[keep, , drop=F]
  } else if(!is.null(colnames(C57_female@libraries))&&all(keep%in%colnames(C57_female@libraries))){
    C57_female@libraries=C57_female@libraries[, keep, drop=F]
  }
}
zyg=getZygosity(C57_adipo)
C57_female@zygosity=zyg[keep, , drop=F]
keep=getSampleNames(C57_female)
if(hasCNV(C57_adipo)){
  cnv_f=getCNV(C57_adipo,samples=keep)
  C57_female=qsea:::setCNV(C57_female,cnv_f)
}

##################################################################################################
#this is to take the enrichment pattern of both males and females ###
# enrichment_c57_female=set_blind
# keep=intersect(females_C57_relational_table$sample_name,getSampleNames(enrichment_c57_female))
# enrichment_c57_female@count_matrix=getCounts(set_blind,samples=keep)
# st=as.data.frame(getSampleTable(set_blind))
# st=st[match(keep,st$sample_name), , drop=F]
# rownames(st)=st$sample_name
# enrichment_c57_female@sampleTable=st
# if(!is.null(enrichment_c57_female@libraries)&&length(enrichment_c57_female@libraries)>0){
#   if(!is.null(rownames(enrichment_c57_female@libraries))&&all(keep%in%rownames(enrichment_c57_female@libraries))){
#     enrichment_c57_female@libraries=enrichment_c57_female@libraries[keep, , drop=F]
#   } else if(!is.null(colnames(enrichment_c57_female@libraries))&&all(keep%in%colnames(enrichment_c57_female@libraries))){
#     enrichment_c57_female@libraries=enrichment_c57_female@libraries[, keep, drop=F]
#   }
# }
# zyg=getZygosity(set_blind)
# enrichment_c57_female@zygosity=zyg[keep, , drop=F]
# keep=getSampleNames(enrichment_c57_female)
# if(hasCNV(set_blind)){
#   cnv_f=getCNV(set_blind,samples=keep)
#   enrichment_c57_female=qsea:::setCNV(enrichment_c57_female,cnv_f)
# }
# if(qsea:::hasEnrichment(set_blind)){
#   enr=set_blind@enrichment
#   if(!is.null(enr$parameters)) enr$parameters =enr$parameters[keep, , drop=F]
#   if(!is.null(enr$factors)) enr$factors =enr$factors[, keep, drop=F]
#   
#   enrichment_c57_female=qsea:::setEnrichment(enrichment_c57_female,enr)
# }
#######################################################################
#recalculate the enrichment pattern in females ####
ee=which(getRegions(C57_female)$CpG_density>1 & getRegions(C57_female)$CpG_density<15)
signal=(15-getRegions(C57_female)$CpG_density[ee])*0.55/15+0.25
enrichment_c57_female=addEnrichmentParameters(C57_female,enrichmentPattern = "CpG",windowIdx = ee, signal = signal)
save(enrichment_c57_female,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_females.rda")

pca_female_C57=getPCA(enrichment_c57_female,norm_method = "rpm")
save(pca_female_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_female_object.rda")
# groups_female_C57=getSampleTable(C57_female)
# groups_female_C57$group=ifelse(groups_female_C57$group == "C", "#0072B2", ifelse(groups_female_C57$group == "SL","#D55E00", NA))
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_1_2_only_females_by_treatment_plot.tiff",
#      height = 1200, width = 2500, res = 150)
# plotPCA(pca_female_C57,plotComponents=c(1,2),bg=groups_female_C57$group,main="Female C57 adypocites")
# dev.off()

design_female_C57=model.matrix(~group*generation*family,getSampleTable(C57_female))
glm_female_c57=fitNBglm(enrichment_c57_female,design_female_C57,norm_method = "rpm") #you need to run this with the object from addEnrichmentParameters()
save(glm_female_c57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_C57_females.rda")
glm_female_c57=addContrast(enrichment_c57_female,glm_female_c57,coef = 2, name = "treatment")
glm_female_c57=addContrast(enrichment_c57_female,glm_female_c57,coef = 3, name = "generation")
glm_female_c57=addContrast(enrichment_c57_female,glm_female_c57,coef = 4, name = "interaction_treatment_and_generation")
save(glm_female_c57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_females_C57_with_contrast.rda")

# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_females_C57_with_contrast.rda")
# significant_female_generation=isSignificant(glm_female_c57,contrast = "generation",fdr_th = .1)
# significant_female_treatment=isSignificant(glm_female_c57,contrast = "treatment",fdr_th = .1)
# significant_female_interaction=isSignificant(glm_female_c57,contrast = "interaction_treatment_and_generation",fdr_th = .1)

# significant_peaks_C57=makeTable(enrichment_c57_female,glm = glm_female_c57,
#                                 groupMeans = getSampleGroups(enrichment_c57_female),keep = XXXXX, norm_methods = "beta")

####################################################################################################################
#In the males we need to delete sample /proj/naiss2024-23-57/C57_female_lineage_adipocytes/aligned/27_unique_sorted.bam
males_C57_relational_table=relational_table_C57%>%filter(sex=="male")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_added.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized.rda")
#we need to subset CNV, sampleTable, libraries and enrichment from the set_blind object
C57_male=C57_adipo
keep=intersect(males_C57_relational_table$sample_name,getSampleNames(C57_adipo))
C57_male@count_matrix=getCounts(C57_adipo,samples=keep)
st=as.data.frame(getSampleTable(C57_adipo))
st=st[match(keep,st$sample_name), , drop=F]
rownames(st)=st$sample_name
C57_male@sampleTable=st
if(!is.null(C57_male@libraries)&&length(C57_male@libraries)>0){
  if(!is.null(rownames(C57_male@libraries))&&all(keep%in%rownames(C57_male@libraries))){
    C57_male@libraries=C57_male@libraries[keep, , drop=F]
  } else if(!is.null(colnames(C57_male@libraries))&&all(keep%in%colnames(C57_male@libraries))){
    C57_male@libraries=C57_male@libraries[, keep, drop=F]
  }
}
zyg=getZygosity(C57_adipo)
C57_male@zygosity=zyg[keep, , drop=F]
keep=getSampleNames(C57_male)
if(hasCNV(C57_adipo)){
  cnv_f=getCNV(C57_adipo,samples=keep)
  C57_male=qsea:::setCNV(C57_male,cnv_f)
}

##########################################################################################
# enrichment_c57_male=set_blind
# keep=intersect(males_C57_relational_table$sample_name,getSampleNames(enrichment_c57_male))
# enrichment_c57_male@count_matrix=getCounts(set_blind,samples=keep)
# st=as.data.frame(getSampleTable(set_blind))
# st=st[match(keep,st$sample_name), , drop=F]
# rownames(st)=st$sample_name
# enrichment_c57_male@sampleTable=st
# if(!is.null(enrichment_c57_male@libraries)&&length(enrichment_c57_male@libraries)>0){
#   if(!is.null(rownames(enrichment_c57_male@libraries))&&all(keep%in%rownames(enrichment_c57_male@libraries))){
#     enrichment_c57_male@libraries=enrichment_c57_male@libraries[keep, , drop=F]
#   } else if(!is.null(colnames(enrichment_c57_male@libraries))&&all(keep%in%colnames(enrichment_c57_male@libraries))){
#     enrichment_c57_male@libraries=enrichment_c57_male@libraries[, keep, drop=F]
#   }
# }
# zyg=getZygosity(set_blind)
# enrichment_c57_male@zygosity=zyg[keep, , drop=F]
# keep=getSampleNames(enrichment_c57_male)
# if(hasCNV(set_blind)){
#   cnv_f=getCNV(set_blind,samples=keep)
#   enrichment_c57_male=qsea:::setCNV(enrichment_c57_male,cnv_f)
# }
# if(qsea:::hasEnrichment(set_blind)){
#   enr=set_blind@enrichment
#   if(!is.null(enr$parameters)) enr$parameters =enr$parameters[keep, , drop=F]
#   if(!is.null(enr$factors)) enr$factors =enr$factors[, keep, drop=F]
#   
#   enrichment_c57_male=qsea:::setEnrichment(enrichment_c57_male,enr)
# }
###############################################################################################
#recalculate the enrichment pattern in males ####
ee=which(getRegions(C57_male)$CpG_density>1 & getRegions(C57_male)$CpG_density<15)
signal=(15-getRegions(C57_male)$CpG_density[ee])*0.55/15+0.25
enrichment_c57_male=addEnrichmentParameters(C57_male,enrichmentPattern = "CpG",windowIdx = ee, signal = signal)
save(enrichment_c57_male,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_males.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_enrichment_males.rda")

pca_male_C57=getPCA(enrichment_c57_male,norm_method = "rpm")
save(pca_male_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_male_object.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_male_object.rda")
# groups_male_C57=getSampleTable(C57_male)
# groups_male_C57$group=ifelse(groups_male_C57$group == "C", "#0072B2", ifelse(groups_male_C57$group == "SL","#D55E00", NA))
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_1_2_only_males_by_treatment_plot.tiff",
#      height = 1200, width = 2500, res = 150)
# plotPCA(pca_male_C57,plotComponents=c(1,2),bg=groups_male_C57$group,main="male C57 adypocites")
# dev.off()

design_male_C57=model.matrix(~group*generation*family,getSampleTable(C57_male))
glm_male_c57=fitNBglm(enrichment_c57_male,design_male_C57,norm_method = "rpm") #you need to run this with the object from addEnrichmentParameters()
save(glm_male_c57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_C57_males.rda")
glm_male_c57=addContrast(enrichment_c57_male,glm_male_c57,coef = 2, name = "treatment")
glm_male_c57=addContrast(enrichment_c57_male,glm_male_c57,coef = 3, name = "generation")
glm_male_c57=addContrast(enrichment_c57_male,glm_male_c57,coef = 4, name = "interaction_treatment_and_generation")
save(glm_male_c57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_males_C57_with_contrast.rda")

# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/glm_males_C57_with_contrast.rda")
# significant_male_generation=isSignificant(glm_male_c57,contrast = "generation",fdr_th = .1)
# significant_male_treatment=isSignificant(glm_male_c57,contrast = "treatment",fdr_th = .1)
# significant_male_interaction=isSignificant(glm_male_c57,contrast = "interaction_treatment_and_generation",fdr_th = .1)
# significant_peaks_male_C57=makeTable(enrichment_c57_male,glm = glm_male_c57,
#                                 groupMeans = getSampleGroups(enrichment_c57_male),keep = significant_male_interaction, norm_methods = "beta")



###############################
# enr=qsea:::estimateEnrichmentLM(C57_adipo,windowIdx = ee,signal = signal,pattern_name = "CpG")
# apply(enr$factors,2, function(x) sum(is.finite(x)))
# upperQ=apply(enr$factors,2,quantile,p=0.75,na.rm=T)
# which(!is.finite(upperQ))
# bad_samples=c("2_6_2","2_4","4_2","11_3")
# raw_bad=getCounts(C57_adipo,windows=ee,samples=bad_samples)
# cbind(
#   sumNA=colSums(is.na(raw_bad)),
#   sumCounts=colSums(raw_bad,na.rm = T),
#   nFinite=colSums(is.finite(raw_bad))
# )
# nm=normMethod("beta")
# nm$beta=nm$beta[!grepl("enrichment",nm$beta)]
# vals_bad=qsea:::getNormalizedValues(C57_adipo,methods = nm,windows = ee,samples = bad_samples)
# 
# getOffset(C57_adipo,bad_samples,scale="fraction")
# getLibSize(C57_adipo,bad_samples,normalized=F)
# getCounts(C57_adipo,samples=bad_samples)
# library(Rsamtools)
# st=getSampleTable(C57_adipo)
# file.exists(st[bad_samples, "file_name"])
# bam=as.character(st[bad_samples, "file_name",drop=T])
# cb=lapply(bam,function(f) countBam(f))

#so the bams were truncated, I will rerun only those but in the mean time I will drop those and keep the script building
# C57_clean=C57_adipo
# keep=setdiff(getSampleNames(C57_adipo),bad_samples)
# C57_clean@count_matrix=getCounts(C57_adipo,samples=keep)
# st=as.data.frame(getSampleTable(C57_adipo))
# st=st[match(keep,st$sample_name), , drop=F]
# rownames(st)=st$sample_name
# C57_clean@sampleTable=st
# if(!is.null(C57_clean@libraries)&&length(C57_clean@libraries)>0){
#   if(!is.null(rownames(C57_clean@libraries))&&all(keep%in%rownames(C57_clean@libraries))){
#     C57_clean@libraries=C57_clean@libraries[keep, , drop=F]
#   } else if(!is.null(colnames(C57_clean@libraries))&&all(keep%in%colnames(C57_clean@libraries))){
#     C57_clean@libraries=C57_clean@libraries[, keep, drop=F]
#   }
# }
# zyg=getZygosity(C57_adipo)
# C57_clean@zygosity=zyg[keep, , drop=F]
# C57_clean@enrichment=list()
# C57_clean@cnv=GenomicRanges::GRanges()

# C57_clean=addCNV(C57_clean,file_name = "file_name",paired = T,MeDIP = T)
# save(C57_clean,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/cnv_normalized_wihtout_problematic_samples.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/cnv_normalized_wihtout_problematic_samples.rda")
# nrow(C57_clean@sampleTable)
# kp_69=getSampleNames(C57_clean)
# length(kp_69)
# nrow(getSampleTable(C57_clean))
# length(getSampleNames(C57_clean))
# head(C57_clean@libraries[["file_name"]])
# lib=C57_clean@libraries[["file_name"]]
# samples=getSampleNames(C57_clean)
# C57_clean@libraries[["file_name"]]=lib[samples, , drop=F]
# C57_clean=addLibraryFactors(C57_clean)
# #estimation of CpG density
# C57_clean=addPatternDensity(C57_clean,"CG",name = "CpG")
# C57_clean=addOffset(C57_clean,enrichmentPattern = "CpG")
# save(C57_clean,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized_wihtout_problematic_samples.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_normalized_wihtout_problematic_samples.rda")
# #enrichment efficiency, this can be supplied with validated data or to try and estimate from the data
# ee=which(getRegions(C57_clean)$CpG_density>1 & getRegions(C57_clean)$CpG_density<15)
# signal=(15-getRegions(C57_clean)$CpG_density[ee])*0.55/15+0.25
# set_blind=addEnrichmentParameters(C57_clean,enrichmentPattern = "CpG",windowIdx = ee, signal = signal)
# 
# #from here is new code, copy from here
# summary(getOffset(C57_clean,scale="fraction"))
# #plotEPmatrix(set_blind) #don't run yourself
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/subset_CNV_plot.tiff", 
#      height = 1200, width = 2500, res = 150)
# plotCNV(set_blind)
# dev.off()
# # pca_subset=getPCA(set_blind,norm_method = "beta")
# # save(pca_subset,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_object.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/PCA_object.rda")
# variance_PC=pca_subset@svd
# prop_var=variance_PC$d^2/sum(variance_PC$d^2)
# percen_var=100*prop_var
# groups_C57=getSampleTable(C57_clean)
# groups_C57$color=ifelse(groups_C57$group == "C", "#0072B2", ifelse(groups_C57$group == "SL","#D55E00", NA))
# #do also colors for the sex, IT WAS SEPARATED BY SEX
# groups_C57$color_sex=ifelse(groups_C57$sex == "male", "#009E73", ifelse(groups_C57$sex == "female","#CC79A7", NA))
# tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/subset_PCA_1_2_BY_SEX_plot.tiff", 
#      height = 1200, width = 2500, res = 150)
# pc=c(1,2)
# plotPCA(pca_subset,plotComponents=c(1,2),bg=groups_C57$color_sex,main="C57 adipocites")
# dev.off()
# 
# design_C57=model.matrix(~group*sex,getSampleTable(C57_clean))
# glm_c57=fitNBglm(set_blind,design_C57,norm_method = "beta") #you need to run this with the object from addEnrichmentParameters()
# glm_c57=addContrast(C57_clean,glm_c57,coef = c(2,3),name = "treatment")



################################


















