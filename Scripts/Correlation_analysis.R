
# Gene expression correlation analysis with drug efficacy:
# Correlations with gene expression were analyzed for the top 50 drugs with the most varying AUC values 
# for all non-hematological samples within our cohort with matching RNA-seq data based on highest coefficient variations

# Inputs: 
# New_cohort_ID.csv
# InVitro_DrugResponseDatabase.xlsx
# GeneExpression_TPM_Counts.txt

# Outputs:
# CoV_Top50drugs.csv
# r_Gene_Drug.txt

library(xlsx)
library(readxl)
library(tidyverse)
library(dplyr)

## Cohort info ##
id<-read.csv("New_cohort_ID.csv",check.names = F)
id<-id[which(id$RNAbatch =="zero" | is.na(id$RNAbatch )),] 
cohort_solid<-id[id$`Tumor type` !="HM",] 
# cohort with matching WGS and RNAseq data
cohort_solid_rna<-cohort_solid[which(grepl("RNA",cohort_solid$`Molecular profiling`) & cohort_solid$`RNA Dataset` == "ExpressionAnalysis"),]

## Drug efficacy data ##
FileA <- "InVitro_DrugResponseDatabase.xlsx"
results_AUC <- read.xlsx(FileA,"AUC",check.names =F)[3:128,]
rownames(results_AUC)<-results_AUC$`Drug name`
results_AUC<-results_AUC[,-c(1,2)]
AUC <- as.matrix(t(results_AUC))

# AUC Z-score calculated based on solid tumor cohort
AUC<-AUC[rownames(AUC) %in% cohort_solid$`HTS ID`,]
AUC_Z<-scale(AUC,center = T,scale=T)

AUC_Z<-AUC_Z[rownames(AUC_Z) %in% cohort_solid_rna$`HTS ID` ,]
AUC<-AUC[rownames(AUC) %in% cohort_solid_rna$`HTS ID`,]

# Cov of drugs based on AUC
cv_list_drug<-apply(AUC,2, function(x) sd(x,na.rm = T) / mean(x,na.rm = T))
cv_list_drug<-cv_list_drug[order(cv_list_drug,decreasing = T)]
# Top drugs with highest Cov
cv_list_drug_50<-cv_list_drug[1:50]
AUC<-AUC[,colnames(AUC) %in% names(cv_list_drug_50)]
`%notin%` <- Negate(`%in%`)
AUC<-AUC[,colnames(AUC)  %notin% c("Trametinib","Cobimetinib","Selumetinib")]  

AUC_Z<-AUC_Z[,colnames(AUC_Z) %in% colnames(AUC)]
AUC_Z<-AUC_Z[rownames(AUC),] 

## Gene expression TPM ##
tpm<-read.delim("GeneExpression_TPM_Counts.txt",check.names = F)
rownames(tpm)<-tpm$gene_id
tpm<-tpm[,-(1:2)]
tpm<-tpm[,colnames(tpm) %in% cohort_solid_rna$`RNAseq ID`]

# Cov of genes based on TPM
cv_list<-apply(tpm,1, function(x) sd(x) / abs(mean(x)))
cv_list<-cv_list[order(cv_list,decreasing = T)]
# Top 75% most variable genes
tpm<-tpm[rownames(tpm) %in% names(cv_list[1:19750]),] # 0.75 x 26334 = 19750
# Only select genes for which TPM>1 in >=5 samples
tpm_fil<-tpm[rowSums(tpm > 1) >= 5,]
tpm_log<- log(tpm_fil+0.01, base=2)
tpm_log<- t(tpm_log)

## Correlation matrix ## 
AUC_Z<-AUC_Z[cohort_solid_rna$`HTS ID`,]
tpm_log<-tpm_log[cohort_solid_rna$`RNAseq ID`,]

cor_matrix<- stats::cor(AUC_Z,tpm_log,method = 'pearson',use = "pairwise")
cor_matrix <- as.data.frame(cor_matrix)
flat_cor_matrix<-cor_matrix %>%
  rownames_to_column(var = "drug") %>%  
  gather(gene, r, -drug) 


