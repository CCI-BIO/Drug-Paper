# Association analysis of tumor type and drug: 
# The Wilcoxon rank-sum test was adopted to test significant differences in drug response profiles 
# between a tumor type of interest and all other types of tumors,
# followed by the Benjamini-Hochberg (BH) procedure for multiple tests correction.  
# The association pairs with P-value < 0.05 and Z-score difference > 0.5 or < -0.5 
# suggest specific tumor type was significantly sensitive or resistant to the corresponding drug. 

library(xlsx)
library(readxl)

# drugs metadata
FileC <- "2020_InVitroDrugInformationWorksheet_updated_CCI_PMC_2.xlsx"
DrugInfo <- read.xlsx(FileC, "Compounds_MergedOverview")
DrugInfo <- DrugInfo[1:126,]
rownames(DrugInfo) <- DrugInfo[,2]
DrugInfo <- DrugInfo[,-(2:7)]

target_drug<-rownames(DrugInfo)[DrugInfo$Category == "Targeted"]
chemo_drug<-rownames(DrugInfo)[DrugInfo$Category == "Chemotherapeutic"]

# drug efficacy data
auc<-as.data.frame(read_excel("InVitro_DrugResponseDatabase_MASTER_20220812_Complete.xlsx", sheet="AUC")[-c(1:2),])
auc<-auc[match(DrugInfo$DrugID_CCI,auc$`Compound ID`),]
rownames(auc)<-rownames(DrugInfo)

auc<-as.data.frame(t(auc[,-c(1:2)]))
auc <- data.frame(lapply(auc, function(x) as.numeric(as.character(x))),
                  check.names=F, row.names = rownames(auc))
auc_z<-as.data.frame(scale(auc,center = T,scale = T))

# cohort info and Tumor type selection
id<-read.csv("New_cohort_ID_2_9.csv",check.names = F)
id<-id[which(id$RNAbatch =="zero" | is.na(id$RNAbatch )),]

id<-id[match(rownames(auc_z),id$`HTS ID`),]
id$`Tumor type`[which(id$`Tumor type` == "ST")]<-"SAR"
id$`Tumor type`[which(id$`Tumor type` == "KT")]<-"SO"
id$`Tumor type`[which(id$`Tumor type` == "CNS")]<-"BT"

auc_z<-as.data.frame(cbind(id$`Tumor type`,auc_z))
colnames(auc_z)[1]<-"cancer_type"

# Create a table for association analysis 
n=126
association<-data.frame(drug=rep("",n*5),cancer=rep("",n*5),Pvalue=rep("",n*5),`z-score difference`=rep("",n*5),
                        `AvgZ in cancer`=rep("",n*5),`AvgZ in others`=rep("",n*5),
                        `Cancer cases`=rep("",n*5),`Other cancer cases`=rep("",n*5),check.names = F)

# Run association analysis
count=1
for (i in 2:ncol(auc_z)){
  for (j in unique(auc_z$cancer_type)){

    association$drug[count]<-colnames(auc_z)[i]
    association$cancer[count]<-j
    
    cancer.list<-auc_z[[colnames(auc_z)[i]]][auc_z$cancer_type == j]
    others.list<-auc_z[[colnames(auc_z)[i]]][auc_z$cancer_type != j]
    cancer.list<-as.numeric(cancer.list[!is.na(cancer.list)])
    others.list<-as.numeric(others.list[!is.na(others.list)])
    
    w.test<-wilcox.test(cancer.list,others.list, alternative = "two.sided",exact = F)
    
    association$Pvalue[count]<-w.test$p.value
    association$`z-score difference`[count]<-mean(others.list)-mean(cancer.list)
    association$`AvgZ in cancer`[count]<-mean(cancer.list)
    association$`AvgZ in others`[count]<-mean(others.list)
    association$`Cancer cases`[count]<-length(cancer.list)
    association$`Other cancer cases`[count]<-length(others.list)
    count=count+1
  }
}

association$Pvalue<-as.numeric(association$Pvalue)
association$`z-score difference`<-as.numeric(association$`z-score`)

# Multiple tests correction
association$Adj.Pvalue<-""
association$Adj.Pvalue[association$cancer =="BT"]<-p.adjust(association$Pvalue[association$cancer =="BT"],method="BH")
association$Adj.Pvalue[association$cancer =="HM"]<-p.adjust(association$Pvalue[association$cancer =="HM"],method="BH")
association$Adj.Pvalue[association$cancer =="NB"]<-p.adjust(association$Pvalue[association$cancer =="NB"],method="BH")
association$Adj.Pvalue[association$cancer =="SAR"]<-p.adjust(association$Pvalue[association$cancer =="SAR"],method="BH")
association$Adj.Pvalue[association$cancer =="SO"]<-p.adjust(association$Pvalue[association$cancer =="SO"],method="BH")

# Annotate sensitivity
association$sensitivity<-"No"
association$sensitivity[association$`z-score difference`> 0.5 & association$Adj.Pvalue < 0.01 ]<-"Relatively sensitive"
association$sensitivity[association$`z-score difference`< -0.5 & association$Adj.Pvalue < 0.01 ]<-"Relatively resistant"

# label outliers (target drugs only)
association$delabel <- NA
for (i in 1:630){
  if (association$sensitivity[i] != "No" & association$drug[i] %in% target_drug){
    association$delabel[i]<- paste(association$drug[i],association$cancer[i],sep="-")
  }
}