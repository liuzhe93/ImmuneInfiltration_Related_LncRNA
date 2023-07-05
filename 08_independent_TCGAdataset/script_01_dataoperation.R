setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
rm(list=ls())
#安装TCGAbiolinks
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("TCGAbiolinks")
#######################################download data######################################

library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)
TCGAbiolinks::getGDCprojects()$project_id
cancer_type="TCGA-COAD"
clinical<-GDCquery_clinic(project=cancer_type,type="clinical")
dim(clinical)
# 461  70
View(clinical)
save(clinical,file="COAD_clinical.Rdata")
write.csv(clinical, file="TCGAbiolinks-COAD-clinical.csv")
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
query <- GDCquery(project = cancer_type,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = "STAR - Counts")
samplesDown <- getResults(query,cols=c("cases"))  

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
write.csv(dataSmNT,"Normal_barcode.csv",quote=F,row.names=F)
write.csv(dataSmTP,"Cancer_barcode.csv",quote=F,row.names=F)
queryDown <- GDCquery(project = "TCGA-COAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts", 
                      barcode = dataSmTP)
######################################GDCdownload() and GDCquery（）#######################
setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset/download_data")

GDCdownload(queryDown,method = "api", directory = "GDCdata",files.per.chunk = 20)

GDCdownload(query = queryDown)

########create SE（SummarizedExperiment）file######
dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename =
                          "COAD_case.rda")
########TCGAanalyze_Preprocessing()preprocess clean data, remove outliers##########
testData<-load("COAD_case.rda")
testData
data
assays(data)
names(assays(data))[5]
dataPrep2 <- TCGAanalyze_Preprocessing(object = data,
                                       cor.cut = 0.6,
                                       datatype = "fpkm_unstrand")
#Number of outliers: 2
#将预处理后的数据dataPrep2，写入新文件“COAD_dataPrep.csv”
write.csv(dataPrep2,file = "COAD_dataPrep.csv",quote = FALSE)
#########################TCGAtumor_purity（）########################

purityDATA <- TCGAtumor_purity(colnames(data), 0, 0, 0, 0, 0.6)

Purity_COAD<-purityDATA$pure_barcodes
dataPrep2<-as.data.frame(dataPrep2)
overlap<-intersect(colnames(dataPrep2), Purity_COAD)
puried_data<-subset(dataPrep2, select = overlap)
###################################gene annotation of expression matrix################################################
library("SummarizedExperiment")
anno_data<-rowData(data)   #传入数据dataPrep1必须为SummarizedExperiment对象
anno_data<-as.data.frame(anno_data)
anno_data<-subset(anno_data, select = c("gene_id", "gene_name"))
puried_data$gene_id<-rownames(puried_data)
merged_data<-merge(anno_data, puried_data, by = "gene_id")
write.csv(merged_data,file = "puried.COAD.csv",quote = FALSE, row.names = F)
merged_data_uniq<-aggregate(merged_data[,3:ncol(merged_data)], list(merged_data$gene_name), mean)
dim(merged_data)
#[1] 60660   453
rownames(merged_data_uniq)<-merged_data_uniq$Group.1
merged_data_uniq$Group.1<-NULL
dim(merged_data_uniq)
#[1] 59427   451
write.csv(merged_data_uniq,file = "puried.COAD.uniq.csv",quote = FALSE)


#########################normalization and filter##########################

dataFilt_COAD_final <- read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset/download_data/puried.COAD.uniq.csv",
                                header = T,check.names = FALSE,row.names = 1)


data_sel<-rbind(dataFilt_COAD_final["CYB561D2",], dataFilt_COAD_final["LINC00638",], dataFilt_COAD_final["DANCR",])
data_sel_t<-t(data_sel)
data_sel_t<-as.data.frame(data_sel_t)
data_sel_t$id<-substring(rownames(data_sel_t),1,12)


setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
cli_data<-read.csv("TCGAbiolinks-COAD-clinical.csv",header=T,row.names=1)
cli_data$os<-ifelse(cli_data$vital_status=='Alive',cli_data$days_to_last_follow_up,cli_data$days_to_death)
library(dplyr)
cli_data_sel<-cli_data %>%
  select(submitter_id,os,vital_status,gender,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n,age_at_index)
#0=alive, 1=dead.
cli_data_sel$fustat<-ifelse(cli_data_sel$vital_status=='Alive',0,1)
cli_data_sel$futime<-cli_data_sel$os
results_input<-cli_data_sel %>%
  select(submitter_id,futime,fustat,age_at_index,gender,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n)
colnames(results_input)<-c("id","futime","fustat","age","gender","AJCC_stage","AJCC_T","AJCC_M","AJCC_N")
write.csv(results_input,"indepInput.csv",quote=F,row.names = F)

exp_cli<-merge(results_input,data_sel_t,by="id")
exp_cli<-subset(exp_cli, select = c("id", "futime", "fustat", "CYB561D2", "LINC00638", "DANCR"))

for(i in 1:nrow(exp_cli)){
  exp_cli$id[i]<-paste(exp_cli$id[i],i,sep = "_")
}
write.csv(exp_cli, "TCGA_exp_surv.csv", row.names = F, quote = F)


#cli_data_sel$Age_ope<-ifelse(cli_data_sel$age_at_index>=60,">=60","<60")

