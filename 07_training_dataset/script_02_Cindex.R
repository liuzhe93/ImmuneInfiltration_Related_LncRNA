setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/07_training_dataset")
rm(list=ls())

###############################################GSE C-index calculation##########################################
com_auc<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/05_Gaussian_mixture_mode/cluster_results.csv",
                  header=T)
com_auc_sorted=com_auc[order(com_auc$auc,decreasing = T),]
genes<-com_auc_sorted$gene[1]
gene_list<-strsplit(genes,split="\\ \\+\\ ")[[1]]

uniSigExp<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/04_univariateCox_regression/uniSigExp.txt",
                      sep="\t",header=T)
genes<-c("CYB561D2", "LINC00638", "DANCR")
uniSigExp_filtered<-uniSigExp[,c("id","futime","fustat",genes)]
dim(uniSigExp_filtered)
#[1] 566   6
mydata <- na.omit(uniSigExp_filtered) 
dim(mydata)
#[1] 562   14
row.names(mydata)<-mydata$id
mydata<-mydata[,-1]
dim(mydata)
#[1] 562  5

library("survival")
fit <- coxph(Surv(futime, fustat) ~ CYB561D2 + LINC00638 + DANCR, data = mydata)
sum.surv <- summary(fit)
c_index <- sum.surv$concordance
c_index
#         C      se(C) 
#0.59942799 0.02182945 

###############################################TCGA C-index calculation##########################################

