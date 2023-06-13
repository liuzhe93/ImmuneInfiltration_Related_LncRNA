setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())

merged_data<-read.csv("ExpPro_rowGene_colImmuneCell.csv", header = T, row.names = 1)
top30_genes<-read.csv("top30_expressed_lncRNA_in_each_ImmuCell.csv", header = T, row.names = 1)
idx<-top30_genes$top30_merge_uniq
data_top30<-merged_data[idx,]

rmax <- vector()
rmin <- vector()
for (i in 1:nrow(data_top30)){
  rmax = c(rmax,max(data_top30[i,]))
  rmin = c(rmin,min(data_top30[i,]))
}
temp <- cbind(data_top30,rmax,rmin)

data_norm<-data_top30
for(i in 1:nrow(data_norm)){
  for(j in 1:ncol(data_norm)){
    data_norm[i,j] = (data_norm[i,j]-temp[i,21])/(temp[i,20]-temp[i,21])
  }
}
write.csv(data_norm, "data_norm.csv", quote = F)

n_celltype = ncol(data_norm)
data_spec<-matrix(data=rep(0,n_celltype), nrow = n_celltype, ncol = 1)
data_spec<-as.data.frame(data_spec)
colnames(data_spec)<-"TissSpecScor"
for(i in 1:nrow(data_norm)){
  data_spec[i,1] = 0
  for(j in 1:ncol(data_norm)){
    temp = (1-data_norm[i,j])/(n_celltype-1)
    data_spec[i,1] = temp + data_spec[i,1]
  }
}
rownames(data_spec)<-row.names(data_norm)
write.csv(data_spec, "data_tissueSpecScor.csv", quote = F)

data_spec_sort<-data_spec
data_spec_sort$LncRNA<-rownames(data_spec_sort)
data_spec_sort<-data_spec_sort[order(-data_spec_sort$TissSpecScor),]

# R checks whether the data conforms to the Normal distribution
#1  histogram
hist(data_spec_sort$TissSpecScor)
#2 Q-Q plot
qqnorm(data_spec_sort$TissSpecScor, main ="QQ")
qqline(data_spec_sort$TissSpecScor)
#3 shapiro.test()
nortest1<-shapiro.test(data_spec_sort$TissSpecScor)
nortest1
#Shapiro-Wilk normality test
#data:  data_spec_sort$TissSpecScor
#W = 0.99547, p-value = 0.0296
#p-value反应服从正态分布的概率，值越小越小的概率符合，通常0.05做标准，大于0.05则表示符合正态分布（此处为0.2542），故符合正态分布

library("dplyr")
hklncRNA<-filter(data_spec_sort, data_spec_sort$TissSpecScor < 0.5)
write.csv(hklncRNA, "HouseKeepingLncRNA.csv", quote = F)
  
  