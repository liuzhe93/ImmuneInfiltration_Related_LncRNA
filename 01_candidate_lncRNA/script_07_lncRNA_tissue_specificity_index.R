#reference：https://blog.csdn.net/tommyhechina/article/details/80356110

###############################microarry data organization############################################################
setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())

data_norm<-read.csv("data_norm.csv", header = T, row.names = 1)
n = ncol(data_norm)
data_spec<-matrix(data=rep(0,n), nrow = n, ncol = 1)
data_spec<-as.data.frame(data_spec)
colnames(data_spec)<-"TissSpecScor"
for(i in 1:nrow(data_norm)){
  for(j in 1:ncol(data_norm)){
    temp = (1-data_norm[i,j])/(n-1)
    data_spec[i,1] = temp + data_spec[i,1]
  }
}
rownames(data_spec)<-row.names(data_norm)
write.csv(data_spec, "data_tissueSpecScor.csv", quote = F)

data_spec<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA/data_tissueSpecScor.csv")
data_spec_sort<-data_spec[order(-data_spec$TissSpecScor),]
num_sel1<-round(dim(data_spec_sort)[1]*0.3)
num_sel1
#[1] 221
speclncRNA<-data_spec_sort[1:num_sel1,]
num_sel2<-(dim(data_spec_sort)[1]-num_sel1+1)
num_sel2
hklncRNA<-data_spec_sort[num_sel2:dim(data_spec_sort)[1],]

write.csv(hklncRNA, "HouseKeepingLncRNA.csv", quote = F, row.names = F)
write.csv(speclncRNA, "SpecificityLncRNA.csv", quote = F, row.names = F)



