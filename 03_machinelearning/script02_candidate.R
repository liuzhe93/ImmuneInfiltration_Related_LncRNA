setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning")
rm(list=ls())

lasso<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/LassoLR/genelist_Lasso.csv")
rf<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/RandomForest/genelist_RandomForest.csv")
boruta<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/Boruta/genelist_Boruta.csv")
xgboost<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/Xgboost/genelist_Xgboost.csv")

LASSO<- lasso$x
RandomForest<- rf$x
Boruta<-boruta$x
Xgboost<-xgboost$x

length(LASSO)
#[1] 14
length(RandomForest)
#[1] 43
length(Boruta)
#[1] 6
length(Xgboost)
#[1] 46

list<-list(LASSO=LASSO,RandomForest=RandomForest,Boruta=Boruta,Xgboost=Xgboost)
library("ggplot2")
library("Vennerable")
library("ggVennDiagram")

venn = Venn(list)
data = process_data(venn)
pdf("Venn.pdf")
ggplot() +
  geom_sf(aes(fill = id), data = venn_region(data)) +
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()+
  theme(legend.position = 'none')+
  ggtitle("Veen plot for LncRNAs")
dev.off()

library("UpSetR")
pdf("UpsetPlot.pdf", width = 7, height = 4)
upset(fromList(list),nsets = 6,order.by = "freq")
dev.off()

merged_features<-c(LASSO,RandomForest,Boruta,Xgboost)
mydata<-table(merged_features)
mydata<-as.data.frame(mydata)
mydata_sel<-mydata[mydata$Freq>2,]
dim(mydata)
#[1] 60  2
dim(mydata_sel)
#[1] 12  2

write.csv(mydata_sel$merged_features, "value_lncRNAs.csv", quote = F, row.names = F)
