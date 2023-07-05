setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())

lncRNA<-read.table("LncRNAannotation/GENCODE/gencode.v21.long_noncoding_RNAs.gtf", header = F, sep = "\t")
colnames(lncRNA)<-c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
library("dplyr")
library("tidyr")
lncRNA_sel<-lncRNA %>%
  separate(col = attributes,into = c("gene_id","transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type",
                                     "transcript_status", "transcript_name", "level", "havana_gene"), sep = ";")
lncRNA_sel$gene_name_cor<-lncRNA_sel$gene_name
for(i in 1:dim(lncRNA_sel)[1]){
  lncRNA_sel$gene_name_cor[i]<-gsub(" gene_name ","",lncRNA_sel$gene_name[i])
}
head(lncRNA_sel$gene_name_cor)
lncRNA_sel_gen<-subset(lncRNA_sel, lncRNA_sel$feature == "gene")

hklncRNA<-read.csv("HouseKeepingLncRNA.csv", header = T)
speclncRNA<-read.csv("SpecificityLncRNA.csv", header = T)

#xxx<-xx[order(-xx$Freq),]

hklncRNA_inf<-merge(hklncRNA, lncRNA_sel_gen, by.x = "X", by.y = "gene_name_cor")
write.csv(hklncRNA_inf, "hklncRNA_inf.csv", quote = F, row.names = F)

speclncRNA_inf<-merge(speclncRNA, lncRNA_sel_gen, by.x = "X", by.y = "gene_name_cor")
data_norm<-read.csv("data_norm.csv", header = T, row.names = 1)
data_norm_sel<-data_norm[speclncRNA_inf$X,]
data_norm_sel$Max<-rep("NA",dim(data_norm_sel)[1])
for(i in 1:dim(data_norm_sel)[1]){
  temp<-which.max(data_norm_sel[i,1:19])
  value<-colnames(t(temp))
  data_norm_sel[i,20] <- value
}
data_norm_sel$X<-rownames(data_norm_sel)
speclncRNA_fin<-merge(speclncRNA_inf, data_norm_sel, by = "X")
write.csv(speclncRNA_fin, "speclncRNA_inf.csv", quote = F, row.names = F)

