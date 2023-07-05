#reference：https://blog.csdn.net/tommyhechina/article/details/80356110

###############################microarry data organization############################################################
setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA")
rm(list=ls())
### raw data download
####obtain raw data by using getGEOSUppFiles() function
library("GEOquery")
library("affy")
#GSEName="GSE3910"
#rawdata = getGEOSuppFiles(GSEName)
#dir <- choose.dir(caption = "Select folder")
cel.files <- list.files(path = "/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA/immunecells_selected", 
                        pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
data.raw <- ReadAffy(filenames = cel.files)
sampleNames(data.raw)
est <- rma(data.raw)
expmtx_expso <- exprs(est)
write.csv(expmtx_expso, file="expmtx_expso.csv")

expmtx_expso<-read.csv("expmtx_expso.csv", header = T)
label_anno<-read.csv("LncRNA_probe.csv", header = T)
length(unique(label_anno$gene))
#[1] 1422
# the number of lncRNA genes: 1422
merged<-merge(label_anno, expmtx_expso, by.x = "probe", by.y = "X")
write.csv(merged, file="exprSet_probe_gene.csv", quote = F, row.names = F)

merged<-read.csv("exprSet_probe_gene.csv", header = T)
exprSet_symbol1 <- aggregate(x = merged[,3:ncol(merged)],
                             by = list(merged$gene),
                             FUN = median)
rownames(exprSet_symbol1)<-exprSet_symbol1$Group.1
dim(exprSet_symbol1)
exprSet_symbol1[1:5,1:5]
exprSet_symbol1$Group.1<-NULL
write.csv(exprSet_symbol1, file="exprSet_uniq.csv", quote = F)
sample_names<-read.csv("SamplesName.csv",header=T)
colnames(sample_names)<-c("sample","ImmuneCellType")
exprSet_symbol1_t<-t(exprSet_symbol1)
exprSet_symbol1_t<-as.data.frame(exprSet_symbol1_t)
exprSet_symbol1_t$sample<-row.names(exprSet_symbol1_t)
sample_names$sample_cor<-sample_names$sample
for(i in 1:dim(sample_names)[1]){
  sample_names$sample_cor[i]<-gsub("-","\\.",sample_names$sample[i])
}
exp_pro<-merge(sample_names,exprSet_symbol1_t,by.x="sample_cor", by.y="sample" )
exp_pro$sample_cor<-NULL
write.csv(exp_pro, "exp_trans_with_celltype.csv", quote = F, row.names = F)

exp_pro_t<-t(exp_pro)
write.csv(exp_pro_t, "exp_with_celltype.csv", quote = F)

