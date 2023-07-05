setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/09_ImmuneInfiltration")
rm(list=ls())
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(msigdbr) 
library(org.Hs.eg.db)
library(GSVA)

markergenes<-read.csv("markergene.csv",header=T)
markergenes$term<-markergenes$Cell.type
markergenes$gene<-markergenes$Metagene
markergenes<-subset(markergenes,select=c("term","gene"))

GeneID2kegg_list<- tapply(markergenes$term,as.factor(markergenes$gene),function(x) x)
kegg2GeneID_list<- tapply(markergenes$gene,as.factor(markergenes$term),function(x) x)

geneSet=kegg2GeneID_list
gmt_file='kegg2symbol.gmt'
write.gmt <- function(geneSet,gmt_file){
  sink( gmt_file )
  lapply(names(geneSet), function(i){
    cat( paste(c(i,"NA",geneSet[[i]]),collapse = "\t"))
    cat("\n")
  })
  sink()
}
write.gmt(geneSet,gmt_file)



cancer_exp<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/cancer_exp.csv")
immune_exp<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/immune_exp.csv")
merged_data<-merge(immune_exp, cancer_exp, by = "X")
rownames(merged_data)<-merged_data$X


normal_labels<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/immune_label.csv")
normal_labels<-normal_labels[,-1]
dim(normal_labels)
#[1] 115   2
cancer_labels<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/coloncancer_label.csv")
cancer_labels<-cancer_labels[,-1]
dim(cancer_labels)
#[1] 566   2
normal_labels$samplename<-gsub("-",".",normal_labels$colnames.expr.rma.)
normal<-merged_data[,c("X",normal_labels$samplename)]
rownames(normal)<-normal$X
normal<-normal[,-1]
cancer_labels$samplename<-gsub("-",".",cancer_labels$colnames.expr.rma.)
cancer<-merged_data[,c("X",cancer_labels$samplename)]
rownames(cancer)<-cancer$X
cancer<-cancer[,-1]

normal<-as.matrix(normal)
cancer<-as.matrix(cancer)
exp<-cbind(normal,cancer)

library(GSEABase)
geneSet<-getGmt("kegg2symbol.gmt")
gsva_matrix = gsva(exp, geneSet, method = "ssgsea", kcdf = "Gaussian", abs.ranking=F) 
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}
nor_gsva_matrix1<-normalization(gsva_matrix)
write.table(nor_gsva_matrix1,file="norm_hallmark_ssGSEA_Score.txt",sep="\t",quote=F)

risk_score<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",sep="\t",header=T)
rownames(risk_score)<-risk_score$id
risk_score<-subset(risk_score,select=c("riskScore"))
AS<-risk_score
SF<-nor_gsva_matrix1
for(i in 1:ncol(SF)){
  colnames(SF)[i] <- strsplit(colnames(SF),"_")[[i]][1]
}

sameSample<-intersect(colnames(SF),rownames(AS))
AS_t<-t(AS)
SF1<-SF
SF1<-SF1[,sameSample]
AS1<-AS_t

outTab=data.frame()
for(i in row.names(SF1)){
  for(j in row.names(AS1)){
    x=as.numeric(SF1[i,])
    y=as.numeric(AS1[j,])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if(cor>0){
      outTab=rbind(outTab,cbind(SF=i,AS=j,cor,pvalue,Regulation="1_positive"))
    }
    if(cor<0){
      outTab=rbind(outTab,cbind(SF=i,AS=j,cor,pvalue,Regulation="2_negative"))
    }
  }
}
outTab$y_value=-log10(as.numeric(outTab$pvalue))
outTab$x_value=as.numeric(outTab$cor)
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names = F)
write.csv(file="corResults.csv",outTab,quote=F,row.names = F)

library(ggrepel)
pdf("Volcano.pdf",width=6,height = 8)
ggplot(data=outTab)+
  geom_point(aes(x_value,y_value,colour = Regulation))+
  geom_text_repel(aes(x_value,y_value,label = SF))+
  theme_classic(base_size = 16)

dev.off()

