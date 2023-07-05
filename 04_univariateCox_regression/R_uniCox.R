setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/04_univariateCox_regression")
rm(list=ls())

pFilter=0.05                                                      #定义单因素显著性
library(survival)                                                 #引用包
feature<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/value_lncRNAs.csv")
lncRNAs<-feature$x
lncRNAs<-gsub("\\-","\\.",lncRNAs)

mydata<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/exp_cli.csv")
mydata_sel<-subset(mydata, select = c("Sample_geo_accession", "os.delay..months.", "os.event",lncRNAs))
colnames(mydata_sel)[1:3]<-c("GEO","futime","fustat")
rownames(mydata_sel)<-mydata_sel$GEO
mydata_sel$GEO<-NULL

outTab=data.frame()
mydata_sel_cp<-apply(mydata_sel, 2, as.numeric)
rownames(mydata_sel_cp)<-row.names(mydata_sel)
rt<-mydata_sel_cp
rt<-as.data.frame(rt)
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
dim(outTab)
#[1] 12   6
dim(sigTab)
#[1] 11  6
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)
uniSigExp<-read.table("uniSigExp.txt",sep="\t",header=T)
dim(uniSigExp)
#[1] 566  14

library("ggplot2")
head(outTab)
outTab$z<-as.numeric(outTab$z)
outTab$pvalue<-as.numeric(outTab$pvalue)
outTab$color<-ifelse(outTab$pvalue<0.05 ,"red","grey")
color<-c(red = "red", grey = "grey")
rownames(outTab)<-outTab$id
library(ggrepel)
p<-ggplot(outTab, aes(z, -log10(pvalue), col = color)) +
  geom_point() +
  geom_text_repel(aes(z, -log10(pvalue), label = rownames(outTab)))+
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="Univariate Cox coefficient", y = "-log10(P-Value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdf("UnivariateCOXRegressionAnalysis.pdf", width = 4, height = 5)
p
dev.off()





