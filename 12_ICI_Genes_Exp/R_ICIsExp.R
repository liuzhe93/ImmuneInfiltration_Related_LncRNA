setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/12_ICI_Genes_Exp")
rm(list=ls())

library("survival")
library("survminer")


rt=read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",header=T,sep="\t",row.names = 1)
rt$sample<-rownames(rt)
#rt$sample<-rt$id

microarrary<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/cancer_exp.csv")
rownames(microarrary)<-microarrary$X
microarrary<-microarrary[,-1]
microarrary_t<-t(microarrary)
microarrary_t<-as.data.frame(microarrary_t)
microarrary_t$SampleName<-row.names(microarrary_t)
exp_PD_1<-microarrary_t[,"PDCD1"]
exp_PD_L1<-microarrary_t[,"CD274"]
exp_CTLA4<-microarrary_t[,"CTLA4"]
exp<-cbind(microarrary_t$SampleName,exp_PD_1,exp_PD_L1,exp_CTLA4)
colnames(exp)<-c("SampleName","PD_1","PD_L1","CTLA_4")
exp<-as.data.frame(exp)
exp$sample<-substr(exp$SampleName,1,9)
array_data<-merge(exp,rt,by="sample")
array_data<-subset(array_data,select = c("PD_1","PD_L1","CTLA_4","riskScore","risk"))
array_data$PD_1<-as.numeric(array_data$PD_1)
array_data$PD_L1<-as.numeric(array_data$PD_L1)
array_data$CTLA_4<-as.numeric(array_data$CTLA_4)


library(ggpubr)
pdf("PD-1_immune.pdf")
ggboxplot(array_data, x="risk", y="PD_1",
          color="risk", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
  #  geom_hline(yintercept = mean(merged_risk_filtered$Score),
  #             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "wilcox.test", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                     # Pairwise comparison against all
dev.off()

pdf("PD-L1_immune.pdf")
ggboxplot(array_data, x="risk", y="PD_L1",
          color="risk", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
  #  geom_hline(yintercept = mean(merged_risk_filtered$Score),
  #             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "wilcox.test", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                     # Pairwise comparison against all
dev.off()

pdf("CTLA-4_immune.pdf")
ggboxplot(array_data, x="risk", y="CTLA_4",
          color="risk", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
  #  geom_hline(yintercept = mean(merged_risk_filtered$Score),
  #             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "wilcox.test", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                     # Pairwise comparison against all
dev.off()

library("ggExtra")
library("ggplot2")
library("ggpubr")
library("ggpmisc")
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
xval<-array_data$riskScore
yval<-array_data$PD_1
data_PD_1<-cbind(xval,yval)
colnames(data_PD_1)<-c("riskScore","PD_1")
data_PD_1<-as.data.frame(data_PD_1)
pdf("PD_1_cor.pdf")
p<-ggscatter(data_PD_1,x="riskScore",y="PD_1",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()


xval<-array_data$riskScore
yval<-array_data$PD_L1
data_PD_L1<-cbind(xval,yval)
colnames(data_PD_L1)<-c("riskScore","PD_L1")
data_PD_L1<-as.data.frame(data_PD_L1)
pdf("PD_L1_cor.pdf")
p<-ggscatter(data_PD_L1,x="riskScore",y="PD_L1",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()


xval<-array_data$riskScore
yval<-array_data$CTLA_4
data_CTLA_4<-cbind(xval,yval)
colnames(data_CTLA_4)<-c("riskScore","CTLA_4")
data_CTLA_4<-as.data.frame(data_CTLA_4)
pdf("CTLA_4_cor.pdf")
p<-ggscatter(data_CTLA_4,x="riskScore",y="CTLA_4",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()

