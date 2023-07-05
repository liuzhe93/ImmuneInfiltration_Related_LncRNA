setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
rm(list=ls())

library(survival)
library("survminer")
rt=read.table("risk.txt",header=T,sep="\t",row.names = 1)

diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)


fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

summary(fit)    
###################################AUC###################################################
library(survivalROC)
pdf(file="ROC.pdf",width=6,height=6)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)


roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

###################################risk plot###################################################
library(pheatmap)
rt=read.table("risk.txt",header=T,sep="\t",row.names = 1)
rt=rt[order(rt$riskScore),]                                     

#risk curve
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
#line[line>10]=10
pdf(file="riskScore.pdf",width = 6,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
dev.off()

#survival status
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat.pdf",width = 6,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#risk heatmap
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 6,height = 3)
pheatmap(rt1, scale = "row",
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()

