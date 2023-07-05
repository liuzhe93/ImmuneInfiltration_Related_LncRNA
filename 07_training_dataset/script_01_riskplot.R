setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/07_survival_analysis")
rm(list=ls()) 
library(survival)
library("survminer")
rt=read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",header=T,sep="\t",row.names = 1)
rt$futime<-rt$futime/12

diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)


fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
###################################绘制生存曲线###################################################
pdf(file="survival.pdf",onefile = FALSE,
    width = 7,             #图片的宽度
    height =6)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()

summary(fit)    #查看五年生存率
###################################模型评估：AUC###################################################
library(survivalROC)
pdf(file="ROC.pdf",width=6,height=6)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)


uniSigExp<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/04_univariateCox_regression/uniSigExp.txt",
                      sep="\t",header=T)
genes<-c("CYB561D2", "LINC00638", "DANCR")
uniSigExp_filtered<-uniSigExp[,c("id","futime","fustat",genes)]
dim(uniSigExp_filtered)
#[1] 566   6
mydata <- na.omit(uniSigExp_filtered) # 删除缺失值
dim(mydata)
#[1] 562   14
row.names(mydata)<-mydata$id
mydata<-mydata[,-1]
dim(mydata)
#[1] 562  5

rt=mydata
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
rt<-cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk))
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

###################################模型评估：风险图###################################################
library(pheatmap)
rt=read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",header=T,sep="\t",row.names = 1)
rt=rt[order(rt$riskScore),]                                     #按照riskScore对样品排序

#绘制风险曲线
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

#绘制生存状态图
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

#绘制风险热图
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

