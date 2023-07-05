setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
rm(list=ls())

library(timeROC)
library(survival)
library(survivalROC)


mydata<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset/TCGA_exp_surv.csv")
mydata$futime<-mydata$futime/365
rt=mydata
rownames(rt)<-rt$id
rt<-rt[,-1]
#rt<-na.omit(rt)
#COX模型构建
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
#multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型参数
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab<-as.data.frame(outTab)
outTab$id<-rownames(outTab)
outTab<-subset(outTab, select = c("id","coef","HR","HR.95L","HR.95H","pvalue"))
write.table(outTab,file="multiCox.xls",sep="\t",row.names = F,quote=F)

#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)

risk<-read.table("risk.txt",header = T, row.names = 1)
time_roc_res <- timeROC(
  T = risk$futime,
  delta = risk$fustat,
  marker = risk$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(7, 8, 9),
  #times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)

time_roc_res$AUC
#      t=7       t=8       t=9 
#0.7027902 0.7795187 0.7377628 

confint(time_roc_res, level = 0.95)$CI_AUC
#     2.5% 97.5%
#t=7 56.56 84.00
#t=8 62.36 93.54
#t=9 56.21 91.35

time_ROC_df <- data.frame(
  TP_7year = time_roc_res$TP[, 1],
  FP_7year = time_roc_res$FP[, 1],
  TP_8year = time_roc_res$TP[, 2],
  FP_8year = time_roc_res$FP[, 2],
  TP_9year = time_roc_res$TP[, 3],
  FP_9year = time_roc_res$FP[, 3]
)
library(ggplot2)
pdf("ROC_curve_year789_TCGA.pdf")
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_7year, y = TP_7year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_8year, y = TP_8year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_9year, y = TP_9year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 7 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 8 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 9 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()


