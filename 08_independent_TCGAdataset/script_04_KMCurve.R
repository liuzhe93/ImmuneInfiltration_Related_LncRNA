setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
rm(list=ls())

risk<-read.table("risk.txt", sep="\t",header=T, row.names = 1)
head(risk)

library("survival")
library("survminer")
rt=risk
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

#plot survival curve
pdf(file="survival.pdf",onefile = FALSE,
    width = 5.5,             
    height =5)             
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


