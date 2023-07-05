setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/LassoLR")
rm(list=ls())

#####################################1. LassoLR######################################################
library("glmnet")
library("survival")
set.seed(1)
rt=read.table("../uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)       
rt$futime[rt$futime<=0]=1
rt<-na.omit(rt)
rt_cp<-apply(rt, 2, as.numeric)
rownames(rt_cp)<-rownames(rt)
rt<-rt_cp
rt<-as.data.frame(rt)
rt<-na.omit(rt)
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(as.numeric(rt$futime),as.numeric(rt$fustat)))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names=F,quote=F)
genelist_lasso<-colnames(lassoSigExp)[c(-1,-2,-3)]
length(genelist_lasso)
#[1] 13
write.csv(genelist_lasso, "genelist_Lasso.csv", quote = F, row.names = F)

