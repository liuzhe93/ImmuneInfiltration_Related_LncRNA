setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/RandomForest")
rm(list=ls())


#################################5. RandomForest#####################################################
library("survival")
library("randomForestSRC")
set.seed(1)
uniSigExp<-read.table("../uniSigExp.txt",sep="\t",header=T,row.names = 1)
uniSigExp<-cbind(uniSigExp$fustat,uniSigExp$futime,uniSigExp[,3:dim(uniSigExp)[2]])
colnames(uniSigExp)<-c("fustat","futime",colnames(uniSigExp)[3:dim(uniSigExp)[2]])
pbc<-uniSigExp
pbc<-apply(uniSigExp, 2, as.numeric)
rownames(pbc)<-rownames(uniSigExp)
pbc<-as.data.frame(pbc)
pbc.obj <- rfsrc(Surv(futime, fustat) ~ ., pbc, nodesize = 20, importance = TRUE)
vs.pbc <- var.select(object = pbc.obj, conservative = "high")
topvars <- vs.pbc$topvars
if (library("survival", logical.return = TRUE))
{
  cox.weights <- function(rfsrc.f, rfsrc.data) {
    event.names <- all.vars(rfsrc.f)[1:2]
    p <- ncol(rfsrc.data) - 2
    event.pt <- match(event.names, names(rfsrc.data))
    xvar.pt <- setdiff(1:ncol(rfsrc.data), event.pt)
    sapply(1:p, function(j) {
      cox.out <- coxph(rfsrc.f, rfsrc.data[, c(event.pt, xvar.pt[j])])
      pvalue <- summary(cox.out)$coef[5]
      if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
    })
  }
  data_rf=pbc  #输入数据
  rfsrc.f <- as.formula(Surv(futime, fustat) ~ .)  
  cox.wts <- cox.weights(rfsrc.f, data_rf)
  b.obj <- rfsrc(rfsrc.f, data_rf , nsplit = 10, xvar.wt = cox.wts,importance = "random",na.action ="na.impute",ntree = 1000)
  vh.breast.cox <- var.select(rfsrc.f, data_rf, method = "vh", nstep = 5,nrep = 10, xvar.wt = cox.wts,conservative = "high")
}
png(file="rf.png",width=3000,height=3000)
random.plot=plot(b.obj)
vimp=as.data.frame(b.obj$importance)
vimp$symbol=rownames(vimp)
vimp=vimp[order(vimp$`b.obj$importance`,decreasing = T),]
head(vimp)
dev.off()
write.csv(vimp,"gene_important.csv",quote=F)
pdf("gene_important.pdf", width=24,height = 12)
random.plot=plot(b.obj)
dev.off()
genelist_RandomForest<-vimp[which(vimp$`b.obj$importance`>0),]
genelist_RandomForest<-genelist_RandomForest$symbol
length(genelist_RandomForest)
#[1] 43
genelist_RandomForest<-gsub("\\.","-",genelist_RandomForest)
write.csv(genelist_RandomForest, "genelist_RandomForest.csv", quote = F, row.names = F)

