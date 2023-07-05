setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression")
rm(list=ls())

com_auc<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/05_Gaussian_mixture_mode/cluster_results.csv",
                  header=T)
com_auc_sorted=com_auc[order(com_auc$auc,decreasing = T),]
genes<-com_auc_sorted$gene[1]
gene_list<-strsplit(genes,split="\\ \\+\\ ")[[1]]


uniSigExp<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/04_univariateCox_regression/uniSigExp.txt",
                      sep="\t",header=T)
genes<-colnames(uniSigExp)[4:ncol(uniSigExp)]
uniSigExp_filtered<-uniSigExp[,c("id","futime","fustat",genes)]
dim(uniSigExp_filtered)
#[1] 566   14
mydata <- na.omit(uniSigExp_filtered) 
dim(mydata)
#[1] 562   14
dim(mydata)
#[1] 562  14
uniSigExp_filter<-subset(mydata, select = c("id", "futime","fustat", gene_list))

library("survival")
rt=uniSigExp_filter
rownames(rt)<-rt$id
rt<-rt[,-1]

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
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)

