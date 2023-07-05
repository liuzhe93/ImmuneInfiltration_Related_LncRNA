setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/09_IndependentAnalysis")
rm(list=ls())

clinical<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/09_IndependentAnalysis/SampleInfor.csv")
clinical<-subset(clinical,select=c("Sample_geo_accession","Gender","Age","tnm.stage","tnm.t","tnm.n","tnm.m","os.event","os.delay"))
colnames(clinical)<-c("Sample","Gender","Age","TNM_stage","TNM_t","TNM_n","TNM_m","fustat","futime")
clinical$Gender<-gsub("Sex: ","",clinical$Gender)
clinical$Age<-gsub("age.at.diagnosis \\(year\\): ","",clinical$Age)
clinical$TNM_stage<-gsub("tnm.stage: ","",clinical$TNM_stage)
clinical$TNM_t<-gsub("tnm.t: ","",clinical$TNM_t)
clinical$TNM_n<-gsub("tnm.n: ","",clinical$TNM_n)
clinical$TNM_m<-gsub("tnm.m: ","",clinical$TNM_m)
clinical$fustat<-gsub("os.event: ","",clinical$fustat)
clinical$futime<-gsub("os.delay \\(months\\): ","",clinical$futime)
clinical$futime<-as.numeric(clinical$futime)/12

risk<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",sep="\t",header=T)

merged_data<-merge(risk,clinical,by.x="id",by.y="Sample")
merged_data<-subset(merged_data,select=c("id","futime.y","fustat.x","Age","Gender","TNM_stage","TNM_t","TNM_n","TNM_m","riskScore"))
colnames(merged_data)<-c("id","futime","fustat","Age","Gender","TNM_stage","TNM_t","TNM_n","TNM_m","riskScore")
merged_data$Gender<-ifelse(merged_data$Gender=="Female",0,1)
merged_data$TNM_t<-ifelse(merged_data$TNM_t=='T0',1,
                          ifelse(merged_data$TNM_t=='T1',2,
                                 ifelse(merged_data$TNM_t=='T2',3,
                                        ifelse(merged_data$TNM_t=='T3',4,
                                               ifelse(merged_data$TNM_t=='T4',5,
                                                      "NA")))))
merged_data$TNM_n<-ifelse(merged_data$TNM_n=='N0',1,
                          ifelse(merged_data$TNM_n=='N1',2,
                                 ifelse(merged_data$TNM_n=='N2',3,
                                        ifelse(merged_data$TNM_n=='N3',4,
                                               "NA"))))
merged_data$TNM_m<-ifelse(merged_data$TNM_m=='M0',1,
                          ifelse(merged_data$TNM_m=='M1',2,
                                 ifelse(merged_data$TNM_m=='MX',3,
                                        "NA")))

merged_data<-na.omit(merged_data)
merged_data$Age<-as.numeric(merged_data$Age)
merged_data$Gender<-as.numeric(merged_data$Gender)
merged_data$TNM_stage<-as.numeric(merged_data$TNM_stage)
merged_data$TNM_t<-as.numeric(merged_data$TNM_t)
merged_data$TNM_n<-as.numeric(merged_data$TNM_n)
merged_data$TNM_m<-as.numeric(merged_data$TNM_m)
merged_data<-na.omit(merged_data)
write.table(merged_data,"indepInput.txt",sep="\t",row.names = F)

###########################################univariate test########################################################################
library(survival)
library(forestplot)

clrs <- fpColors(box="green",line="darkblue", summary="royalblue")             
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

rt=read.table("uniCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
pdf(file="forest_uni.pdf",onefile = FALSE,
    width = 6,             
    height = 4,            
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()


###########################################multivariate test########################################################################
library(survival)
library(forestplot)

clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
pdf(file="forest_multi.pdf",onefile = FALSE,
    width = 6,             
    height = 4,            
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()
