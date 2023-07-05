setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/05_Gaussian_mixture_mode")
rm(list=ls())

library("gtools")
cmat_1 <- combinations(11,1)
cmat_2 <- combinations(11,2)
cmat_3 <- combinations(11,3)
cmat_4 <- combinations(11,4)
cmat_5 <- combinations(11,5)
cmat_6 <- combinations(11,6)
cmat_7 <- combinations(11,7)
cmat_8 <- combinations(11,8)
cmat_9 <- combinations(11,9)
cmat_10 <- combinations(11,10)
cmat_11 <- combinations(11,11)


uniSigExp<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/04_univariateCox_regression/uniSigExp.txt",
                      sep="\t",header=T)
genes<-colnames(uniSigExp)[4:ncol(uniSigExp)]
uniSigExp_filtered<-uniSigExp[,c("id","futime","fustat",genes)]
dim(uniSigExp_filtered)
#[1] 566   14
mydata <- na.omit(uniSigExp_filtered) # 删除缺失值
dim(mydata)
#[1] 562   14
row.names(mydata)<-mydata$id
mydata<-mydata[,-1]
dim(mydata)
#[1] 562  13


library(survival)
library(survivalROC)

temp<-mydata[,c("futime","fustat")]
mydata_case1=data.frame()
result=data.frame()
for(i in 1:11){
  j = i+2
  mydata_case1<-mydata[,j]
  gene_name<-colnames(mydata)[j]
  dat<-cbind(temp,mydata_case1)
  colnames(dat)[3]<-gene_name
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=gene_name, auc=roc$AUC))
}
write.csv(result,"case1_AUC.csv",quote = F, row.names = F)

mydata_case2=data.frame()
result=data.frame()
for(i in 1:10){
  j = i+2
  mydata_case2_1<-mydata[,j]
  gene_name_1<-colnames(mydata)[j]
  dat<-cbind(temp,mydata_case2_1)
  colnames(dat)[3]<-gene_name_1
  a = 1:11
  if(i==1){
    b = a[-1]
  }else if(i==2){
    b = a[-1]
    b = b[-1]
  }else if(i==3){
    b = a[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==4){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==5){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==6){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==7){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==8){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==9){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }else if(i==10){
    b = a[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
    b = b[-1]
  }
  for(m in b){
    n = m+2
    dat<-dat[,1:3]
    mydata_case2_2<-mydata[,n]
    gene_name_2<-colnames(mydata)[n]
    dat<-cbind(dat,mydata_case2_2)
    colnames(dat)[4]<-gene_name_2
    rt=dat
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
    result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,sep=" + "), auc=roc$AUC))
  }
}
write.csv(result,"case2_AUC.csv",quote = F, row.names = F)

cmat_3 <- combinations(11,3)
temp<-mydata[,c("futime","fustat")]
mydata_case3=data.frame()
result=data.frame()
for(i in 1:165){
  loc_1 <- cmat_3[i,1]
  loc_2 <- cmat_3[i,2]
  loc_3 <- cmat_3[i,3]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  mydata_case3<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3])
  colnames(mydata_case3)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3)
  dat<-mydata_case3
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case3_AUC.csv",quote = F, row.names = F)

cmat_4 <- combinations(11,4)
temp<-mydata[,c("futime","fustat")]
mydata_case4=data.frame()
result=data.frame()
for(i in 1:330){
  loc_1 <- cmat_4[i,1]
  loc_2 <- cmat_4[i,2]
  loc_3 <- cmat_4[i,3]
  loc_4 <- cmat_4[i,4]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  mydata_case4<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4])
  colnames(mydata_case4)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4)
  dat<-mydata_case4
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case4_AUC.csv",quote = F, row.names = F)

cmat_5 <- combinations(11,5)
temp<-mydata[,c("futime","fustat")]
mydata_case5=data.frame()
result=data.frame()
for(i in 1:462){
  loc_1 <- cmat_5[i,1]
  loc_2 <- cmat_5[i,2]
  loc_3 <- cmat_5[i,3]
  loc_4 <- cmat_5[i,4]
  loc_5 <- cmat_5[i,5]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  gene_name_5 <- genes[loc_5]
  mydata_case5<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4],mydata[,gene_name_5])
  colnames(mydata_case5)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5)
  dat<-mydata_case5
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case5_AUC.csv",quote = F, row.names = F)


cmat_6 <- combinations(11,6)
temp<-mydata[,c("futime","fustat")]
mydata_case6=data.frame()
result=data.frame()
for(i in 1:462){
  loc_1 <- cmat_6[i,1]
  loc_2 <- cmat_6[i,2]
  loc_3 <- cmat_6[i,3]
  loc_4 <- cmat_6[i,4]
  loc_5 <- cmat_6[i,5]
  loc_6 <- cmat_6[i,6]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  gene_name_5 <- genes[loc_5]
  gene_name_6 <- genes[loc_6]
  mydata_case6<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4],
                      mydata[,gene_name_5],mydata[,gene_name_6])
  colnames(mydata_case6)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,gene_name_6)
  dat<-mydata_case6
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,
                                                 gene_name_6,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case6_AUC.csv",quote = F, row.names = F)


cmat_7 <- combinations(11,7)
temp<-mydata[,c("futime","fustat")]
mydata_case7=data.frame()
result=data.frame()
for(i in 1:330){
  loc_1 <- cmat_7[i,1]
  loc_2 <- cmat_7[i,2]
  loc_3 <- cmat_7[i,3]
  loc_4 <- cmat_7[i,4]
  loc_5 <- cmat_7[i,5]
  loc_6 <- cmat_7[i,6]
  loc_7 <- cmat_7[i,7]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  gene_name_5 <- genes[loc_5]
  gene_name_6 <- genes[loc_6]
  gene_name_7 <- genes[loc_7]
  mydata_case7<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4],
                      mydata[,gene_name_5],mydata[,gene_name_6],mydata[,gene_name_7])
  colnames(mydata_case7)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,gene_name_6,
                            gene_name_7)
  dat<-mydata_case7
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,
                                                 gene_name_6,gene_name_7,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case7_AUC.csv",quote = F, row.names = F)

cmat_8 <- combinations(11,8)
temp<-mydata[,c("futime","fustat")]
mydata_case8=data.frame()
result=data.frame()
for(i in 1:165){
  loc_1 <- cmat_8[i,1]
  loc_2 <- cmat_8[i,2]
  loc_3 <- cmat_8[i,3]
  loc_4 <- cmat_8[i,4]
  loc_5 <- cmat_8[i,5]
  loc_6 <- cmat_8[i,6]
  loc_7 <- cmat_8[i,7]
  loc_8 <- cmat_8[i,8]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  gene_name_5 <- genes[loc_5]
  gene_name_6 <- genes[loc_6]
  gene_name_7 <- genes[loc_7]
  gene_name_8 <- genes[loc_8]
  mydata_case8<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4],
                      mydata[,gene_name_5],mydata[,gene_name_6],mydata[,gene_name_7],mydata[,gene_name_8])
  colnames(mydata_case8)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,gene_name_6,
                            gene_name_7,gene_name_8)
  dat<-mydata_case8
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,
                                                 gene_name_6,gene_name_7,gene_name_8,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case8_AUC.csv",quote = F, row.names = F)


cmat_9 <- combinations(11,9)
temp<-mydata[,c("futime","fustat")]
mydata_case9=data.frame()
result=data.frame()
for(i in 1:55){
  loc_1 <- cmat_9[i,1]
  loc_2 <- cmat_9[i,2]
  loc_3 <- cmat_9[i,3]
  loc_4 <- cmat_9[i,4]
  loc_5 <- cmat_9[i,5]
  loc_6 <- cmat_9[i,6]
  loc_7 <- cmat_9[i,7]
  loc_8 <- cmat_9[i,8]
  loc_9 <- cmat_9[i,9]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  gene_name_5 <- genes[loc_5]
  gene_name_6 <- genes[loc_6]
  gene_name_7 <- genes[loc_7]
  gene_name_8 <- genes[loc_8]
  gene_name_9 <- genes[loc_9]
  mydata_case9<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4],
                      mydata[,gene_name_5],mydata[,gene_name_6],mydata[,gene_name_7],mydata[,gene_name_8],
                      mydata[,gene_name_9])
  colnames(mydata_case9)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,gene_name_6,
                            gene_name_7,gene_name_8,gene_name_9)
  dat<-mydata_case9
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,
                                                 gene_name_6,gene_name_7,gene_name_8,gene_name_9,sep=" + "), auc=roc$AUC))
}
write.csv(result,"case9_AUC.csv",quote = F, row.names = F)


cmat_10 <- combinations(11,10)
temp<-mydata[,c("futime","fustat")]
mydata_case10=data.frame()
result=data.frame()
for(i in 1:11){
  loc_1 <- cmat_10[i,1]
  loc_2 <- cmat_10[i,2]
  loc_3 <- cmat_10[i,3]
  loc_4 <- cmat_10[i,4]
  loc_5 <- cmat_10[i,5]
  loc_6 <- cmat_10[i,6]
  loc_7 <- cmat_10[i,7]
  loc_8 <- cmat_10[i,8]
  loc_9 <- cmat_10[i,9]
  loc_10 <- cmat_10[i,10]
  gene_name_1 <- genes[loc_1]
  gene_name_2 <- genes[loc_2]
  gene_name_3 <- genes[loc_3]
  gene_name_4 <- genes[loc_4]
  gene_name_5 <- genes[loc_5]
  gene_name_6 <- genes[loc_6]
  gene_name_7 <- genes[loc_7]
  gene_name_8 <- genes[loc_8]
  gene_name_9 <- genes[loc_9]
  gene_name_10 <- genes[loc_10]
  mydata_case10<-cbind(temp,mydata[,gene_name_1],mydata[,gene_name_2],mydata[,gene_name_3],mydata[,gene_name_4],
                      mydata[,gene_name_5],mydata[,gene_name_6],mydata[,gene_name_7],mydata[,gene_name_8],
                      mydata[,gene_name_9],mydata[,gene_name_10])
  colnames(mydata_case10)<-c("futime","fustat",gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,gene_name_6,
                            gene_name_7,gene_name_8,gene_name_9,gene_name_10)
  dat<-mydata_case10
  rt=dat
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
  result<-rbind(result,cbind(model=i, gene=paste(gene_name_1,gene_name_2,gene_name_3,gene_name_4,gene_name_5,
                                                 gene_name_6,gene_name_7,gene_name_8,gene_name_9,gene_name_10,
                                                 sep=" + "), auc=roc$AUC))
}
write.csv(result,"case10_AUC.csv",quote = F, row.names = F)


dat<-mydata
rt=dat
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
gene_names<-paste("CYB561D2","LINC00638","LINC01119", "ADARB2.AS1", "GABARAPL3", "PRR34.AS1", "OVCH1.AS1","DANCR","DSCR10",
                  "DSCR9","LINC01208",sep=" + ")
result<-cbind(model=i, gene=gene_names, auc=roc$AUC)
write.csv(result,"case11_AUC.csv",quote = F, row.names = F)


library(survival)                                        #引用包
rt=mydata  #读取输入文件
#COX模型构建
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型参数
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

#输出病人风险值
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

########################################最后整合以上的结果############################################################
result_1<-read.csv("case1_AUC.csv")
result_2<-read.csv("case2_AUC.csv")
result_3<-read.csv("case3_AUC.csv")
result_4<-read.csv("case4_AUC.csv")
result_5<-read.csv("case5_AUC.csv")
result_6<-read.csv("case6_AUC.csv")
result_7<-read.csv("case7_AUC.csv")
result_8<-read.csv("case8_AUC.csv")
result_9<-read.csv("case9_AUC.csv")
result_10<-read.csv("case10_AUC.csv")
result_11<-read.csv("case11_AUC.csv")

cluster_data<-rbind(result_1,result_2,result_3,result_4,result_5,result_6,result_7,result_8,result_9,result_10,result_11)
for(i in 1:nrow(cluster_data)){
  cluster_data$model[i]=i
}

mydata<-cluster_data
mydata$auc <- scale(mydata$auc) # 数据标准化
# Model Based Clustering
library(mclust)
fit <- Mclust(mydata)
plot(fit) # plot results
##输入"1"：代表BIC可视化
##输入"0"：代表退出
summary(fit) # display the best model
#Gaussian finite mixture model fitted by EM algorithm 
#---------------------------------------------------- 
#  
#  Mclust VVV (ellipsoidal, varying volume, shape, and orientation) model with 9 components: 
#  
#  log-likelihood    n df       BIC       ICL
#       -33054.86 2047 89 -66788.26 -67134.47
#
#Clustering table:
#  1   2   3   4   5   6   7   8   9 
#327 334 120 247 296 109 261 132 221

cluster_data$cluster<-1:2047
cluster_data$numGenes<-1:2047
for(i in 1:2047){
  cluster_data$numGenes[i]<-length(strsplit(cluster_data$gene[i],split="\\ \\+\\ ")[[1]])
  cluster_data$cluster[i]<-paste0("cluster_",fit$classification[i])
}
write.csv(cluster_data,"cluster_results.csv",quote=F,row.names = F)
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
pdf("CombinationGenesAndAUC.pdf")
ggplot(data=cluster_data, aes(x=model, y=auc, color=cluster))+geom_point(size=3)+ 
  stat_ellipse(aes(fill=cluster),level = 0.95,alpha=0.2,geom = "polygon")
dev.off()



