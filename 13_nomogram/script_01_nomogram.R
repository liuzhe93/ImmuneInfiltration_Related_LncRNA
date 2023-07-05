setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/13_nomogram")
rm(list=ls())

library(rms)
library(timeROC)
library(survival)
library(regplot)

risk<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",
                 sep="\t",header=T, row.names = 1)
risk$Sample<-row.names(risk)
head(risk)


clinical<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/10_IndependentAnalysis/SampleInfor.csv")
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


merged_data<-merge(risk, clinical, by = "Sample")
merged_data<-subset(merged_data, select = c("Sample","futime.y","fustat.x","CYB561D2","LINC00638","DANCR","riskScore",
                                            "risk","Gender","Age","TNM_stage", "TNM_t", "TNM_n", "TNM_m"))
colnames(merged_data)<-c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore","risk","Gender",
                         "Age","Stage", "TNM_t", "TNM_n", "TNM_m")
merged_data$Gender_num<-ifelse(merged_data$Gender=="Female","2",
                               ifelse(merged_data$Gender=="Male","1",
                                      "NA"))
merged_data$Age_num<-ifelse(merged_data$Age>=60,"1",
                            ifelse(merged_data$Age<60,"2",
                                   "NA"))
merged_data$Stage_num<-merged_data$Stage
merged_data$TNM_t_num<-ifelse(merged_data$TNM_t=="T0", "0",
                              ifelse(merged_data$TNM_t=="T1", "1",
                                     ifelse(merged_data$TNM_t=="T2", "2",
                                            ifelse(merged_data$TNM_t=="T3", "3",
                                                   ifelse(merged_data$TNM_t=="T4", "4",
                                                          "NA")))))
merged_data$TNM_n_num<-ifelse(merged_data$TNM_n=="N0", "0",
                              ifelse(merged_data$TNM_n=="N1", "1",
                                     ifelse(merged_data$TNM_n=="N2", "2",
                                            ifelse(merged_data$TNM_n=="N3", "3",
                                                   "NA"))))
merged_data$TNM_m_num<-ifelse(merged_data$TNM_m=="M0", "0",
                              ifelse(merged_data$TNM_m=="M1", "1",
                                     ifelse(merged_data$TNM_m=="MX", "2",
                                            "NA")))
head(merged_data)
data_selected<-subset(merged_data, select = c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore",
                                              "risk","Gender","Age","Stage","TNM_t", "TNM_n", "TNM_m" ))
colnames(data_selected)<-c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore","risk","Gender",
                           "Age","Stage","TNM_t","TNM_n","TNM_m")
rt<-data_selected
head(rt)
rownames(rt)<-rt$sample_name
rt<-subset(rt, select = c("futime", "Status", "Age", "Gender", "Stage", "TNM_t","TNM_n","TNM_m", "riskScore"))
colnames(rt)<-c("futime", "fustat", "Age", "Gender", "Stage", "TNM_t","TNM_n","TNM_m", "riskScore")
rt$Stage<-ifelse(rt$Stage=="0","Stage 0",
                 ifelse(rt$Stage=="1","Stage I",
                        ifelse(rt$Stage=="2","Stage II",
                               ifelse(rt$Stage=="3","Stage III",
                                      ifelse(rt$Stage=="4","Stage IV",
                                             "NA")))))
head(rt)
write.csv(rt, "nomogram_risk.csv", quote = F)

rt$Age<-as.numeric(rt$Age)
ddist <- datadist(rt)
options(datadist='ddist')

f_cph <- cph(Surv(futime,fustat) ~ Age+Gender+Stage+riskScore,
             x=T, y=T, surv=T,
             data=rt)
f_cph

print(f_cph)
ddist <- datadist(rt)
options(datadist='ddist')
med  <- Quantile(f_cph)
surv <- Survival(f_cph) 

pdf("nomogram_distribution.pdf")
regplot(f_cph,        #对观测2的六个指标在列线图上进行计分展示        
        #        observation=rt[6,], #也可以不展示        
        points=TRUE,        
        plots=c("density","no plot"),        #预测1年和2年的死亡风险，此处单位是day        
        failtime = c(3,5,10),
        odds=F,        
        droplines=F,        
        leftlabel=T,        
        prfail = TRUE, #cox回归中需要TRUE        
        showP = T, #是否展示统计学差异        
        #droplines = F,#观测2示例计分是否画线        #    
        #        colors = mycol, #用前面自己定义的颜色        
        rank="range", #根据统计学差异的显著性进行变量的排序        
        interval="confidence",        
        title="Cox regression") #展示观测的可信区间## [[1]]##   
dev.off()



##############################################################################################
risk<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
result_time<-clinical
head(result_time)

risk$sample_name<-row.names(risk)
head(risk)
merged_data<-merge(risk, result_time, by.x = "sample_name", by.y = "Sample")
merged_data<-subset(merged_data, select = c("sample_name","futime.y","fustat.x","CYB561D2","LINC00638","DANCR","riskScore",
                                            "risk","Gender","Age","TNM_stage"))
colnames(merged_data)<-c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore","risk","Gender",
                         "Age","Stage")
merged_data$Gender_num<-ifelse(merged_data$Gender=="Female","2",
                               ifelse(merged_data$Gender=="Male","1",
                                      "NA"))
merged_data$Age_num<-ifelse(merged_data$Age>=60,"1",
                            ifelse(merged_data$Age<60,"2",
                                   "NA"))
merged_data$Stage_num<-merged_data$Stage
rownames(merged_data)<-merged_data$sample_name
head(merged_data)
merged_data_selected<-subset(merged_data,select = c("futime", "Status", "Age_num", "Gender_num", "Stage_num", "riskScore"))
merged_data_selected$futime<-as.numeric(merged_data_selected$futime)
merged_data_selected$Status<-as.numeric(merged_data_selected$Status)
merged_data_selected$Age_num<-as.numeric(merged_data_selected$Age_num)
merged_data_selected$Gender_num<-as.numeric(merged_data_selected$Gender_num)
merged_data_selected$Stage_num<-as.numeric(merged_data_selected$Stage_num)
merged_data_selected_rmNA<-na.omit(merged_data_selected)
f<-coxph(Surv(futime,Status)~Age_num+Gender_num+Stage_num+riskScore,merged_data_selected_rmNA)
merged_data_selected_rmNA$pr_failure3<-c(1-(summary(survfit(f,newdata=merged_data_selected_rmNA),times=3)$surv))
merged_data_selected_rmNA$pr_failure5<-c(1-(summary(survfit(f,newdata=merged_data_selected_rmNA),times=5)$surv))
merged_data_selected_rmNA$pr_failure10<-c(1-(summary(survfit(f,newdata=merged_data_selected_rmNA),times=10)$surv))
head(merged_data_selected_rmNA)

library(dplyr)
library(dcurves)
pdf("DCA_year3.pdf")
dca(Surv(futime,Status)~pr_failure3,
    data = merged_data_selected_rmNA,
    time = 3,
    thresholds = 1:50/100) %>%
  plot(smooth = T)
dev.off()

pdf("DCA_year5.pdf")
dca(Surv(futime,Status)~pr_failure5,
    data = merged_data_selected_rmNA,
    time = 5,
    thresholds = 1:50/100) %>%
  plot(smooth = T)
dev.off()

pdf("DCA_year10.pdf")
dca(Surv(futime,Status)~pr_failure10,
    data = merged_data_selected_rmNA,
    time = 10,
    thresholds = 1:50/100) %>%
  plot(smooth = T)
dev.off()

