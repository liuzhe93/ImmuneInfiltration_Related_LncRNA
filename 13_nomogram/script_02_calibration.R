setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/13_nomogram")
rm(list = ls())

library(rms)
library(timeROC)
risk<-read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",
                 sep="\t",header=T, row.names = 1)
risk$sample_name<-row.names(risk)
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
result_time<-clinical
head(result_time)


merged_data<-merge(risk, result_time, by.x = "sample_name", by.y = "Sample")
merged_data<-subset(merged_data, select = c("sample_name","futime.y","fustat.x","CYB561D2","LINC00638",
                                            "DANCR","riskScore","risk","Gender","Age","TNM_stage"))
colnames(merged_data)<-c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore","risk","Gender",
                         "Age","Stage")
merged_data$Gender_num<-ifelse(merged_data$Gender=="Female","2",
                               ifelse(merged_data$Gender=="Male","1",
                                      "NA"))
merged_data$Age_num<-ifelse(merged_data$Age>=60,"1",
                            ifelse(merged_data$Age<60,"2",
                                   "NA"))
merged_data$Stage_num<-merged_data$Stage
head(merged_data)
data_selected<-subset(merged_data, select = c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR",
                                              "riskScore","risk","Gender","Age","Stage"))
colnames(data_selected)<-c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore",
                           "risk","Gender","Age","Stage")
rt<-data_selected
head(rt)
rownames(rt)<-rt$sample_name
rt<-subset(rt, select = c("futime", "Status", "Age", "Gender", "Stage", "riskScore"))
colnames(rt)<-c("futime", "fustat", "Age", "Gender", "Stage", "riskScore")
rt$Stage<-ifelse(rt$Stage=="0","Stage 0",
                 ifelse(rt$Stage=="1","Stage I",
                        ifelse(rt$Stage=="2","Stage II",
                               ifelse(rt$Stage=="3","Stage III",
                                      ifelse(rt$Stage=="4","Stage IV",
                                             "NA")))))
head(rt)
ddist <- datadist(rt)
options(datadist='ddist')
f_cph <- cph(Surv(futime,fustat) ~ Age+Gender+Stage+riskScore,
             x=T, y=T, surv=T,
             data=rt)
#查看多因素Cox分析结果，最下方可见其对应的Coef、p值
print(f_cph)
ddist <- datadist(rt)
options(datadist='ddist')
med  <- Quantile(f_cph)
surv <- Survival(f_cph) 
cal_1<-calibrate(f_cph,u=3,cmethod='KM',m=15,B=200)# usually B=200 or 300
cal_2<-calibrate(f_cph,u=5,cmethod='KM',m=15,B=200)# usually B=200 or 300
cal_3<-calibrate(f_cph,u=10,cmethod='KM',m=15,B=200)# usually B=200 or 300
#par(mar=c(7,4,4,3),cex=1.0)
pdf("calibrate_3_years.pdf",width=6,height=6) 
plot(cal_1,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1 years OS',#便签
     ylab='Actual 3 years OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0)) #控制坐标轴的位置
dev.off()
pdf("calibrate_5_years.pdf",width=6,height=6) 
plot(cal_2,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5 years OS',#便签
     ylab='Actual 5 years OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0)) #控制坐标轴的位置
dev.off()
pdf("calibrate_10_years.pdf",width=6,height=6) 
plot(cal_3,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 10 years OS',#便签
     ylab='Actual 10 years OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0)) #控制坐标轴的位置
dev.off()


