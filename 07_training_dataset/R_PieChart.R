setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/14_data_distribution")
rm(list=ls())

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
                                            "risk","Gender","Age","TNM_stage"))
colnames(merged_data)<-c("sample_name","futime","Status","CYB561D2","LINC00638","DANCR","riskScore","risk","Gender",
                         "Age","Stage")
library("ggplot2")
library("ggforce")
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
###############################################Status#################################################################
########high risk group##############
merged_data_high_status<-merged_data_high[,"Status"]
table(merged_data_high$Status)
#  0   1 
#169 112 
num<-c(169,112)
type<-c("Alive","Dead")
status_high<-data.frame(num,type)
pdf("GEO_high_status.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = status_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_status<-merged_data_low[,"Status"]
table(merged_data_low$Status)
#  0   1 
#202  79
num<-c(202,79)
type<-c("Alive","Dead")
status_low<-data.frame(num,type)
pdf("GEO_low_status.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = status_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_status_data<-cbind(status_high,status_low)
row.names(merged_status_data)<-merged_status_data$type
merged_status_data<-merged_status_data[,-2]
merged_status_data<-merged_status_data[,-3]
colnames(merged_status_data)<-c("high", "low")
chisq.test(merged_status_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_status_data
#X-squared = 8.1214, df = 1, p-value = 0.004375

###############################################Stage#################################################################
merged_data_cp<-merged_data
merged_data$Stage_ope<-ifelse(merged_data$Stage=="0", "Stage 0",
                              ifelse(merged_data$Stage=="1", "Stage I",
                                     ifelse(merged_data$Stage=="2", "Stage II",
                                            ifelse(merged_data$Stage=="3", "Stage III",
                                                   ifelse(merged_data$Stage=="4", "Stage IV",
                                                          "NA")))))
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
########high risk group##############
merged_data_high_stage<-merged_data_high[,"Stage_ope"]
table(merged_data_high$Stage_ope)
#Stage 0   Stage I  Stage II Stage III  Stage IV 
#      2        13       115       107        44 
num<-c(2,13,115,107,44)
type<-c("Stage 0","Stage I","Stage II","Stage III","Stage IV")
stage_high<-data.frame(num,type)
pdf("GEO_high_stage.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("orange","red","green","blue","black"))+
  geom_arc_bar(data = stage_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_stage<-merged_data_low[,"Stage_ope"]
table(merged_data_low$Stage_ope)
#Stage 0   Stage I  Stage II Stage III  Stage IV 
#      2        19       147        97        16 
num<-c(2,19,147,97,16)
type<-c("Stage 0","Stage I","Stage II","Stage III","Stage IV")
stage_low<-data.frame(num,type)
pdf("GEO_low_stage.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("orange","red","green","blue","black"))+
  geom_arc_bar(data = stage_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_stage_data<-cbind(stage_high,stage_low)
row.names(merged_stage_data)<-merged_stage_data$type
merged_stage_data<-merged_stage_data[,-2]
merged_stage_data<-merged_stage_data[,-3]
colnames(merged_stage_data)<-c("high", "low")
chisq.test(merged_stage_data)
#Pearson's Chi-squared test
#
#data:  merged_stage_data
#X-squared = 18.59, df = 4, p-value = 0.0009458

###############################################Gender#################################################################
########high risk group##############
merged_data_high_Gender<-merged_data_high[,"Gender"]
table(merged_data_high$Gender)
#Female   Male 
#  118    163
num<-c(118,163)
type<-c("Female","Male")
gender_high<-data.frame(num,type)
pdf("GEO_high_gender.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = gender_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_Gender<-merged_data_low[,"Gender"]
table(merged_data_low$Gender)
#Female   Male 
#   135    146
num<-c(135,146)
type<-c("Female","Male")
gender_low<-data.frame(num,type)
pdf("GEO_low_gender.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = gender_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_gender_data<-cbind(gender_high,gender_low)
row.names(merged_gender_data)<-merged_gender_data$type
merged_gender_data<-merged_gender_data[,-2]
merged_gender_data<-merged_gender_data[,-3]
colnames(merged_gender_data)<-c("high", "low")
chisq.test(merged_gender_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_gender_data
#X-squared = 1.8403, df = 1, p-value = 0.1749

###############################################TCGA_Age#################################################################
merged_data$Age_ope<-ifelse(merged_data$Age>=60,">=60","<60")
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
########high risk group##############
merged_data_high_age<-merged_data_high[,"Age_ope"]
table(merged_data_high$Age_ope)
#<60 >=60 
#73  186 
num<-c(73,186)
type<-c("<60",">=60")
age_high<-data.frame(num,type)
pdf("TCGA_high_age.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = age_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_age<-merged_data_low[,"Age_ope"]
table(merged_data_low$Age_ope)
#<60 >=60 
#63  197
num<-c(63,197)
type<-c("<60",">=60")
age_low<-data.frame(num,type)
pdf("TCGA_low_age.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = age_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_age_data<-cbind(age_high,age_low)
row.names(merged_age_data)<-merged_age_data$type
merged_age_data<-merged_age_data[,-2]
merged_age_data<-merged_age_data[,-3]
colnames(merged_age_data)<-c("high", "low")
chisq.test(merged_age_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_age_data
#X-squared = 0.85476, df = 1, p-value = 0.3552

###############################################Age#################################################################
merged_data$Age_ope<-ifelse(merged_data$Age>=60,">=60","<60")
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
########high risk group##############
merged_data_high_age<-merged_data_high[,"Age_ope"]
table(merged_data_high$Age_ope)
#<60 >=60 
# 75  206 
num<-c(75,206)
type<-c("<60",">=60")
age_high<-data.frame(num,type)
pdf("GEO_high_age.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = age_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_age<-merged_data_low[,"Age_ope"]
table(merged_data_low$Age_ope)
#<60 >=60 
# 75  206
num<-c(75,206)
type<-c("<60",">=60")
age_low<-data.frame(num,type)
pdf("GEO_low_age.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = age_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_age_data<-cbind(age_high,age_low)
row.names(merged_age_data)<-merged_age_data$type
merged_age_data<-merged_age_data[,-2]
merged_age_data<-merged_age_data[,-3]
colnames(merged_age_data)<-c("high", "low")
chisq.test(merged_age_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_age_data
#X-squared = 0, df = 1, p-value = 1



