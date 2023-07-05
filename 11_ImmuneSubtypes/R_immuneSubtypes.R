setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/11_ImmuneSubtypes")
rm(list=ls())

data<-read.csv("data.csv")
data<-subset(data,select = c("TCGA.Participant.Barcode","Immune.Subtype"))
data<-na.omit(data)
dataFilt_COAD_final <- read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset/download_data/puried.COAD.uniq.csv",
                                row.names = 1)
dataFilt_COAD_final[1:5,1:5]
head(data)
dataFilt_COAD_final_t<-t(dataFilt_COAD_final)
dataFilt_COAD_final_t<-as.data.frame(dataFilt_COAD_final_t)
dataFilt_COAD_final_t$Sample<-rownames(dataFilt_COAD_final_t)
exp_CYB561D2<-dataFilt_COAD_final_t[,"CYB561D2"]
exp_LINC00638<-dataFilt_COAD_final_t[,"LINC00638"]
exp_DANCR<-dataFilt_COAD_final_t[,"DANCR"]
exp<-cbind(dataFilt_COAD_final_t$Sample, exp_CYB561D2, exp_LINC00638, exp_DANCR)
colnames(exp)<-c("SampleName","CYB561D2","LINC00638","DANCR")
exp<-as.data.frame(exp)
exp$SampleName<-substr(exp$SampleName,1,12)
exp$SampleName<-gsub("\\.","-",exp$SampleName)
data1<-merge(data,exp,by.x="TCGA.Participant.Barcode",by.y="SampleName")

library(survival)                                         #引用包


rt=read.table("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/06_multivariateCox_regression/risk.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取输入文件

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

result_time<-read.csv("result_time.csv")

input_data<-merge(result_time,data1,by.x = "Id", by.y = "TCGA.Participant.Barcode")
input_data<-na.omit(input_data)
input_data_C1<-subset(input_data,Immune.Subtype=="C1")
input_data_C1$futime<-as.numeric(input_data_C1$futime)
input_data_C1$fustat<-as.numeric(input_data_C1$fustat)
input_data_C1$CYB561D2<-as.numeric(input_data_C1$CYB561D2)
input_data_C1$LINC00638<-as.numeric(input_data_C1$LINC00638)
input_data_C1$DANCR<-as.numeric(input_data_C1$DANCR)

#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=input_data_C1)
input_data_C1$Score<-riskScore

input_data_C2<-subset(input_data,Immune.Subtype=="C2")
input_data_C2$futime<-as.numeric(input_data_C2$futime)
input_data_C2$fustat<-as.numeric(input_data_C2$fustat)
input_data_C2$CYB561D2<-as.numeric(input_data_C2$CYB561D2)
input_data_C2$LINC00638<-as.numeric(input_data_C2$LINC00638)
input_data_C2$DANCR<-as.numeric(input_data_C2$DANCR)
#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=input_data_C2)
input_data_C2$Score<-riskScore

input_data_C3<-subset(input_data,Immune.Subtype=="C3")
input_data_C3$futime<-as.numeric(input_data_C3$futime)
input_data_C3$fustat<-as.numeric(input_data_C3$fustat)
input_data_C3$CYB561D2<-as.numeric(input_data_C3$CYB561D2)
input_data_C3$LINC00638<-as.numeric(input_data_C3$LINC00638)
input_data_C3$DANCR<-as.numeric(input_data_C3$DANCR)
#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=input_data_C3)
input_data_C3$Score<-riskScore

input_data_C4<-subset(input_data,Immune.Subtype=="C4")
input_data_C4$futime<-as.numeric(input_data_C4$futime)
input_data_C4$fustat<-as.numeric(input_data_C4$fustat)
input_data_C4$CYB561D2<-as.numeric(input_data_C4$CYB561D2)
input_data_C4$LINC00638<-as.numeric(input_data_C4$LINC00638)
input_data_C4$DANCR<-as.numeric(input_data_C4$DANCR)
#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=input_data_C4)
input_data_C4$Score<-riskScore

input_data_C5<-subset(input_data,Immune.Subtype=="C5")
input_data_C5$futime<-as.numeric(input_data_C5$futime)
input_data_C5$fustat<-as.numeric(input_data_C5$fustat)
input_data_C5$CYB561D2<-as.numeric(input_data_C5$CYB561D2)
input_data_C5$LINC00638<-as.numeric(input_data_C5$LINC00638)
input_data_C5$DANCR<-as.numeric(input_data_C5$DANCR)
#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=input_data_C5)
input_data_C5$Score<-riskScore

input_data_C6<-subset(input_data,Immune.Subtype=="C6")
input_data_C6$futime<-as.numeric(input_data_C6$futime)
input_data_C6$fustat<-as.numeric(input_data_C6$fustat)
input_data_C6$CYB561D2<-as.numeric(input_data_C6$CYB561D2)
input_data_C6$LINC00638<-as.numeric(input_data_C6$LINC00638)
input_data_C6$DANCR<-as.numeric(input_data_C6$DANCR)
#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=input_data_C6)
input_data_C6$Score<-riskScore

merged_risk<-rbind(input_data_C1,input_data_C2,input_data_C3,input_data_C4,input_data_C6)
dim(input_data_C1)
#[1] 331  8
dim(input_data_C2)
#[1] 82 8
dim(input_data_C3)
#[1]  7 8
dim(input_data_C4)
#[1] 12 8
dim(input_data_C6)
#[1]  1 8
table(input_data$Immune.Subtype)
# C1  C2  C3  C4  C6 
#331  82   7  12   1
merged_risk$CellType<-c(rep("C1",331),rep("C2",82),rep("C3",7),rep("C4",12),rep("C6",1))
merged_risk_filtered<-subset(merged_risk,Score < 1)
library(ggpubr)
pdf("RiskScore_immune.pdf",width=4,height = 4)
ggboxplot(merged_risk_filtered, x="CellType", y="Score",
          color="CellType", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
#  geom_hline(yintercept = mean(merged_risk_filtered$Score),
#             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "anova", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")+                      # Pairwise comparison against all
  ylim(0,1)
dev.off()

pdf("CYB561D2_immune.pdf",width=4,height = 4)
ggboxplot(merged_risk_filtered, x="CellType", y="CYB561D2",
          color="CellType", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
  #  geom_hline(yintercept = mean(merged_risk_filtered$Score),
  #             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "anova", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                     # Pairwise comparison against all
dev.off()


pdf("LINC00638_immune.pdf",width=4,height = 4)
ggboxplot(merged_risk_filtered, x="CellType", y="LINC00638",
          color="CellType", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
  #  geom_hline(yintercept = mean(merged_risk_filtered$Score),
  #             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "anova", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                     # Pairwise comparison against all
dev.off()

pdf("DANCR_immune.pdf",width=4,height = 4)
ggboxplot(merged_risk_filtered, x="CellType", y="DANCR",
          color="CellType", add="jitter",
          legend="none") + 
  rotate_x_text(angle = 45) + 
  #  geom_hline(yintercept = mean(merged_risk_filtered$Score),
  #             linetype=2) + # 添加base mean的水平线
  stat_compare_means(method = "anova", label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                     # Pairwise comparison against all
dev.off()
