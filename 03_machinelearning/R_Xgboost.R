setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/03_machinelearning/Xgboost")
rm(list=ls())

#ref: https://zhuanlan.zhihu.com/p/112009086

#################################3. Xgboost#####################################################
library('survival')
library('xgboost')
library('gbm')
# set random state
set.seed(1)
uniSigExp<-read.table("../uniSigExp.txt",sep="\t",header=T,row.names = 1)
uniSigExp$fustat<-as.numeric(uniSigExp$fustat)
x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
y=data.matrix(uniSigExp[,c(1:2)])
rownames(y)<-rownames(x)
mydata<-cbind(y,x)
mydata<-na.omit(mydata)
mydata<-as.data.frame(mydata)
colnames(mydata)[1:2]<-c("futime","fustat")
dim(mydata)
#[1] 562  73

index<-sort(sample(nrow(mydata),nrow(mydata)*.9))
train<-mydata[index,]
test<-mydata[-index,]

x_train=as.matrix(train[,c(3:ncol(train))])
y_train=data.matrix(train$fustat)
rownames(y_train)<-rownames(x_train)
y_train<-as.data.frame(y_train)

x_test=as.matrix(test[,c(3:ncol(test))])
y_test=data.matrix(test$fustat)
rownames(y_test)<-rownames(x_test)
y_test<-as.data.frame(y_test)


train_fin <- list(data=x_train,label=y_train$V1)
test_fin <- list(data=x_test,label=y_test$V1)
dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label) 
dtest <- xgb.DMatrix(data = test_fin$data, label = test_fin$label)
xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5,  
               objective='binary:logistic', nround=25)

importance <- xgb.importance(colnames(x), model = xgb)  
head(importance)
pdf("Imporece.pdf")
xgb.ggplot.importance(importance)
dev.off()

#混淆矩阵
pre_xgb = round(predict(xgb,newdata = dtest))
table(y_test$V1,pre_xgb,dnn=c("true","pre"))
#    pre
#true  0  1
#   0 24 10
#   1 15 8
#ROC曲线
library(ROCR) 
library(pROC)
xgboost_roc <- roc(y_test$V1,as.numeric(pre_xgb))
pdf("ROCC.pdf")
plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2),grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,main='ROC curve')
dev.off()

impo_score<-importance
impo_score_fi<-impo_score[impo_score$Importance>0.01,]

genelist_xgboost<-impo_score_fi$Feature
length(genelist_xgboost)
#[1] 46
genelist_xgboost<-gsub("\\.","-",genelist_xgboost)
write.csv(genelist_xgboost, "genelist_Xgboost.csv", quote = F, row.names = F)



