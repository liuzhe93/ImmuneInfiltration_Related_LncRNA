setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
rm(list=ls())
#安装TCGAbiolinks
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("TCGAbiolinks")
#######################################一、数据下载阶段######################################
# 第一步：GDCquery（）筛选我们需要的数据，TCGAbiolinks包下载TCGA数据进行表达差异分析-肝癌案例
library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)
#选定要下载的cancer类型
TCGAbiolinks::getGDCprojects()$project_id
cancer_type="TCGA-COAD"
#选择下载你想要的数据类型
clinical<-GDCquery_clinic(project=cancer_type,type="clinical")
# 查看下载的数据
dim(clinical)
# 461  70
View(clinical)
save(clinical,file="COAD_clinical.Rdata")
write.csv(clinical, file="TCGAbiolinks-COAD-clinical.csv")
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
query <- GDCquery(project = cancer_type,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = "STAR - Counts")
samplesDown <- getResults(query,cols=c("cases"))  
#getResults(query, rows, cols)根据指定行名或列名从query中获取结果,此处用来获得样本的barcode
# 此处共检索出521个barcodes
# 从samplesDown中筛选出TP（实体肿瘤）样本的barcodes
# TCGAquery_SampleTypes(barcode, typesample)
# TP代表PRIMARY SOLID TUMOR；NT-代表Solid Tissue Normal（其他组织样本可参考学习文档）
##此处共检索出478个TP样本barcodes 41个NT样本barcode
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
write.csv(dataSmNT,"Normal_barcode.csv",quote=F,row.names=F)
write.csv(dataSmTP,"Cancer_barcode.csv",quote=F,row.names=F)
queryDown <- GDCquery(project = "TCGA-COAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "STAR - Counts", 
                      barcode = dataSmTP)
#barcode参数：根据传入barcodes进行数据过滤
######################################第二步：GDCdownload()下载GDCquery（）得到的结果#######################
# 下载数据，默认存放位置为当前工作目录下的GDCdata文件夹中。
setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset/download_data")
#GDCdownload(queryDown,method = "api", directory = "GDCdata",
#            files.per.chunk = 10)
GDCdownload(queryDown,method = "api", directory = "GDCdata",files.per.chunk = 20)
#method ；"API"或者"client"。"API"速度更快，但是容易下载中断。
#directory：下载文件的保存地址。Default: GDCdata。
#files.per.chunk = NULL:使用API下载大文件的时候，可以把文件分成几个小文件来下载，可以解决下载容易中断的问题。
GDCdownload(query = queryDown)
#读取下载的数据并将其准备到R对象中，在工作目录生成（save=TRUE）COAD_case.rda文件
# GDCprepare():Prepare GDC data,准备GDC数据，使其可用于R语言中进行分析
########第三步：GDCprepare()将前面GDCquery（）的结果准备成R语言可处理的SE（SummarizedExperiment）文件。######
dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename =
                          "COAD_case.rda")
########第四步：TCGAanalyze_Preprocessing()对数据进行预处理：使用spearman相关系数去除数据中的异常值##########
# 去除dataPrep1中的异常值，dataPrep1数据中含有肿瘤组织和正常组织的数据
testData<-load("COAD_case.rda")
testData
data
assays(data)
names(assays(data))[5]
dataPrep2 <- TCGAanalyze_Preprocessing(object = data,
                                       cor.cut = 0.6,
                                       datatype = "fpkm_unstrand")
#Number of outliers: 2
#将预处理后的数据dataPrep2，写入新文件“COAD_dataPrep.csv”
write.csv(dataPrep2,file = "COAD_dataPrep.csv",quote = FALSE)
#########################第五步：TCGAtumor_purity（）筛选肿瘤纯度大于60%的肿瘤barcodes########################
# TCGAtumor_purity(barcodes, estimate, absolute, lump, ihc, cpe)，使用来自5种方法的5个估计值作为阈值对TCGA样本进行过滤，这5个值是estimate, absolute, lump, ihc, cpe，这里设置cpe=0.6（cpe是派生的共识度量，是将所有方法的标准含量归一化后的均值纯度水平，以使它们具有相等的均值和标准差）
#筛选肿瘤纯度大于等于60%的样本数据
purityDATA <- TCGAtumor_purity(colnames(data), 0, 0, 0, 0, 0.6)
# filtered 为被过滤的数据， pure_barcodes是我们要的肿瘤数据
#the following TCGA barcodes do not have info on tumor purity:
#[1] "TCGA-A6-2672-01B-03R-2302-07"
Purity_COAD<-purityDATA$pure_barcodes
dataPrep2<-as.data.frame(dataPrep2)
overlap<-intersect(colnames(dataPrep2), Purity_COAD)
puried_data<-subset(dataPrep2, select = overlap)
###################################第七步：进行表达矩阵基因注释################################################
#基因注释,需要加载“SummarizedExperiment”包，“SummarizedExperiment container”每个由数字或其他模式的类似矩阵的对象表示。行通常表示感兴趣的基因组范围和列代表样品。
library("SummarizedExperiment")
anno_data<-rowData(data)   #传入数据dataPrep1必须为SummarizedExperiment对象
anno_data<-as.data.frame(anno_data)
anno_data<-subset(anno_data, select = c("gene_id", "gene_name"))
puried_data$gene_id<-rownames(puried_data)
merged_data<-merge(anno_data, puried_data, by = "gene_id")
write.csv(merged_data,file = "puried.COAD.csv",quote = FALSE, row.names = F)
merged_data_uniq<-aggregate(merged_data[,3:ncol(merged_data)], list(merged_data$gene_name), mean)
dim(merged_data)
#[1] 60660   453
rownames(merged_data_uniq)<-merged_data_uniq$Group.1
merged_data_uniq$Group.1<-NULL
dim(merged_data_uniq)
#[1] 59427   451
write.csv(merged_data_uniq,file = "puried.COAD.uniq.csv",quote = FALSE)


#########################第八步：进行表达矩阵标准化和过滤，得到用于差异分析的表达矩阵##########################
#library("EDASeq")
#dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(tabDF = merged_data_uniq,
#                                                    geneInfo = geneInfo,
#                                                    method = "gcContent")
#将标准化后的数据再过滤，去除掉表达量较低（count较低）的基因，得到最终的数据
#dataNorm[is.na(dataNorm)] <- 0
#dim(dataNorm)
#[1] 15901   451
#将标准化后的数据再过滤，去除掉表达量较低（count较低）的基因，得到最终的数据
#dataNorm<-na.omit(dataNorm)
#dim(dataNorm)
#[1] 13940   451
#dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                  method = "quantile", 
#                                  qnt.cut =  0.25)
#str(dataFilt)
#num [1:12925, 1:390] 0 30 803 0 0 2 0 263 6 4 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:12925] "A1BG" "A1CF" "A2M" "A4GALT" ...
#..$ : chr [1:390] "TCGA-ZP-A9CY-01A-11R-A38B-07" "TCGA-DD-AAVZ-01A-11R-A41C-07" "TCGA-DD-AADN-01A-11R-A41C-07" "TCGA-DD-A1EB-01A-11R-A131-07" ...
#write.csv(dataFilt,file = "TCGA_COAD_final.csv",quote = FALSE)  

#读入表达矩阵文件
dataFilt_COAD_final <- read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset/download_data/puried.COAD.uniq.csv",
                                header = T,check.names = FALSE,row.names = 1)
# 先看一下矩阵长啥样，心里有个数：每一行是一个基因，每一列是一个样本
#View(dataFilt_COAD_final)

data_sel<-rbind(dataFilt_COAD_final["CYB561D2",], dataFilt_COAD_final["LINC00638",], dataFilt_COAD_final["DANCR",])
data_sel_t<-t(data_sel)
data_sel_t<-as.data.frame(data_sel_t)
data_sel_t$id<-substring(rownames(data_sel_t),1,12)


setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/08_independent_TCGAdataset")
cli_data<-read.csv("TCGAbiolinks-COAD-clinical.csv",header=T,row.names=1)
cli_data$os<-ifelse(cli_data$vital_status=='Alive',cli_data$days_to_last_follow_up,cli_data$days_to_death)
library(dplyr)
cli_data_sel<-cli_data %>%
  select(submitter_id,os,vital_status,gender,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n,age_at_index)
#0=alive, 1=dead.
cli_data_sel$fustat<-ifelse(cli_data_sel$vital_status=='Alive',0,1)
cli_data_sel$futime<-cli_data_sel$os
results_input<-cli_data_sel %>%
  select(submitter_id,futime,fustat,age_at_index,gender,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n)
colnames(results_input)<-c("id","futime","fustat","age","gender","AJCC_stage","AJCC_T","AJCC_M","AJCC_N")
write.csv(results_input,"indepInput.csv",quote=F,row.names = F)

exp_cli<-merge(results_input,data_sel_t,by="id")
exp_cli<-subset(exp_cli, select = c("id", "futime", "fustat", "CYB561D2", "LINC00638", "DANCR"))

for(i in 1:nrow(exp_cli)){
  exp_cli$id[i]<-paste(exp_cli$id[i],i,sep = "_")
}
write.csv(exp_cli, "TCGA_exp_surv.csv", row.names = F, quote = F)


#cli_data_sel$Age_ope<-ifelse(cli_data_sel$age_at_index>=60,">=60","<60")

