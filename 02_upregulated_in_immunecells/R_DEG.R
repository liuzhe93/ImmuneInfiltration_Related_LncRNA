setwd("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells")
rm(list=ls())

### colon cancer samples
library("GEOquery")
library("affy")
library(samr)

fdrcutoff<-0.05
normal_labels<-read.csv("immune_label.csv")
normal_labels<-normal_labels[,-1]
dim(normal_labels)
#[1] 115   2
cancer_labels<-read.csv("coloncancer_label.csv")
cancer_labels<-cancer_labels[,-1]
dim(cancer_labels)
#[1] 566   2
# 585-566=19 non-tumor samples

library("GEOquery")
library("affy")
#GSEName="GSE3910"
#rawdata = getGEOSuppFiles(GSEName)
#dir <- choose.dir(caption = "Select folder")

##############################cancer data analysis##################################################
cel.files <- list.files(path = "/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/02_upregulated_in_immunecells/coloncancer_selected", 
                        pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
data.raw <- ReadAffy(filenames = cel.files)
sampleNames(data.raw)
est <- rma(data.raw)
expmtx_expso <- exprs(est)
write.csv(expmtx_expso, file="expmtx_expso.csv", quote = F)
library("hgu133plus2.db")
ids<-toTable(hgu133plus2SYMBOL)
ids
length(unique(ids$symbol))
#[1] 20857 the number of probe in the annotation package
library("stringr")
table(rownames(expmtx_expso) %in% ids$probe_id)
#FALSE  TRUE 
#11537 43138
data<-expmtx_expso[rownames(expmtx_expso) %in% ids$probe_id,]
dim(data)
#[1] 43138   566
ids<-ids[match(rownames(data),ids$probe_id),]
tmp<-by(data, ids$symbol, function(x) rownames(x)[which.max(rowMeans(x))])
probes<-as.character(tmp)
data<-data[rownames(data) %in% probes,]
ids<-ids[ids$probe_id %in% probes,]
row.names(data)<-ids$symbol
write.csv(data, "cancer_exp.csv", quote = F)
cancer_exp<-data


##############################normal data analysis##################################################
cel.files <- list.files(path = "/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA/immunecells_selected", 
                        pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)
basename(cel.files)
data.raw <- ReadAffy(filenames = cel.files)
sampleNames(data.raw)
est <- rma(data.raw)
expmtx_expso <- exprs(est)
write.csv(expmtx_expso, file="expmtx_expso_immune.csv", quote = F)
ids<-toTable(hgu133plus2SYMBOL)
ids
length(unique(ids$symbol))
#[1] 20857 the number of probe in the annotation package
library("stringr")
table(rownames(expmtx_expso) %in% ids$probe_id)
#FALSE  TRUE 
#11537 43138
data<-expmtx_expso[rownames(expmtx_expso) %in% ids$probe_id,]
dim(data)
#[1] 43138   115
ids<-ids[match(rownames(data),ids$probe_id),]
tmp<-by(data, ids$symbol, function(x) rownames(x)[which.max(rowMeans(x))])
probes<-as.character(tmp)
data<-data[rownames(data) %in% probes,]
ids<-ids[ids$probe_id %in% probes,]
row.names(data)<-ids$symbol
write.csv(data, "immune_exp.csv", quote = F)
immune_exp<-data


cancer_exp<-as.data.frame(cancer_exp)
cancer_exp$gene_name<-rownames(cancer_exp)

immune_exp<-as.data.frame(immune_exp)
immune_exp$gene_name<-rownames(immune_exp)

hklncRNA<-read.csv("/Users/liuzhe/Desktop/cityu/LncRNA_CRC/analysis/01_candidate_lncRNA/hklncRNA_inf.csv")
genelist<-unique(hklncRNA$X)

cancer_exp_sel<-cancer_exp[genelist,]
immune_exp_sel<-immune_exp[genelist,]

normal<-immune_exp_sel[,c("gene_name",normal_labels$colnames.expr.rma.)]
cancer<-cancer_exp_sel[,c("gene_name",cancer_labels$colnames.expr.rma.)]
normal<-normal[,-1]
cancer<-cancer[,-1]

normal<-as.matrix(normal)
cancer<-as.matrix(cancer)
gid<-as.matrix(genelist)

label<-c(rep(1,115),rep(2,566));
exp<-cbind(normal,cancer)
exp_fil<-na.omit(exp)
samfit<-SAM(exp_fil,label,resp.type="Two class unpaired",geneid=gid,genenames=gid,nperms=1000,logged2=T,fdr.output=fdrcutoff)

siggene.up<-cbind(samfit$siggenes.table$genes.up[,1],as.numeric(samfit$siggenes.table$genes.up[,7])/100)
siggene.down<-cbind(samfit$siggenes.table$genes.lo[,1],as.numeric(samfit$siggenes.table$genes.lo[,7])/100)
sig.up<-siggene.up[siggene.up[,2]<fdrcutoff,1];
sig.down<-siggene.down[siggene.down[,2]<fdrcutoff,1];
sig.all<-c(sig.up,sig.down)
write.table(sig.all,"all_005.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(sig.up,"up_005.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(sig.down,"down_005.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

#We focus on genes that up regulated in immune cells compared to cancer cells
res<-samfit$siggenes.table$genes.lo
res<-as.data.frame(res)
dim(res)
#[1] 87   7
res$`Gene ID`<-NULL
write.csv(res, "Upreg_in_immune.csv", quote = F, row.names = F)



############################################heatmap###########################################
##lncRNAs that up-regulated in immune cells compared to cancer cells
library("pheatmap")
cancer_exp<-read.csv("cancer_exp.csv")
immune_exp<-read.csv("immune_exp.csv")
up<-read.table("up_005.txt", header = T)
down<-read.table("down_005.txt", header = T)
downgenes<-up$x
upgenes<-down$x
merged_data<-merge(immune_exp, cancer_exp, by="X")
rownames(merged_data)<-merged_data$X
merged_data<-merged_data[,-1]
exp_up<-merged_data[upgenes,]
exp_down<-merged_data[downgenes,]
exp_up_fil<-na.omit(exp_up)
exp_down_fil<-na.omit(exp_down)
exp_data<-rbind(exp_up_fil, exp_down_fil)

annotation_col = data.frame(Type=c(rep("immune_cell", 115), rep("cancer_cell", 566)))
annotation_col$Type<-as.factor(annotation_col$Type)
rownames(annotation_col) = colnames(exp_data)
rownam<-c(rep("Up_regulated_lncRNAs",71),rep("down_regulated_lncRNAs",29))
annotation_row = matrix(rownam,nrow = 100,ncol=1,byrow=TRUE)
rownames(annotation_row) = rownames(exp_data)
colnames(annotation_row)<-"gene_type"
annotation_row<-as.data.frame(annotation_row)
ann_colors = list(Type = c(immune_cell = "#ed1299", cancer_cell = "#09f9f5"),
                  gene_type = c(Up_regulated_lncRNAs = "red", down_regulated_lncRNAs  = "blue"))
pdf(file="heatmap_DEGs.pdf")
pheatmap(exp_data, scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F,
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col,
         annotation_row = annotation_row, gaps_row = 71, gaps_col = 115,
         angle_col = "45", annotation_colors = ann_colors,  
         show_rownames = F, show_colnames = F, main = "Title", fontsize = 6)
dev.off()




